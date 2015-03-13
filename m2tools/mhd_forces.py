import checkpoint
import command
from autolog import logmethod


class HydromagneticForceCalculator(object):
    @logmethod
    def __init__(self, filename):
        self.checkpoint = checkpoint.Checkpoint(filename)

    @logmethod
    def magnetic_field(self):
        import numpy as np
        B = np.array([self.checkpoint.cell_primitive['B1'][:],
                      self.checkpoint.cell_primitive['B2'][:],
                      self.checkpoint.cell_primitive['B3'][:]])
        return B

    @logmethod
    def velocity_field(self):
        import numpy as np
        v = np.array([self.checkpoint.cell_primitive['v1'][:],
                      self.checkpoint.cell_primitive['v2'][:],
                      self.checkpoint.cell_primitive['v3'][:]])
        return v

    @logmethod
    def vector_potential(self):
        import numpy as np
        Bx = self.magnetic_field()
        Bk = np.fft.fftn(Bx, axes=[1,2,3])
        Ks = self.get_ks(Bx.shape)
        K2 = Ks[0]**2 + Ks[1]**2 + Ks[2]**2
        K2[0,0,0] = 1
        Ak = np.zeros_like(Bk)
        Ak[0] = 1.j*(Ks[1]*Bk[2] - Ks[2]*Bk[1])/K2
        Ak[1] = 1.j*(Ks[2]*Bk[0] - Ks[0]*Bk[2])/K2
        Ak[2] = 1.j*(Ks[0]*Bk[1] - Ks[1]*Bk[0])/K2
        Ax = np.fft.ifftn(Ak, axes=[1,2,3]).real
        return Ax

    @logmethod
    def helicity(self):
        import numpy as np
        A = self.vector_potential()
        B = self.magnetic_field()
        return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]

    def get_ks(self, shape):
        """
        Get the wavenumber array, assuming all dimensions are length 1, while
        the number of cells along each axis may be unique
        """
        import numpy as np
        N0, N1, N2 = shape[1:]
        Ks = np.zeros(shape)
        Ks[0] = np.fft.fftfreq(N0)[:,None,None]*2*N0*np.pi
        Ks[1] = np.fft.fftfreq(N1)[None,:,None]*2*N1*np.pi
        Ks[2] = np.fft.fftfreq(N2)[None,None,:]*2*N2*np.pi
        return Ks

    @logmethod
    def magnetic_gradient_tensor(self, method='spectral', return_B=False):
        """
        Form the gradient tensor gradB[i,j] := dB_j/dx_i using either spectral
        or finite-difference method
        """
        import numpy as np
        Bx = self.magnetic_field()
        Bk = np.fft.fftn(Bx, axes=[1,2,3])
        Ks = self.get_ks(Bx.shape)
        N = Ks.shape[1:]
        d = five_point_deriv
        gradB = np.zeros((3,3) + N)
        for i in range(3):
            for j in range(3):
                if method == 'spectral':
                    gradB[i,j] = np.fft.ifftn(1.j * Ks[i] * Bk[j]).real
                elif method == 'finite-difference':
                    gradB[i,j] = d(Bx[j], i, h=1.0/N[i])
                else:
                    raise ValueError("bad argument value method='%s'" % method)
        if return_B:
            return gradB, Bx
        else:
            return gradB

    @logmethod
    def magnetic_pressure_gradient_force(self, method='spectral'):
        """ contract B_j with dB_j/dx_i; return -grad(B^2)/2 """
        import numpy as np
        gradB, B = self.magnetic_gradient_tensor(method=method, return_B=True)
        F = np.zeros_like(B)
        for i in range(3):
            for j in range(3):
                F[i] -= B[j] * gradB[i,j]
        return F

    @logmethod
    def magnetic_tension(self, method='spectral'):
        """ contract B_i with dB_j/dx_i; return B dot grad(B) """
        import numpy as np
        gradB, B = self.magnetic_gradient_tensor(method=method, return_B=True)
        F = np.zeros_like(B)
        for i in range(3):
            for j in range(3):
                F[j] += B[i] * gradB[i,j]
        return F

    @logmethod
    def curlB(self):
        return curl_of_field(self.magnetic_field())

    @logmethod
    def power_spectrum(self, which="magnetic", bins=256):
        import numpy as np
        if which == "magnetic":
            fx = self.magnetic_field()
        elif which == "velocity":
            fx = self.velocity_field()
        else:
            raise ValueError()

        Ks = self.get_ks(fx.shape)
        Bk = np.fft.fftn(fx, axes=[1,2,3])
        B2 = (Bk[0]*Bk[0].conj() + Bk[1]*Bk[1].conj() + Bk[2]*Bk[2].conj()).real
        K = (Ks[0]**2 + Ks[1]**2 + Ks[2]**2)**0.5
        bins = np.linspace(K[K!=0].min(), K.max(), bins)
        dk, bins = np.histogram(K, bins=bins)
        dP, bins = np.histogram(K, bins=bins, weights=B2, normed=True)
        dk[dk==0] = 1
        x = 0.5*(bins[1:] + bins[:-1])
        y = dP / dk
        return x, y


def curl_of_field(B):
    import numpy as np
    d = lambda f, a: five_point_deriv(f, a, h=1.0/B.shape[a+1])
    curl = np.zeros_like(B)
    curl[0] = d(B[2], 1) - d(B[1], 2)
    curl[1] = d(B[0], 2) - d(B[2], 0)
    curl[2] = d(B[1], 0) - d(B[0], 1)
    return curl


def five_point_deriv(f, axis=0, h=1.0):
    import numpy as np
    return (-1*np.roll(f, -2, axis) +
            +8*np.roll(f, -1, axis) + 
            -8*np.roll(f, +1, axis) + 
            +1*np.roll(f, +2, axis)) / (12.0 * h)



class CalculateForces(command.Command):
    _alias = "calc-forces"
    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+')
        parser.add_argument("--fields", default='A')
        parser.add_argument("--diff")
        parser.add_argument("--noplot", action="store_true")

    def run(self, args):
        import matplotlib.pyplot as plt
        self.helicity = [ ]
        for filename in args.filenames:
            self.run_file(filename, args)
        if not args.noplot:
            plt.show()
        else:
            import pickle
            outf = open("helicity.pk", 'w')
            pickle.dump(self.helicity, outf)

    def run_file(self, filename, args):
        import matplotlib.pyplot as plt
        calc = HydromagneticForceCalculator(filename)
        diff = HydromagneticForceCalculator(args.diff) if args.diff else None

        def fields(chk, key):
            F = {
                'B' : lambda : chk.magnetic_field(),
                'J' : lambda : chk.curlB(),
                'A' : lambda : chk.vector_potential(),
                'H' : lambda : chk.helicity(),
                'curlA' : lambda : curl_of_field(chk.vector_potential())}
            return F[key]()

        if args.noplot: return

        for f in args.fields.split(','):

            j = calc.checkpoint.domain_resolution[1]/2
            data = fields(calc, f)[0][:,j,:] if f is not 'H' else fields(calc, f)[:,j,:]
            if diff is not None:
                data -= fields(calc, f)[0][:,j,:] if f is not 'H' else fields(calc, f)[:,j,:]

            if f is 'H':
                print "total helicity is", data.mean()

            plt.figure()
            plt.imshow(data, interpolation='nearest')
            plt.title(f)
            plt.colorbar()



class CalculateForceFreePotentialChange(command.Command):
    _alias = "ff-decay"
    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+')
        parser.add_argument("--diff")

    def run(self, args):
        import numpy as np
        import matplotlib.pyplot as plt
        
        alpha = 3**0.5 * 2 * np.pi

        f0 = HydromagneticForceCalculator(args.filenames[0])
        A0 = f0.vector_potential()
        B0 = f0.magnetic_field()
        J0 = f0.curlB()
        ff_violation0 = J0 - alpha * B0

        for filename in args.filenames[1:]:
            f = HydromagneticForceCalculator(filename)
            A = f.vector_potential()
            B = f.magnetic_field()
            J = f.curlB()
            ff_violation = J - alpha * B

            #plt.hist(ff_violation.flat, histtype='step', bins=256, log=True, label=filename)

            plt.figure()
            plt.suptitle("J - J0")
            plt.imshow((J - J0)[2][:,:,0])
            plt.colorbar()


            plt.figure()
            plt.suptitle("B - B0")
            plt.imshow((B - B0)[2][:,:,0])
            plt.colorbar()

            plt.figure()
            plt.suptitle("v3")
            plt.imshow(f.checkpoint.cell_primitive["v3"][:,:,0])
            plt.colorbar()


            plt.figure()
            plt.suptitle("B")
            plt.imshow(B[2][:,:,0])
            plt.colorbar()

        #plt.legend()
        plt.show()

        # dA = A1 - A0
        # dB = B1 - B0

        # quantile = np.percentile(ff_violation0, 99.9)
        # indices = np.where(ff_violation0 > quantile)
        # i0 = indices[1][0]
        # print i0


        # plt.figure()
        # plt.imshow(ff_violation0[1][i0,:,:])
        # plt.colorbar()

        # plt.figure()
        # plt.imshow(ff_violation0[1][0,:,:])
        # plt.colorbar()


        # d2M = self.dot(dB, dB - alpha * dA)
        # print d2M.mean() / d2M.std()


    def dot(self, A, B):
        return A[0]*B[0] + A[1]*B[1] + A[2]*B[2]

    def mod(self, A):
        return self.dot(A, A)**0.5



class CalculatePowerSpectrum(command.Command):
    _alias = "power-spectrum"
    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+', help="chkpt.h5:[velocity|magnetic]")
        parser.add_argument("--bins", type=int, default=256)

    def run(self, args):
        import matplotlib.pyplot as plt
        for filename in args.filenames:
            if ':' in filename:
                filename, args.which = filename.split(':')
            else:
                args.which = 'magnetic'
            x, y = self.run_file(filename, args)
        plt.loglog(x, y[0]*(x/x[0])**(-5./3), label=r'$k^{-5/3}$')
        plt.legend(loc='best')
        plt.show()

    def run_file(self, filename, args):
        import matplotlib.pyplot as plt
        calc = HydromagneticForceCalculator(filename)
        x, y = calc.power_spectrum(which=args.which, bins=args.bins)
        plt.xlim(x[0], x[-1])
        plt.loglog(x, y, label=filename+':'+args.which)
        return x, y


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    n = 5
    L = 3.0
    N = 256
    k = np.fft.fftfreq(N)*2*N*np.pi/L
    x = np.linspace(0.0, L*(1.0 - 1.0/N), N)
    y   = np.sin(2 * n * np.pi * x / L)
    yp1 = np.cos(2 * n * np.pi * x / L) * (2 * n * np.pi / L)
    yp2 = five_point_deriv(y, h=L/N)
    yp3 = np.fft.ifft(1.j * k * np.fft.fft(y)).real
    plt.plot(x, yp1, label="$y_1'(x)$")
    plt.plot(x, yp2, label="$y_2'(x)$")
    plt.plot(x, yp3, label="$y_3'(x)$")
    plt.legend()
    plt.show()
