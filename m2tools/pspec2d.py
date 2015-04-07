
from m2tools import checkpoint
from m2tools import command
from smooth import smooth


def power_spectrum(Bx, bins):
    import numpy as np
    N0, N1 = Bx.shape[1:]
    Ks = np.zeros(Bx.shape)
    Ks[0] = np.fft.fftfreq(N0)[:,None]*N0
    Ks[1] = np.fft.fftfreq(N1)[None,:]*N1
    Bk = np.fft.fftn(Bx, axes=[1,2])
    B2 = (Bk[0]*Bk[0].conj() + Bk[1]*Bk[1].conj() + Bk[2]*Bk[2].conj()).real
    K = (Ks[0]**2 + Ks[1]**2 + Ks[2]**2)**0.5
    dP, bins = np.histogram(K, bins=bins, weights=B2)
    dk = 1.0*(bins[1:] - bins[:-1])
    x  = 0.5*(bins[1:] + bins[:-1])
    y  = dP / dk
    return x, y



class PowerSpectrum2D(command.Command):
    _alias = "pspec2d"
    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+')
        parser.add_argument("--time-series", action="store_true", default=False)
        parser.add_argument("--nocache", action="store_true", default=False)
        parser.add_argument("--recalc", action="store_true", default=False)
        parser.add_argument("--noshow", action="store_true", default=False,
                            help="don't plot anything, just cache the answer")
        parser.add_argument("-s", "--smooth", default=0, type=int)
        parser.add_argument("-c", "--compensate", default=0, type=float)
        parser.add_argument("--bins", type=str, default='1:2048:512',
                            help="if calculating, k0:k1:num_bin_edges; for time-series k0,k1,k2,...")
        parser.add_argument("--title", default=None, type=str)
        parser.add_argument("--fit-window", default=None, help="k0:k1 or t0:t1")
        parser.add_argument("--save-series", default=None, type=str,
                            help="name of a file to write the time series")


    def run(self, args):
        import matplotlib.pyplot as plt
        import numpy as np

        plt.figure(figsize=[10,8])

        if args.time_series:
            self.do_time_series(args)
        else:
            self.do_Pofk(args)
    
        if not args.noshow:
            plt.suptitle(args.title)
            plt.show()


    def make_bins(self, args):
        import numpy as np
        k0, k1, num_bin_edges = args.bins.split(':')
        bins = np.linspace(float(k0), float(k1), int(num_bin_edges))
        return bins


    def do_Pofk(self, args):
        import matplotlib.pyplot as plt
        import numpy as np

        for filename in args.filenames:
            k, P, t = self.load_or_calculate_spectrum(filename, args)
            k = k[:k.shape[0]/2]
            P = P[:P.shape[0]/2]
            y = P * k**args.compensate

            if args.smooth:
                W = args.smooth
                y = smooth(y, window_len=(2*W+1))[W:-W]
            if not args.noshow: plt.loglog(k, y, ls='-', marker=None, mfc='none', label=r'$t=%3.2e$'%t)

        plt.xlabel(r"$k/k_1$", fontsize=18)
        plt.ylabel(r"$k^{%3.2f} P_k$"%args.compensate, fontsize=18)
        plt.xlim(0.66*k.min(), k.max()*1.5)


    def do_time_series(self, args):
        import matplotlib.pyplot as plt
        import numpy as np

        ks, Ps, ts = self.load_or_calculate_time_series(args)

        kimax = np.argmax(Ps, axis=1)
        Pmax = [Ps[ti,ki] for ti, ki in enumerate(kimax)]
        kmax = [ks[ti,ki] for ti, ki in enumerate(kimax)]

        kind = "Pmax"

        if kind == "kmax":
            plt.loglog(ts, kmax)
            plt.ylabel(r"$k_t$", fontsize=16)

        if kind == "Pmax":
            plt.loglog(ts, Pmax)
            plt.ylabel(r"$P(k_t,t)$", fontsize=16)

        if kind == "Pkt":
            try:
               k0s = [int(k) for k in args.bins.split(',')]
            except:
               k0s = range(10)

            for ki in k0s:
                plt.loglog(ts, Ps[:,ki], '-', label=r'$k=%3.2f$' % ks[0][ki])

            plt.ylabel(r"$P_k(t)$", fontsize=16)
            plt.legend(loc='best')

        plt.xlabel(r"$t$", fontsize=16)


    def load_or_calculate_time_series(self, args):
        import numpy as np

        if args.filenames[0].endswith('.npz'):
            # We assume filenames are time-series databases
            print "loading the time series database from", args.filenames[0]
            db = np.load(args.filenames[0])
            return db['ks'], db['Ps'], db['ts']

        elif args.filenames[0].endswith('.h5'):
            # We assume filenames are HDF5 checkpoints with 'spectra' group
            print "calculating the time series of spectra"
            ks, Ps, ts = [], [], []
            for filename in args.filenames:
                k, P, t = self.load_or_calculate_spectrum(filename, args)
                ks.append(k)
                Ps.append(P)
                ts.append(t)

            sorted_ind = np.argsort(ts)
            ks = np.array(ks)[sorted_ind]
            Ps = np.array(Ps)[sorted_ind]
            ts = np.array(ts)[sorted_ind]

            if args.save_series:
                print "saving the time series database as", args.save_series
                np.savez(args.save_series, ts=ts, ks=ks, Ps=Ps)

            return ks, Ps, ts


    def load_or_calculate_spectrum(self, filename, args):
        import h5py
        import numpy as np

        chkpt = checkpoint.Checkpoint(filename, "r+")
        t = chkpt.status['time_simulation']

        if args.recalc and "spectra" in chkpt._file:
            del chkpt._file["spectra"]

        if "spectra" in chkpt._file:
            spectra = chkpt._file["spectra"]
            if "magnetic" in spectra:
                magnetic = spectra["magnetic"]
                print "loading the cached spectrum for", filename
                return magnetic["k"][:], magnetic["P"][:], t

        print "calculating power spectrum for", filename
        B1 = chkpt.cell_primitive['B1'][...]
        B2 = chkpt.cell_primitive['B2'][...]
        B3 = chkpt.cell_primitive['B3'][...]
        B = np.zeros((3,) + B1.shape)
        B[0] = B1
        B[1] = B2
        B[2] = B3
        k, P = power_spectrum(B, self.make_bins(args))
        chkpt.close()

        if not args.nocache:
            print "caching spectrum"
            h5f = h5py.File(filename)
            spectra = h5f.require_group("spectra")
            magnetic = spectra.require_group("magnetic")
            magnetic["k"] = k
            magnetic["P"] = P
            h5f.close()

        return k, P, t

