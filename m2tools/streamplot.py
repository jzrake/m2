import os
import matplotlib.pyplot as plt
import numpy as np
import command
import checkpoint



class Plot2d(command.Command):
    _alias = "plot2d"

    def configure_parser(self, parser):
        parser.add_argument('filenames', nargs='+')
        parser.add_argument('--reductions', action='store_true')
        parser.add_argument('--kind', default="stream", choices=["stream", "lic", "relief", "current"])
        parser.add_argument('--hardcopy', type=str, default=None,
                            help="image format extension to use, onscreen if None")
        parser.add_argument('-o', '--output', default=None, help="override default image name")


    def run(self, args):
        self.frame_number = 0

        def make_axes(reductions=False):
            if reductions:
                fig = plt.figure(figsize=[8,12])
                ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2)
                ax2 = plt.subplot2grid((3,1), (2,0), rowspan=1)
                return ax1, ax2
            else:
                fig = plt.figure(figsize=[8,10])
                ax1 = fig.add_subplot(1,1,1)
                return ax1, None

        # If hardcopy, we will be re-using the same figure instance
        if args.hardcopy:
            ax1, ax2 = make_axes(args.reductions)

        for filename in args.filenames:
            # If onscreen, we will be creating a new figure for each frame
            if not args.hardcopy:
                ax1, ax2 = make_axes(args.reductions)
            self.run_file(args, filename, ax1, ax2)

        if not args.hardcopy:
            plt.show()


    def run_file(self, args, filename, ax1, ax2):
        method = {"stream": self.run_stream,
                  "lic": self.run_lic,
                  "relief": self.run_relief,
                  "current": self.run_current,
                  "current_diff": self.run_current_diff}[args.kind]

        chkpt = checkpoint.Checkpoint(filename)

        B1 = chkpt.cell_primitive['B1'][...].T
        B2 = chkpt.cell_primitive['B2'][...].T
        B3 = chkpt.cell_primitive['B3'][...].T
        x, y = chkpt.cell_center_meshgrid
        X = (1*x + 0*y).T
        Y = (0*x + 1*y).T
        B = (B1**2 + B2**2 + B3**2)**0.5

        method(ax1, B1, B2, B3, B, X, Y, args)

        if args.reductions:
            from m2tools.reductions import Reductions
            args.time_marker = chkpt.status['time_simulation']
            red = Reductions()
            red.run_file(args, args.filenames[-1], ax2)

        if ax2:
            ax2.set_ylim(1e-3, 4e1)
            ax2.legend(loc='lower center')
        plt.subplots_adjust(hspace=0)

        if args.hardcopy:
            self.frame_number += 1
            if args.output:
                imgname = args.output
            else:
                dirname = os.path.dirname(filename)
                imgname = "%s/frame-%04d.%s" % (dirname if dirname else ".",
                                                self.frame_number,
                                                args.hardcopy)

            print "writing", imgname
            plt.savefig(imgname)
            if ax1: ax1.clear()
            if ax2: ax2.clear()

        chkpt.close()


    def run_stream(self, ax, B1, B2, B3, B, X, Y, args):
        ax.set_title('in-page magnetic field')
        ax.streamplot(X, Y, B1, B2, density=2, color=B3, linewidth=B)
        ax.axis('equal')
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)


    def run_lic(self, ax, B1, B2, B3, B, X, Y, args):
        import lic.lic_internal

        Nx, Ny = B.shape
        kernel  = np.linspace(0, 1, 111).astype(np.float32)
        texture = np.random.rand(Nx, Ny).astype(np.float32)
        vectors = np.zeros([Nx, Ny, 2], dtype=np.float32)

        vectors[:,:,0] = B1
        vectors[:,:,1] = B2
        image = lic.lic_internal.line_integral_convolution(vectors, texture, kernel)

        ax.set_title('in-page magnetic field')
        ax.imshow(image, cmap='bone')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)


    def run_relief(self, ax, B1, B2, B3, B, X, Y, args):
        ax.set_title('out-of-page magnetic field')
        ax.imshow(B3, extent=[-1,1,-1,1], cmap='bone')
        ax.axis('equal')
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)


    def run_current_diff(self, ax, B1, B2, B3, B, X, Y, args):

        def five_point_deriv(f, axis=0, h=1.0):
            return (-1*np.roll(f, -2, axis) +
                    +8*np.roll(f, -1, axis) + 
                    -8*np.roll(f, +1, axis) + 
                    +1*np.roll(f, +2, axis)) / (12.0 * h)

        d = lambda f, a: five_point_deriv(f, a)
        jz = d(B2, 1) - d(B1, 0) # I know this looks wrong, but B1,B2,B3 are transposed

        if not hasattr(self, "_jz0"): self._jz0 = jz

        ax.set_title('$\delta j_z$ (perturbed current out-of-page)')
        ax.imshow(jz - self._jz0, extent=[-1,1,-1,1], cmap='bone')
        ax.axis('equal')
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)


    def run_current(self, ax, B1, B2, B3, B, X, Y, args):

        def five_point_deriv(f, axis=0, h=1.0):
            return (-1*np.roll(f, -2, axis) +
                    +8*np.roll(f, -1, axis) + 
                    -8*np.roll(f, +1, axis) + 
                    +1*np.roll(f, +2, axis)) / (12.0 * h)

        d = lambda f, a: five_point_deriv(f, a)
        jz = d(B2, 1) - d(B1, 0) # I know this looks wrong, but B1,B2,B3 are transposed

        ax.set_title('$j_z$ (current out-of-page)')
        ax.imshow(jz, extent=[-1,1,-1,1], cmap='bone')
        ax.axis('equal')
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
