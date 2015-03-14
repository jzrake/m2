import command
import checkpoint






class LineIntegralConvolutionPlot(command.Command):
    _alias="lic"

    def configure_parser(self, parser):
        parser.add_argument('filenames', nargs='+')
        parser.add_argument('--hardcopy', type=str, default=None,
                            help="image format extension to use, onscreen if None")
        parser.add_argument('-o', '--output', default=None)

    def run(self, args):
        import matplotlib.pyplot as plt

        self.frame_number = 0
        for filename in args.filenames:
            self.run_file(args, filename)

        if not args.hardcopy:
            plt.show()

    def run_file(self, args, filename):
        import matplotlib.pyplot as plt
        import numpy as np
        import os
        import lic.lic_internal

        chkpt = checkpoint.Checkpoint(filename)
        B1 = chkpt.cell_primitive['B1'][...].T
        B2 = chkpt.cell_primitive['B2'][...].T
        B3 = chkpt.cell_primitive['B3'][...].T
        x, y = chkpt.cell_center_meshgrid
        X = (1*x + 0*y).T
        Y = (0*x + 1*y).T
        B = (B1**2 + B2**2 + B3**2)**0.5

        Nx, Ny = B.shape
        kernel  = np.linspace(0, 1, 111).astype(np.float32)
        texture = np.random.rand(Nx, Ny).astype(np.float32)
        vectors = np.zeros([Nx, Ny, 2], dtype=np.float32)

        vectors[:,:,0] = B1
        vectors[:,:,1] = B2
        image = lic.lic_internal.line_integral_convolution(vectors, texture, kernel)

        if not args.hardcopy:
            plt.figure()

        plt.imshow(image, cmap='bone')

        if args.hardcopy:
            self.frame_number += 1
            if args.output:
                imgname = args.output
            else:
                imgname = "%s/frame-%04d.%s" % (os.path.dirname(filename),
                                                self.frame_number,
                                                args.hardcopy)

            print "writing", imgname
            plt.savefig(imgname)
            plt.clf()

        chkpt.close()



class ReliefPlot(command.Command):
    _alias = "relief"

    def configure_parser(self, parser):
        parser.add_argument('filenames', nargs='+')
        parser.add_argument('--hardcopy', type=str, default=None,
                            help="image format extension to use, onscreen if None")
        parser.add_argument('-o', '--output', default=None)

    def run(self, args):
        import matplotlib.pyplot as plt

        self.frame_number = 0
        for filename in args.filenames:
            self.run_file(args, filename)

        if not args.hardcopy:
            plt.show()

    def run_file(self, args, filename):
        import matplotlib.pyplot as plt
        import os

        chkpt = checkpoint.Checkpoint(filename)
        B1 = chkpt.cell_primitive['B1'][...].T
        B2 = chkpt.cell_primitive['B2'][...].T
        B3 = chkpt.cell_primitive['B3'][...].T

        if not args.hardcopy:
            plt.figure()

        plt.imshow(B3, extent=[-1,1,-1,1], cmap='bone')
        plt.axis('equal')
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)


        if args.hardcopy:
            self.frame_number += 1
            if args.output:
                imgname = args.output
            else:
                imgname = "%s/frame-%04d.%s" % (os.path.dirname(filename),
                                                self.frame_number,
                                                args.hardcopy)

            print "writing", imgname
            plt.savefig(imgname)
            plt.clf()

        chkpt.close()
class StreamlinePlot(command.Command):
    _alias = "streamplot"

    def configure_parser(self, parser):
        parser.add_argument('filenames', nargs='+')
        parser.add_argument('--hardcopy', type=str, default=None,
                            help="image format extension to use, onscreen if None")
        parser.add_argument('-o', '--output', default=None)

    def run(self, args):
        import matplotlib.pyplot as plt

        self.frame_number = 0
        for filename in args.filenames:
            self.run_file(args, filename)

        if not args.hardcopy:
            plt.show()

    def run_file(self, args, filename):
        import matplotlib.pyplot as plt
        import os

        chkpt = checkpoint.Checkpoint(filename)
        B1 = chkpt.cell_primitive['B1'][...].T
        B2 = chkpt.cell_primitive['B2'][...].T
        B3 = chkpt.cell_primitive['B3'][...].T
        x, y = chkpt.cell_center_meshgrid
        X = (1*x + 0*y).T
        Y = (0*x + 1*y).T
        B = (B1**2 + B2**2 + B3**2)**0.5

        if not args.hardcopy:
            plt.figure()

        plt.streamplot(X, Y, B1, B2, density=2, color=B3, linewidth=B)
        plt.axis('equal')
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)


        if args.hardcopy:
            self.frame_number += 1
            if args.output:
                imgname = args.output
            else:
                imgname = "%s/frame-%04d.%s" % (os.path.dirname(filename),
                                                self.frame_number,
                                                args.hardcopy)

            print "writing", imgname
            plt.savefig(imgname)
            plt.clf()

        chkpt.close()
