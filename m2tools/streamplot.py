import command
import checkpoint



class StreamlinePlot(command.Command):
    _alias = "streamplot"

    def configure_parser(self, parser):
        parser.add_argument('filenames', nargs='+')
        parser.add_argument('--hardcopy', action='store_true')
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
                imgname = "%s/frame-%04d.jpg" % (os.path.dirname(filename), self.frame_number)

            print "writing", imgname
            plt.savefig(imgname)
            plt.clf()

        chkpt.close()
