import command
import checkpoint
import transforms


slicer = checkpoint.Slicer()


class PlotDriver(object):

    _hardcopy = None
    _plot_axes = None
    _show_mode = 'show'

    def show(self):
        import matplotlib.pyplot as plt
        if self._show_mode == 'hold':
            pass
        elif self._hardcopy:
            plt.savefig(self._hardcopy, dpi=600)
            plt.clf()
        else:
            plt.show()

    def set_show_mode(self, mode):
        self._show_mode = mode

    def set_hardcopy(self, hardcopy):
        self._hardcopy = hardcopy

    def get_plot_axes(self, **fig_kwds):
        import matplotlib.pyplot as plt
        ax0 = self._plot_axes
        if ax0 is None:
            fig = plt.figure(**fig_kwds)
            ax0 = fig.add_subplot(111)
        return ax0



class ShocktubePlot1d(PlotDriver):

    def __init__(self, chkpt, args):
        self._chkpt = chkpt
        self._args = args

    def plot(self):
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import AxesGrid

        def trim_axis(ax):
            xl = ax.get_xlim()
            yl = ax.get_ylim()
            ax.set_xlim(xl[0]+0.001, xl[1]-0.001)
            ax.set_ylim(yl[0]+0.001, yl[1]-0.001)

        fig = plt.figure(1, (14, 10))
        fig.subplots_adjust(left=0.05, right=0.98)
        grid = AxesGrid(fig, 111,
                        nrows_ncols = (2, 4),
                        axes_pad = 0.05,
                        share_all=True,
                        add_all=True,
                        #direction='column',
                        label_mode='L')

        dsets = ['p','B1','B2','B3','d','v1','v2','v3']
        ymin, ymax = [], []
        for i, dset in enumerate(dsets):
            y = self._chkpt.cell_primitive[dset][...]
            ymin.append(y.min())
            ymax.append(y.max())

        spread = max(ymax) - min(ymin)

        for i, dset in enumerate(dsets):
            y = self._chkpt.cell_primitive[dset][...]
            x = np.linspace(0.0, spread, y.size)
            grid[i].plot(x, y, c='k', lw=2.0, marker='o', mfc='none')
            grid[i].text(0.1, 0.1, dset, transform=grid[i].transAxes, fontsize=16)
            grid[i].set_xlim(0.0, spread)

        trim_axis(grid[0])
        self.show()



class RectangularPlot2d(PlotDriver):
    
    def __init__(self, chkpt, args):
        self._chkpt = chkpt
        self._args = args

    def plot(self):
        import pprint
        import matplotlib.pyplot as plt
        import numpy as np

        pprint.pprint(self._chkpt.status)
        pprint.pprint(self._chkpt.config)

        args = self._args

        self._chkpt.set_selection(slicer[:,:])
        data = self._chkpt.get_field(args.field)

        self.get_plot_axes()
        plt.imshow(data, interpolation='nearest', origin='lower')
        plt.colorbar()
        plt.axis('equal')
        plt.title(self._chkpt.filename+'/'+args.field)
        self.show()



class RectangularPlot3d(PlotDriver):
    
    def __init__(self, chkpt, args):
        self._chkpt = chkpt
        self._args = args

    def plot(self):
        import pprint
        import matplotlib.pyplot as plt
        import numpy as np

        pprint.pprint(self._chkpt.status)
        pprint.pprint(self._chkpt.config)

        args = self._args

        nz = self._chkpt.domain_resolution[3]
        self._chkpt.set_selection(slicer[:,:,nz/2])
        data = self._chkpt.get_field(args.field)

        self.get_plot_axes()
        plt.imshow(data, interpolation='nearest', origin='lower')
        plt.colorbar()
        plt.axis('equal')
        plt.title(self._chkpt.filename+'/'+args.field)
        self.show()



class PolarPlot(PlotDriver):

    def __init__(self, chkpt, args):
        self._chkpt = chkpt
        self._args = args

    def plot(self):
        import pprint
        import matplotlib.pyplot as plt
        import numpy as np

        pprint.pprint(self._chkpt.status)
        pprint.pprint(self._chkpt.config)

        args = self._args
        self.get_plot_axes()

        if self._chkpt.domain_resolution[3] == 1:
            print "r-theta in axial symmetry"
            self._chkpt.set_selection(slicer[1:,1:])
            data = self._chkpt.get_field(args.field)
            R, T = self._chkpt.cell_edge_meshgrid
            X = +R * np.sin(T)
            Z = +R * np.cos(T)
            plt.pcolormesh(X, Z, data, vmin=args.vmin, vmax=args.vmax)
            X = -R * np.sin(T)
            Z = +R * np.cos(T)
            plt.pcolormesh(X, Z, data, vmin=args.vmin, vmax=args.vmax)
        elif not args.equatorial:
            print "r-theta slice"
            self._chkpt.set_selection(slicer[1:,1:,0])
            data = self._chkpt.get_field(args.field)
            R, T, P = self._chkpt.cell_edge_meshgrid
            R = R[...,0]
            T = T[...,0]
            X = +R * np.sin(T)
            Z = +R * np.cos(T)
            plt.pcolormesh(X, Z, data, vmin=args.vmin, vmax=args.vmax)
            X = -R * np.sin(T)
            Z = +R * np.cos(T)
            plt.pcolormesh(X, Z, data, vmin=args.vmin, vmax=args.vmax)
        else:
            print "r-phi slice"
            nt = self._chkpt.domain_resolution[2]
            self._chkpt.set_selection(slicer[1:,nt/2,:])
            data = self._chkpt.get_field(self._args.field)
            R, T, P = self._chkpt.cell_edge_meshgrid
            X = R * np.sin(T) * np.cos(P)
            Y = R * np.sin(T) * np.sin(P)
            X = X[:,nt/2,:]
            Y = Y[:,nt/2,:]
            plt.pcolormesh(X, Y, data, vmin=args.vmin, vmax=args.vmax)

        plt.title(self._chkpt.filename+'/'+args.field)
        plt.colorbar()
        plt.axis('equal')
        self.show()



class PlotCommand(command.Command):
    _alias = "plot"

    def configure_parser(self, parser):
        parser.add_argument('filenames', nargs='+')
        parser.add_argument('--field', default='d')
        parser.add_argument('--vmin', type=float, default=None)
        parser.add_argument('--vmax', type=float, default=None)
        parser.add_argument('--equatorial', action="store_true")
        parser.add_argument('-o', '--output')

    def run(self, args):
        for filename in args.filenames:
            chkpt = checkpoint.Checkpoint(filename)
            chkpt.add_derived_field('gamma', transforms.LorentzFactor())
            print (filename)
            if chkpt.geometry == 'cartesian':
                if chkpt.dimensionality == 1:
                    plotter = ShocktubePlot1d(chkpt, args)

                elif chkpt.dimensionality == 2:
                    plotter = RectangularPlot2d(chkpt, args)

                elif chkpt.dimensionality == 3:
                    plotter = RectangularPlot3d(chkpt, args)

            elif chkpt.geometry == 'spherical':
                plotter = PolarPlot(chkpt, args)

            plotter.set_show_mode('hold')
            #plotter.set_hardcopy(args.output)
            plotter.plot()

        import matplotlib.pyplot as plt
        plt.show()

