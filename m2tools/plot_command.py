import command
import checkpoint



class PlotDriver(object):

    _hardcopy = None
    _plot_axes = None

    def show(self):
        import matplotlib.pyplot as plt
        if self._hardcopy:
            plt.savefig(self._hardcopy, dpi=600)
            plt.clf()
        else:
            plt.show()

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
        import matplotlib.pyplot as plt
        import numpy as np
        data = self._chkpt.cell_primitive[self._args.field][:,:]
        plt.imshow(data)
        plt.colorbar()
        plt.axis('equal')
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

        if self._chkpt.domain_resolution[3] == 1:
            print "r-theta in axial symmetry"
            data = self._chkpt.cell_primitive[self._args.field][1:,1:]
            R, T = self._chkpt.cell_edge_meshgrid
            X = +R * np.sin(T)
            Z = +R * np.cos(T)
            plt.pcolormesh(X, Z, data, vmin=args.vmin, vmax=args.vmax)
            X = -R * np.sin(T)
            Z = +R * np.cos(T)
            plt.pcolormesh(X, Z, data, vmin=args.vmin, vmax=args.vmax)
        elif not args.equatorial:
            print "r-theta slice"
            data = self._chkpt.cell_primitive[self._args.field][1:,1:,0]
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
            data = self._chkpt.cell_primitive[self._args.field][1:,nt/2,:]
            R, T, P = self._chkpt.cell_edge_meshgrid
            X = R * np.sin(T) * np.cos(P)
            Y = R * np.sin(T) * np.sin(P)
            X = X[:,nt/2,:]
            Y = Y[:,nt/2,:]
            plt.pcolormesh(X, Y, data, vmin=args.vmin, vmax=args.vmax)

        plt.colorbar()
        plt.axis('equal')
        self.show()



class PlotCommand(command.Command):
    _alias = "plot"

    def configure_parser(self, parser):
        parser.add_argument('filename')
        parser.add_argument('--field', default='d')
        parser.add_argument('--vmin', type=float, default=None)
        parser.add_argument('--vmax', type=float, default=None)
        parser.add_argument('--equatorial', action="store_true")
        parser.add_argument('-o', '--output')

    def run(self, args):
        chkpt = checkpoint.Checkpoint(args.filename)
        if chkpt.geometry == 'cartesian':
            if (chkpt.domain_resolution == 1).sum() == 2:
                plotter = ShocktubePlot1d(chkpt, args)
            elif (chkpt.domain_resolution == 1).sum() == 1:
                plotter = RectangularPlot2d(chkpt, args)

        elif chkpt.geometry == 'spherical':
            plotter = PolarPlot(chkpt, args)

        plotter.set_hardcopy(args.output)
        plotter.plot()
