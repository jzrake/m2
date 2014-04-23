import command
import checkpoint



class ShocktubePlot1d(object):

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
        plt.show()



class RectangularPlot2d(object):
    
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
        plt.show()



class PolarPlot2d(object):

    def __init__(self, chkpt, args):
        self._chkpt = chkpt
        self._args = args

    def plot(self):
        import matplotlib.pyplot as plt
        import numpy as np
        R, T = self._chkpt.cell_edge_meshgrid
        X = R * np.sin(T)
        Z = R * np.cos(T)
        data = self._chkpt.cell_primitive[self._args.field][1:,1:].T
        plt.pcolormesh(X, Z, data)
        plt.colorbar()
        plt.axis('equal')
        plt.show()



class PlotCommand(command.Command):
    _alias = "plot"

    def configure_parser(self, parser):
        parser.add_argument('filename')
        parser.add_argument('--field', default='d')

    def run(self, args):
        chkpt = checkpoint.Checkpoint(args.filename)

        if chkpt.geometry == 'cartesian':
            if (chkpt.domain_resolution == 1).sum() == 2:
                plotter = ShocktubePlot1d(chkpt, args)
            elif (chkpt.domain_resolution == 1).sum() == 1:
                plotter = RectangularPlot2d(chkpt, args)

        elif chkpt.geometry == 'spherical':
            if (chkpt.domain_resolution == 1).sum() == 1:
                plotter = PolarPlot2d(chkpt, args)

        plotter.plot()
