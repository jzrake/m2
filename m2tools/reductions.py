import checkpoint
import command
from autolog import logmethod


class ReductionsReader(object):

    def __init__(self, filename):
        self.checkpoint = checkpoint.Checkpoint(filename)

    def energy(self, which):
        if which == 'total':
            return self.checkpoint.get_reduction('total_energy')
        else:
            return self.checkpoint.get_reduction('total_'+which+'_energy')



class Reductions(command.Command):
    _alias = "reductions"


    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+')
        parser.add_argument("--time-marker", default=None, type=float)


    def run(self, args):
        import matplotlib.pyplot as plt

        for filename in args.filenames:
            self.run_file(args, filename)

        plt.legend(loc='best')
        plt.show()


    def run_file(self, args, filename, ax=None):
        import matplotlib.pyplot as plt
        import numpy as np

        reduction = ReductionsReader(filename)
        t, Em = reduction.energy('magnetic')
        t, Ei = reduction.energy('internal')
        t, Ek = reduction.energy('kinetic')
        t, Et = reduction.energy('total')
        t, Mt = reduction.checkpoint.get_reduction('total_mass')

        L = reduction.checkpoint.domain_extent_upper[1] - reduction.checkpoint.domain_extent_lower[1]
        t /= L
            
        if ax is None:
            ax = plt.gca()

        try:
            if args.time_marker is not None:
                t1 = args.time_marker / L
                ax.axvline(t1, ls='--')
        except AttributeError:
            pass

        plot = ax.semilogy
        plot(t, Ek, c='g', lw=1.5, ls='-' , label='kinetic')
        plot(t, Em, c='b', lw=1.5, ls='-', label='magnetic')

        ax.set_xlim(t[0], t[-1])
        ax.set_xlabel(r'$t / L$')
        ax.set_ylabel(r'energy')
