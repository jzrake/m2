import checkpoint
import command
from autolog import logmethod


class ReductionsReader(object):
    @logmethod
    def __init__(self, filename):
        self.checkpoint = checkpoint.Checkpoint(filename)

    @logmethod
    def energy(self, which):
        if which == 'total':
            return self.checkpoint.get_reduction('total_energy')
        else:
            return self.checkpoint.get_reduction('total_'+which+'_energy')


class PlotReductions(command.Command):
    _alias = "reductions"
    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+')

    def run(self, args):
        import matplotlib.pyplot as plt
        for filename in args.filenames:
            self.run_file(filename, args)
        plt.legend(loc='best')
        plt.show()

    def run_file(self, filename, args):
        import matplotlib.pyplot as plt
        reduction = ReductionsReader(filename)
        t, Em = reduction.energy('magnetic')
        t, Ei = reduction.energy('internal')
        t, Ek = reduction.energy('kinetic')
        t, Et = reduction.energy('total')
        t, Mt = reduction.checkpoint.get_reduction('total_mass')
        plt.xlim(t[0], t[-1])
        plt.plot(t, (Ek-Ek[0]), label='kinetic')
        plt.plot(t, (Em-Em[0]), label='magnetic')
        #plt.plot(t, Ei-Ei[0], label='internal')
        #plt.plot(t, (Et-Et[0])/Et[0], label='total')
        #plt.plot(t, (Mt-Mt[0])/Mt[0], label='mass')
