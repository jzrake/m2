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
        self.first = True
        for filename in args.filenames:
            self.fit = filename == args.filenames[0]
            self.run_file(filename, args)
    
        plt.legend(loc='best')
        plt.title(r"Free decay of FF equilibrium: $k_0=11, \rm{model}=2$")
        plt.show()

    def run_file(self, filename, args):
        import matplotlib.pyplot as plt
        import numpy as np

        reduction = ReductionsReader(filename)
        t, Em = reduction.energy('magnetic')
        t, Ei = reduction.energy('internal')
        t, Ek = reduction.energy('kinetic')
        t, Et = reduction.energy('total')
        t, Mt = reduction.checkpoint.get_reduction('total_mass')

        t /= 2.0
        N = reduction.checkpoint.domain_resolution[1]

        plot = plt.semilogy
        plot(t, Ek, c='k', lw=N/64.0, ls='-' , label='kinetic'  if self.first else None)
        plot(t, Em, c='k', lw=N/64.0, ls='-.', label='magnetic' if self.first else None)
        self.first = False

        if self.fit:
            n0 = np.where(t<3)[0][-1]
            n1 = np.where(t<5)[0][-1]
            tw = t [n0:n1]
            Ew = Ek[n0:n1]
            a, b = np.polyfit(tw, np.log(Ew), 1)
            Ef = np.exp(a*t[n0:n1+10] + b)
            plot(t[n0:n1+10]+1, Ef)
            plt.text(tw[0]+2, Ef[0], r"$e^{t/\tau}, \tau^{-1}=%3.2f L/c$"%a, fontsize=16)

        plt.text(t[-1], Ek[-1], r"$%d^{3}$"%N, fontsize=16)
        plt.xlabel('time (light-crossing)')
        plt.ylabel('energy')
