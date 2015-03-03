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
        parser.add_argument("--helicity", type=str, default=None)

    def run(self, args):
        import matplotlib.pyplot as plt
        self.first = True
        self.growth_rates = [ ]
        self.resolutions = [ ]
        for filename in args.filenames:
            self.fit = False
            self.run_file(filename, args)

        if args.helicity:
            import pickle
            tH = pickle.load(open(args.helicity, 'r'))
            plt.semilogy([a[0]/2.0 for a in tH], [a[1]*8*26**0.5 for a in tH], '-o', label=r'$H k_0$')

        try:
            k02 = int(self.last_checkpoint_used.model_parameters['k2'])
        except:
            k02 = 11
        model = int(self.last_checkpoint_used.model_parameters['model'])
        plt.legend(loc='best')
        plt.title(r"$\tilde{\alpha}^2=%d, \ \rm{model}=%d$"%(k02, model))
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

#         twin_vals = {64:[53,58],
#                      96:[56,60],
#                      128:[38,44],
#                      192:[30,37],
#                      256:[26,28]}

        if self.fit:
            #n0 = np.where(t < twin_vals[N][0])[0][-1]
            #n1 = np.where(t < twin_vals[N][1])[0][-1]
            if filename.find('256') == -1:
                n0 = np.where(Ek < 5e-7)[0][-1]
                n1 = np.where(Ek < 1e-5)[0][-1]
            else:
                n0 = np.where(Ek < 2e-8)[0][-1]
                n1 = np.where(Ek < 2e-7)[0][-1]
                
            tw = t [n0:n1]
            Ew = Ek[n0:n1]
            a, b = np.polyfit(tw, np.log(Ew), 1)
            Ef = np.exp(a*t[n0:n1+10] + b)
            plot(t[n0:n1+10]+1, Ef, c='r')
            plt.text(tw[0]+1, Ef[0], r"$\tau^{-1}=%3.2f c/L$"%a, fontsize=16)
            self.resolutions.append(N)
            self.growth_rates.append(a)
            
        plt.text(t[-1], Ek[-1], r"$%d^{3}$"%N, fontsize=16)
        plt.xlabel('time (light-crossing)')
        plt.ylabel('energy')
        self.last_checkpoint_used = reduction.checkpoint
