import command
import checkpoint


light_speed = 2.99792458e10 # cm/s


class StellarModel(object):

    def plot_profile(self):
        import matplotlib.pyplot as plt
        import numpy as np

        def cs(model):
            d = model.get_profile('d', model._radius)
            p = model.get_profile('p', model._radius)
            return (p / (d * light_speed**2))**0.5

        transforms = { 'cs': cs }

        for field in self._args.fields.split(','):
            plt.figure()
            if field in transforms:
                y = transforms[field](self)
            else:
                y = self.get_profile(field, self._radius)
            plt.loglog(self._radius, y, '-o', mfc='none')
            plt.xlabel('radius [cm]')
            plt.ylabel(field + ' [cgs]')
        plt.show()


    def write_profile_to_checkpoint(self):
        import numpy as np
        from scipy.interpolate import interp1d

        chkpt = checkpoint.Checkpoint(self._args.output, writable=True)

        u_length = float(chkpt.model_parameters['r_inner'])
        u_density = 1 # g/cm^3
        u_time = u_length / light_speed
        u_velocity = u_length / u_time
        u_pressure = u_density * u_velocity**2 
        units = {
            'r': u_length,
            'd': u_density,
            'p': u_pressure,
            'v1': u_velocity,
            'v3': u_velocity}
        v1 = float(chkpt.model_parameters['v_star'])
        dc = float(chkpt.model_parameters['d_cavity'])
        pc = float(chkpt.model_parameters['p_cavity'])
        rc = float(chkpt.model_parameters['r_cavity'])
        r0 = chkpt.domain_extent_lower[1]
        r1 = chkpt.domain_extent_upper[1]
        N = chkpt.domain_resolution[1]
        r_samp =  0.5 * (np.logspace(np.log10(r0), np.log10(r1), N+1)[1:] +
                         np.logspace(np.log10(r0), np.log10(r1), N+1)[:-1])

        for field in ('v1', 'v3', 'p', 'd'):
            y = self.get_profile(field, r_samp*u_length) / units[field]

            if field == 'd':
                y[r_samp < rc] = dc
            elif field == 'p':
                if pc < 0.0:
                    y[r_samp < rc] = y[r_samp >= rc][0]
                else:
                    y[r_samp < rc] = pc
            elif field == 'v1':
                y[(r_samp > rc) * (r_samp < 2 * rc)] = v1

            for i,yi in enumerate(y):
                chkpt.cell_primitive[field][i+1,...] = yi

        chkpt['stellar_model_file'] = self._args.filename
        chkpt.close()



class HegerStellarModel(StellarModel):
    """
    column	data	unit
    grid	cell number	---
    cell outer total mass	mass coordinate at the cell top interface   g
    cell outer radius	radius coordinate at the cell top interface	cm
    cell outer velocity	velocity of the cell top interface	cm s-1
    cell density	average density of the cell	g cm-3
    cell temperature	average temperature of the cell	K
    cell pressure	average pressure of the cell	dyn cm-2
    cell specific energy	average specific energy of the cell erg g-1
    cell specific entropy	average specific entropy of cell    kB baryon-1
    cell angular velocity	average angular velocity of the cell    rad s-1
    cell A_bar	average mean mass number of nuclei in the cell	1
    cell Y_e	average electrons per baryon in the cell	1
    stability	hydrodynamic stability of outer cell boundary	---
    B_r	poloidal magnetic field at outer cell boundary	Gauss
    B_phi	azimuthal magnetic field at outer cell boundary	Gauss
    neutrons	neutrons	neutrons	(mass fraction)
    H1	1H	protons	(mass fraction)
    He3	3He	---	(mass fraction)
    He4	4He	1 < A < 6	(mass fraction)
    C12	12C	---	(mass fraction)
    N14	14N	---	(mass fraction)
    O16	16O	16O (unburned yet)	(mass fraction)
    Ne20	20Ne	---	(mass fraction)
    Mg24	24Mg	22 < A < 29, excluding 28Si	(mass fraction)
    Si28	28Si	28Si	(mass fraction)
    S32	32S	28 < A < 36	(mass fraction)
    Ar36	36Ar	35 < A < 40	(mass fraction)
    Ca40	40Ca	39 < A < 44	(mass fraction)
    Ti44	44Ti	43 < A < 48	(mass fraction)
    Cr48	48Cr	47 < A < 52	(mass fraction)
    Fe52	52Fe	---	(mass fraction)
    Fe54	54Fe (+56Fe)	A approximately 2Z+2, Iron Peak (mass fraction)
    Ni56	56Ni	A < 2*Z+2, Iron Peak	(mass fraction)
    Fe56	---	56Fe only	(mass fraction)
    'Fe'	---	A > 2*Z + 3, Iron Peak, excluding 56Fe	(mass fraction)
    """
    def __init__(self, args):
        self._args = args
        self._parse_table()

    def _parse_table(self):
        import numpy as np

        table = {
            'M': dict(i=1, v=[]), # enclosed mass
            'radius': dict(i=2, v=[]),
            'vr': dict(i=3, v=[]),
            'density': dict(i=4, v=[]),
            'pressure': dict(i=6, v=[]),
            'omega': dict(i=9, v=[]),
        }
        table_rows = [ ]

        f = open(self._args.filename)
        nrows = 0
        for line in f.readlines()[2:]:
            cols = line.split()
            for column, field in table.iteritems():
                field['v'].append(float(cols[field['i']]))
            nrows += 1

        for n in range(nrows):
            d = dict()
            for k in table:
                d[k] = table[k]['v'][n]
            table_rows.append(d)

        for name, field in table.iteritems():
            field['v'] = np.array(field['v'])

        def v_phi(row):
            r = row['radius']
            w = row['omega']
            return w*r

        def j_specific(row):
            r = row['radius']
            w = row['omega']
            return w*r*r

        table['v_phi'] = dict(v=np.array([v_phi(r)
                                          for r in table_rows]))
        table['j_specific'] = dict(v=np.array([j_specific(r)
                                               for r in table_rows]))
        self._table = table


    @property
    def _radius(self):
        return self._table['radius']['v']


    def get_profile(self, field, r):
        import numpy as np
        from scipy.interpolate import interp1d
        field = {'d': 'density', 'p': 'pressure', 'v3': 'v_phi'}[field]
        v = self._table[field]['v']
        y = np.log10(v)
        r = np.log10(r)
        interpolant = interp1d(r, y, kind='linear')

        return 10**interpolant(r)



class MesaDuffellStellarModel(StellarModel):
    def __init__(self, args):
        self._args = args
        self._tabulate_model()

    def _tabulate_model(self):
        import numpy as np
        from scipy.interpolate import interp1d
        from scipy.integrate import cumtrapz

        R_outer = 1e10 # cm
        M_star = 15 * 1.9891e33 # g (15 solar mass)

        R = 1.0 # in units of 1e10 cm
        M = 1.0 # in units of ~0.45 * 15 solar masses (total mass is ~0.54M)

        u_velocity = light_speed
        u_mass = M_star / 0.45
        u_length = R_outer
        u_density = u_mass / u_length**3
        u_pressure = u_density * light_speed**2 

        rho0 = 1e-6
        a1 = 23.4*(0.4/R)
        b1 = 0.058*(R/0.4)
        M1 = 11.27/11.2*M
        a2 = 472.*(0.4/R)
        b2 = 0.002*(R/0.4)
        M2 = 3.24/11.2*M

        def f(r,a,b,M):
            return (M / np.pi / a) / ((r-b)**2 + a**-2)

        def density(r):
            return ((f(r,a1,b1,M1) + f(r,a2,b2,M2))/(4*np.pi*r**2) *
                    (1 - r/R)**3.85 + rho0) * u_density

        def pressure(r):
            return 0.00025 * density(r) * light_speed**2

        r = np.logspace(-4, 0.0, 128)
        p = [pressure(ri) for ri in r]
        d = [density(ri) for ri in r]
        dMdr = [density(ri) * 4*np.pi*(ri*u_length)**2 for ri in r]
        m = cumtrapz(dMdr, x=r*u_length)
        m = np.concatenate([[m[0]], m])
        logfuncs = { }
        logfuncs['p'] = interp1d(np.log10(r*u_length), np.log10(p), kind='linear')
        logfuncs['d'] = interp1d(np.log10(r*u_length), np.log10(d), kind='linear')
        logfuncs['m'] = interp1d(np.log10(r*u_length), np.log10(m), kind='linear')

        self._logfuncs = logfuncs
        self._radius = r * u_length


    def get_profile(self, field, r):
        import numpy as np
        if field in ['v1', 'v3']: return np.zeros_like(r)
        else: return 10**self._logfuncs[field](np.log10(r))



class ReadStellarModelCommand(command.Command):
    _alias = "make-star"

    def configure_parser(self, parser):
        parser.add_argument("filename")
        parser.add_argument("--fields", default="d")
        parser.add_argument("--inner", default=1e9, type=float)
        parser.add_argument("--outer", default=1e11, type=float)
        parser.add_argument("--points", default=128, type=int)
        parser.add_argument("-o", "--output", default=None)

    def run(self, args):
        if args.filename == 'mesa-duffell':
            model = MesaDuffellStellarModel(args)
        else:
            model = HegerStellarModel(args)
        if args.output is None:
            model.plot_profile()
        elif checkpoint.is_checkpoint_file(args.output):
            model.write_profile_to_checkpoint()
        else:
            raise RuntimeError("output file must be an m2 checkpoint")

