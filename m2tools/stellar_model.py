import command
import checkpoint


class HegerStellarModel(object):
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


    def plot_profile(self):
        import matplotlib.pyplot as plt

        for field in self._args.fields.split(','):
            plt.figure()
            plt.loglog(self._table['radius']['v'], self._table[field]['v'],
                       '-o', mfc='none')
            plt.xlabel('radius [cm]')
            plt.ylabel(field + ' [cgs]')
        plt.show()


    def write_profile_to_hdf5(self):
        import numpy as np
        import h5py
        from scipy.interpolate import interp1d

        h5file = h5py.File(self._args.output, 'w')
        table = self._table

        r0 = np.log10(self._args.inner)
        r1 = np.log10(self._args.outer)
        N = self._args.points
        r_samp =  0.5 * (np.linspace(r0, r1, N+1)[1:] +
                         np.linspace(r0, r1, N+1)[:-1])
        r = np.log10(table['radius']['v'])

        for field in self._args.fields.split(','):
            v = table[field]['v']
            y = np.log10(v)
            interpolant = interp1d(r, y, kind='linear')

            y_samp = interpolant(r_samp)
            h5file[field] = 10**y_samp

        h5file['radius'] = 10**r_samp
        h5file['stellar_model_file'] = self._args.filename


    def write_profile_to_checkpoint(self):
        import numpy as np
        from scipy.interpolate import interp1d

        chkpt = checkpoint.Checkpoint(self._args.output, writable=True)
        table = self._table

        """ Todo: replace with physics units read from checkpoint """
        light_speed = 2.99792458e10 #cm/s
        u_length = 1e9 # cm
        u_density = 1 # g/cm^3
        u_time = u_length / light_speed
        u_pressure = u_density * light_speed**2
        u_velocity = u_length / u_time
        units = {
            'radius': u_length,
            'v_phi': u_velocity,
            'density': u_density,
            'pressure': u_pressure}
        dc = np.log10(float(chkpt.model_parameters['d_cavity']))
        rc = np.log10(float(chkpt.model_parameters['r_cavity']))
        r0 = np.log10(chkpt.domain_extent_lower[1])
        r1 = np.log10(chkpt.domain_extent_upper[1])
        N = chkpt.domain_resolution[1]
        r_samp =  0.5 * (np.linspace(r0, r1, N+1)[1:] +
                         np.linspace(r0, r1, N+1)[:-1])
        r = np.log10(table['radius']['v'] / units['radius'])

        for field_c, field_t in [('v3', 'v_phi'),
                                 ('d', 'density'),
                                 ('p', 'pressure')]:
            v = table[field_t]['v']
            u = units[field_t]
            y = np.log10(v / u)
            interpolant = interp1d(r, y, kind='linear')

            y_samp = interpolant(r_samp)

            if field_t == 'density':
                y_samp[r_samp < rc] = dc
            elif field_t == 'pressure':
                y_samp[r_samp < rc] = y_samp[r_samp >= rc][0]

            for i,y in enumerate(y_samp):
                chkpt.cell_primitive[field_c][i+1,...] = 10**y
        chkpt['stellar_model_file'] = self._args.filename
        chkpt.close()



class ReadStellarModelCommand(command.Command):
    _alias = "read-star"

    def configure_parser(self, parser):
        parser.add_argument("filename")
        parser.add_argument("--fields", default="density")
        parser.add_argument("--inner", default=1e9, type=float)
        parser.add_argument("--outer", default=1e11, type=float)
        parser.add_argument("--points", default=128, type=int)
        parser.add_argument("-o", "--output", default=None)

    def run(self, args):
        model = HegerStellarModel(args)
        if args.output is None:
            model.plot_profile()
        elif checkpoint.is_checkpoint_file(args.output):
            print "writing profile to checkpoint"
            model.write_profile_to_checkpoint()
        else:
            print "writing profile to hdf5 table"
            model.write_profile_to_hdf5()

