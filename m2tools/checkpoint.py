import command


def is_checkpoint_file(filename):
    import h5py
    try:
        h5file = h5py.File(filename, 'r')
        is_checkpoint = 'version' in h5file
        if is_checkpoint:
            is_checkpoint = h5file['version'].value.startswith('m2')
        h5file.close()
        return is_checkpoint
    except Exception as e:
        return False


class Slicer(object):
    """
    Helper class for creating a multi-dimensional slice
    """
    def __getitem__(self, key): return key


class SingleField(object):
    """
    Field transformation representing a single primitive variable field
    """
    def __init__(self, key):
        self._key = key
    def transform(self, prim_selection):
        return prim_selection[self._key]


class HookedDictionary(dict):
    def __init__(self, hook):
        self._hook = hook
    def __setitem__(self, key, val):
        self._hook(key, val)
        super(HookedDictionary, self).__setitem__(key, val)


class Checkpoint(object):
    """
    Class representing an m2 HDF5 checkpoint file
    """
    def __init__(self, filename, writable=False):
        import h5py
        self._file = None # so that instance has attribute even if open fails
        self._file = h5py.File(filename, None if writable else 'r')
        self._model_parameters = { }
        self._config = { }
        self._status = HookedDictionary(self._status_set_hook)
        self._enable_set_hook = False
        for k,v in self._file['model_parameters'].iteritems():
            self._model_parameters[k] = v.value
        for k,v in self._file['config'].iteritems():
            self._config[k] = v.value
        for k,v in self._file['status'].iteritems():
            self._status[k] = v.value
        self._enable_set_hook = True
        self._scaling = lambda x: x
        self._selection = slice(None)
        self._derived_fields = { }
        for k in self.cell_primitive:
            self._derived_fields[k] = SingleField(k)

    def __del__(self):
        if self._file:
            self._file.close()

    def __setitem__(self, key, val):
        if key in self._file: del self._file[key]
        self._file[key] = val

    def __repr__(self):
        return ("<m2 checkpoint {0}>"
                "\n\tproblem ... {1}"
                "\n\tversion ... {2}"
                "\n\ttime_stamp ... {3}"
                "\n\tresolution ... {4}").format(
                    self.filename,
                    self.problem_name,
                    self.version[4:],
                    self.time_stamp,
                    self.domain_resolution[1:])

    def _status_set_hook(self, key, val):
        if not self._enable_set_hook:
            return
        if key in self._file['status']:
            del self._file['status'][key]
        self._file['status'][key] = val

    def _fromstring(self, s, dtype):
        import numpy as np
        return np.fromstring(s[1:-1], dtype=dtype, sep=' ')

    def close(self):
        self._file.close()
        self._file = None

    @property
    def build_date(self):
        return self._file['build_date'].value

    @property
    def time_stamp(self):
        return self._file['time_stamp'].value

    @property
    def version(self):
        return self._file['version'].value

    @property
    def filename(self):
        return self._file.filename

    @property
    def model_parameters(self):
        model_parameters = { }
        for k,v in self._model_parameters.iteritems():
            model_parameters[k] = v
        return model_parameters

    @property
    def config(self):
        config = { }
        for k,v in self._config.iteritems():
            config[k] = v
        return config

    @property
    def status(self):
        return self._status

    @property
    def problem_name(self):
        return self._file['problem_name'].value

    @property
    def periodic_dimension(self):
        cfg = self._config
        return self._fromstring(self._config['periodic_dimension'], int)

    @property
    def dimensionality(self):
        return (self.domain_resolution > 1).sum()

    @property
    def domain_resolution(self):
        cfg = self._config
        return self._fromstring(self._config['domain_resolution'], int)

    @property
    def domain_extent_upper(self):
        cfg = self._config
        return self._fromstring(self._config['domain_extent_upper'], float)

    @property
    def domain_extent_lower(self):
        cfg = self._config
        return self._fromstring(self._config['domain_extent_lower'], float)

    @property
    def coordinate_scalings(self):
        lookup = {'4': 'linear', '5': 'logarithmic'}
        return (None,
                lookup[self._config['coordinate_scaling1']],
                lookup[self._config['coordinate_scaling2']],
                lookup[self._config['coordinate_scaling3']])

    @property
    def geometry(self):
        lookup = {'0': 'cartesian',
                  '1': 'cylindrical',
                  '2': 'spherical',
                  '3': 'parabolic'}
        return lookup[self._config['geometry']]

    @property
    def cell_edge_coordinates(self):
        import numpy as np
        x0 = self.domain_extent_lower
        x1 = self.domain_extent_upper
        ns = self.domain_resolution
        ss = self.coordinate_scalings
        fs = { 'linear': [np.linspace, lambda x: x],
               'logarithmic': [np.logspace, np.log10]}
        return tuple([
            fs[ss[i]][0](fs[ss[i]][1](x0[i]),
                         fs[ss[i]][1](x1[i]),
                         ns[i] + 1) for i in [1,2,3]])

    @property
    def cell_edge_meshgrid(self):
        import numpy as np
        x1, x2, x3 = self.cell_edge_coordinates
        args = [ ]
        if len(x1) > 2: args.append(x1)
        if len(x2) > 2: args.append(x2)
        if len(x3) > 2: args.append(x3)
        return np.ix_(*args)

    @property
    def cell_primitive(self):
        return dict([k,v] for k,v in self._file['cell_primitive'].iteritems())

    @property
    def face_magnetic_flux(self):
        return dict([int(k),v] for k,v in
                    self._file['face_magnetic_flux'].iteritems())

    def add_derived_field(self, key, transform):
        self._derived_fields[key] = transform

    def set_scaling(self, scaling):
        self._scaling = scaling

    def set_selection(self, selection):
        """
        Set the selection used by the get_field method to retrieve data from the
        HDF5 file. Selection is a list of slices.
        """
        self._selection = selection

    def get_field(self, key):
        transform = self._derived_fields[key]
        prim_selection = { }
        for k in self.cell_primitive:
            prim_selection[k] = self.cell_primitive[k][self._selection]
        return self._scaling(transform.transform(prim_selection))



class CheckpointSummary(command.Command):
    _alias = "summarize"

    def configure_parser(self, parser):
        parser.add_argument('filenames', nargs='+')

    def run(self, args):
        for filename in args.filenames:
            chkpt = Checkpoint(filename)
            print chkpt

