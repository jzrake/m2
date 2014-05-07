

class SingleField(object):
    def __init__(self, key):
        self._key = key
    def transform(self, prim_selection):
        return prim_selection[self._key]


class Checkpoint(object):
    def __init__(self, filename):
        import h5py
        self._file = None # so that instance has attribute even if open fails
        self._file = h5py.File(filename, 'r')
        self._config = { }
        self._status = { }
        for k,v in self._file['config'].iteritems():
            self._config[k] = v.value
        for k,v in self._file['status'].iteritems():
            self._status[k] = v.value
        self._selection = slice(None)
        self._derived_fields = { }
        for k in self.cell_primitive:
            self._derived_fields[k] = SingleField(k)

    def __del__(self):
        if self._file:
            self._file.close()

    def _fromstring(self, s, dtype):
        import numpy as np
        return np.fromstring(s[1:-1], dtype=dtype, sep=' ')

    def close(self):
        self._file.close()
        self._file = None

    @property
    def config(self):
        config = { }
        for k,v in self._config.iteritems():
            config[k] = v
        return config

    @property
    def status(self):
        status = { }
        for k,v in self._status.iteritems():
            status[k] = v
        return status

    @property
    def periodic_dimension(self):
        cfg = self._config
        return self._fromstring(self._config['periodic_dimension'], int)

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
        return transform.transform(prim_selection)


if __name__ == "__main__":
    chkpt = Checkpoint('chkpt.final.h5')

    print chkpt.domain_resolution
    print chkpt.domain_extent_lower
    print chkpt.domain_extent_upper
    print chkpt.coordinate_scalings
    print chkpt.cell_edge_coordinates[:2]
    print chkpt.periodic_dimension
    print chkpt.geometry
