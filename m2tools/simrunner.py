import subprocess
import shlex
import os
import h5py


class SimulationRunner(object):

    def __init__(self, problem, usempi=False, np=1):
        self._m2exe = 'src/m2'
        self._np = np
        self._usempi = usempi
        self._problem = problem
        self._m2args = dict(resolution=64,
                            tmax=0.1,
                            output_path='.',
                            relativistic=False,
                            magnetized=False,
                            hdf5_cadence=0.0,
                            tpl_cadence=0.0,
                            model_parameters="")

    def update_args(self, **kwargs):
        self._m2args.update(kwargs)

    def get_command(self):
        if self._usempi:
            mpi_cmd = "mpiexec -np {0} ".format(self._np)
        else:
            mpi_cmd = ""
        argstr = "{0}{1} main.lua {2}".format(mpi_cmd, self._m2exe, self._problem)
        for k,v in self._m2args.iteritems():
            if k == 'magnetized':
                argstr += ' --magnetized' if v else ' --unmagnetized'
            elif k == 'relativistic':
                argstr += ' --relativistic' if v else ' --newtonian'
            else:
                argstr += ' --{0}={1}'.format(k.replace('_', '-'), v)
        return argstr

    def run(self):
        cmd = self.get_command()        
        arg = shlex.split(cmd)
        p = subprocess.Popen(arg, stdout=subprocess.PIPE)
        p.communicate()

        try:
            h5f = h5py.File('chkpt.final.h5', 'r')
        except IOError:
            return None

        fields = { }
        for f in h5f['prim']:
            fields[f] = h5f['prim'][f][...]
        h5f.close()
        os.remove('chkpt.final.h5')

        return fields
