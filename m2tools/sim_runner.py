import subprocess
import shlex
import os
import command
import checkpoint



class SimulationRunner(object):

    def __init__(self, problem, usempi=False, np=1, suppress_output=True):
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
        self._suppress_output = suppress_output

    def update_args(self, **kwargs):
        self._m2args.update(kwargs)

    def get_command(self):
        if self._usempi:
            mpi_cmd = "mpiexec -np {0} ".format(self._np)
        else:
            mpi_cmd = ""
        argstr = "{0}{1} {2}".format(mpi_cmd, self._m2exe, self._problem)
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
        stdout = subprocess.PIPE if self._suppress_output else None
        p = subprocess.Popen(arg, stdout=stdout)
        p.communicate()

        try:
            chkpt = checkpoint.Checkpoint('chkpt.final.h5')
        except IOError:
            return None

        fields = { }
        for f in chkpt.cell_primitive:
            fields[f] = chkpt.cell_primitive[f][...]
        chkpt.close()
        os.remove('chkpt.final.h5')

        return fields
