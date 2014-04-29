import unittest
import sim_runner

TOL = 1e-13


class SerialToParallel(object):

    def set_up(self, problem, **kwargs):
        runner_s = sim_runner.SimulationRunner(problem, usempi=True, np=1)
        runner_p = sim_runner.SimulationRunner(problem, usempi=True, np=8)
        runner_s.update_args(**kwargs)
        runner_p.update_args(**kwargs)
        self.answer_s = runner_s.run()
        self.answer_p = runner_p.run()

    def test_same_serial_parallel(self):
        self.assertTrue(self.answer_p and
                        self.answer_s, 'm2 run did not complete')
        for f in self.answer_s:
            diff = self.answer_s[f] - self.answer_p[f]
            self.assertAlmostEqual(abs(diff).max(), 0.0, delta=TOL)



class RestartedToContinuous(object):

    def set_up(self, problem, np=8, **kwargs):
        # continuous
        runner = sim_runner.SimulationRunner(problem, usempi=True, np=np)
        runner.update_args(**kwargs)
        runner.update_args(tmax=0.2)
        self.answer1 = runner.run()

        # restarted
        runner = sim_runner.SimulationRunner(problem, usempi=True, np=np)
        runner.update_args(**kwargs)
        runner.update_args(tmax=0.1)
        runner.run(remove_chkpt=False)

        runner = sim_runner.SimulationRunner(problem, usempi=True, np=np)
        runner.update_args(**kwargs)
        runner.update_args(tmax=0.2, restart='chkpt.final.h5')
        self.answer2 = runner.run(remove_chkpt=True)

    def test_same_restarted_to_continuous(self):
        self.assertTrue(self.answer1 and
                        self.answer2, 'm2 run did not complete')
        for f in self.answer1:
            diff = self.answer1[f] - self.answer2[f]
            self.assertAlmostEqual(abs(diff).max(), 0.0, delta=TOL)



class TestSymmetry2D(object):

    def set_up(self, problem, usempi=True, np=1, **kwargs):
        runner = sim_runner.SimulationRunner(problem, usempi=True, np=np)
        runner.update_args(**kwargs)
        self.answer = runner.run()

    def test_symmetry(self):
        self.assertIsNotNone(self.answer, 'm2 run did not complete')
        for f in self.answer:
            N = self.answer[f].shape
            sx = self.answer[f][N[0]/2,:]
            sy = self.answer[f][:,N[1]/2]
            diff = sx - sy
            self.assertAlmostEqual(diff.max(), 0.0, delta=TOL)



class TestAntisymmetry2D(object):

    def set_up(self, problem, usempi=False, np=1, **kwargs):
        runner = sim_runner.SimulationRunner(problem, usempi=usempi, np=np)
        runner.update_args(**kwargs)
        self.answer = runner.run()

    def test_antisymmetry(self):
        self.assertIsNotNone(self.answer, 'm2 run did not complete')
        for f in 'dp':
            N = self.answer[f].shape
            sx = self.answer[f][N[0]/2,:]
            diff = sx - sx[::-1]
            self.assertAlmostEqual(diff.max(), 0.0, delta=TOL)



class TestTopDownSymmetry2D(object):

    def set_up(self, problem, usempi=True, np=1, **kwargs):
        runner = sim_runner.SimulationRunner(problem, usempi=True, np=np)
        runner.update_args(**kwargs)
        self.answer = runner.run()

    def test_symmetry(self):
        self.assertIsNotNone(self.answer, 'm2 run did not complete')
        for f in 'dp':
            sx = self.answer[f][1:,1:] # assumes non-periodic domain
            diff = sx - sx[:,::-1]
            self.assertAlmostEqual(diff.max(), 0.0, delta=TOL)


class BrioWuSerialToParallel(unittest.TestCase, SerialToParallel):
    def setUp(self):
        self.set_up('BrioWu', magnetized=False)

class BrioWuMagnetizedSerialToParallel(unittest.TestCase, SerialToParallel):
    def setUp(self):
        self.set_up('BrioWu', magnetized=True)

class RyuJonesSerialToParallel(unittest.TestCase, SerialToParallel):
    def setUp(self):
        self.set_up('RyuJones')

class DensityWaveSerialParallel(unittest.TestCase, SerialToParallel):
    def setUp(self):
        self.set_up('DensityWave')

class DensityWaveSymmetry2D(unittest.TestCase, TestSymmetry2D):
    def setUp(self):
        self.set_up('DensityWave', model_parameters='dims=2')

class DensityWaveRestartedToContinuous(unittest.TestCase, TestSymmetry2D):
    def setUp(self):
        self.set_up('DensityWave', model_parameters='dims=2')

class DensityWaveMagnetizedSymmetry(unittest.TestCase, TestSymmetry2D):
    def setUp(self):
        self.set_up('DensityWave', model_parameters='dims=2,B0=1.0',
                    magnetized=True)

class BlastMHDAntisymmetry2D(unittest.TestCase, TestAntisymmetry2D):
    def setUp(self):
        self.set_up('BlastMHD', usempi=True, np=8, resolution=32)

class BlastMHDRestartedToContinuous2D(unittest.TestCase,
                                        RestartedToContinuous):
    def setUp(self):
        self.set_up('BlastMHD', resolution=128)

class BlastMHDRestartedToContinuous3D(unittest.TestCase, RestartedToContinuous):
    def setUp(self):
        self.set_up('BlastMHD', resolution=16, model_parameters='three_d=true')

class BlastMHDSerialToParallel3D(unittest.TestCase, SerialToParallel):
    def setUp(self):
        self.set_up('BlastMHD', resolution=16, model_parameters='three_d=true')

class MagnetarWindTopDownSymmetry(unittest.TestCase, TestTopDownSymmetry2D):
    """
    This test currently fails because of inaccuracies in the quartic solver, so
    it passes when simple eigenvalues are used.
    """
    def setUp(self):
        self.set_up('MagnetarWind', usempi=False,
                    np=1, resolution=32, rkorder=1, tmax=0.5)

class MagnetarWindSerialToParallel2D(unittest.TestCase, SerialToParallel):
    def setUp(self):
        self.set_up('MagnetarWind', resolution=64, tmax=0.1,
                    model_parameters='three_d=false,B_wind=0.1')

class MagnetarWindRestartedToContinuous2D(unittest.TestCase, RestartedToContinuous):
    def setUp(self):
        self.set_up('MagnetarWind', resolution=64,
                    model_parameters='three_d=false,B_wind=0.1')

class MagnetarWindSerialToParallel3D(unittest.TestCase, SerialToParallel):
    def setUp(self):
        self.set_up('MagnetarWind', resolution=16, tmax=0.1,
                    model_parameters='three_d=true')

class MagnetarWindRestartedToContinuous3D(unittest.TestCase, RestartedToContinuous):
    def setUp(self):
        self.set_up('MagnetarWind', resolution=16, np=1,
                    model_parameters='three_d=true')



if __name__ == '__main__':
    unittest.main(verbosity=2)
