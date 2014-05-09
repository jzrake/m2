import command
import checkpoint


class ThreeDimensionalizer(command.Command):
    _alias = "to3d"

    def configure_parser(self, parser):
        parser.add_argument('input_checkpoint')
        parser.add_argument('output_checkpoint')
        parser.add_argument('--density-perturbation', type=float, default=None)

    def run(self, args):
        c2d = checkpoint.Checkpoint(args.input_checkpoint)
        c3d = checkpoint.Checkpoint(args.output_checkpoint, writable=True)

        if c2d.dimensionality != 2 or c2d.domain_resolution[3] != 1:
            raise TypeError("input checkpoint must be 2D with symmetry along "
                            "the last axis")
        if c3d.dimensionality != 3:
            raise TypeError("output checkpoint must be 3D")

        if (c2d.domain_resolution[1:3] != c3d.domain_resolution[1:3]).all():
            raise TypeError("domain resolution must agree along first two axes")

        if c2d.problem_name != c3d.problem_name:
            raise TypeError("both checkpoints must run the same problem")

        for k in c2d.cell_primitive:
            print "writing {0}/{1}/{2}".format(
                c3d.filename, 'cell_primitive', k)
            indata = c2d.cell_primitive[k][:][:,:,None]
            if args.density_perturbation and k == 'd':
                from numpy.random import random
                A = args.density_perturbation
                pert = A * (random(c3d.cell_primitive[k].shape) - 0.5)
            else:
                pert = 0.0
            c3d.cell_primitive[k][:] = indata * (1.0 + pert)

        for k in c2d.face_magnetic_flux:
            print "writing {0}/{1}/{2}".format(
                c3d.filename, 'face_magnetic_flux', k)
            indata = c2d.face_magnetic_flux[k][:][:,:,None]
            if k in [1,2]:
                indata /= c3d.domain_resolution[3]
            c3d.face_magnetic_flux[k][:] = indata

        for k in c2d.status:
            print "writing {0}/{1}/{2}".format(c3d.filename, 'status', k)
            c3d.status[k] = c2d.status[k]

