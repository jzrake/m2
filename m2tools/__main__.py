import command
import plot_command
import sim_runner
import to3d
import tovtk
import stellar_model
import mhd_forces
import reductions
import streamplot

try:
    import local_tools
except ImportError:
    pass

parser = command.create_parser()
args = parser.parse_args()
cmd = args.command
cmd.run(args)
