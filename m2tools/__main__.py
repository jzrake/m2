import command
import plot_command
import sim_runner
import to3d
import stellar_model

parser = command.create_parser()
args = parser.parse_args()
cmd = args.command
cmd.run(args)
