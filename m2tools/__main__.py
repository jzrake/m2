import command
import plot_command

parser = command.create_parser()
args = parser.parse_args()
cmd = args.command
cmd.run(args)
