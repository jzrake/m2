import command
import polar_plot
import shocktube_plot


parser = command.create_parser()
args = parser.parse_args()


cmd = args.command

cmd.run(args)
