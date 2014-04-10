local argparse = require 'argparse'
local problems = require 'problems'
local m2app    = require 'm2app'
local logger   = require 'logger'
local hdf5     = require 'hdf5'


local function main()
   -----------------------------------------------------------------------------
   -- Load system options from .m2rc file
   -----------------------------------------------------------------------------
   local err, m2rc = pcall(function() return loadfile '.m2rc' end)
   if m2rc then
      local rc = { }
      rc.CommandLineLogger = logger.CommandLineLogger._rc
      m2rc(rc)
   end

   -----------------------------------------------------------------------------
   -- Configure command line parser
   -----------------------------------------------------------------------------
   local parser = argparse()
   :name 'm2'
   parser:flag '-v' '--version'
   :action(
      function()
	 m2app.print_splash()
	 os.exit()
      end)
   parser:flag '--vis'
   parser:flag '-e' '--explain'
   parser:option '--tmax' :convert(tonumber)
   parser:option '--output-path' :default '.'
   parser:option '-N' '--resolution' :convert(tonumber)
   parser:option '--rkorder' :convert(tonumber)
   parser:option '-m' '--model-parameters' :convert(
      function(mp)
	 return load('return {'..mp..'}')()
      end)
   parser:option '--riemann-solver' :convert(
      function(rs)
	 return m2app.to_enum('riemann_'..rs)
      end)
   parser:option '--restart'
   parser:mutex(parser:flag '--relativistic',
		parser:flag '--newtonian')
   parser:mutex(parser:flag '--magnetized',
		parser:flag '--unmagnetized')
   parser:argument 'ProblemClass' :convert(problems) :args '?'

   local args = parser:parse()

   if args.restart then
      local h5f = hdf5.File(args.restart, 'r')
      args.ProblemClass = problems[h5f['problem_name']:value()]
      h5f:close()
   end

   if not args.ProblemClass then
      local choices = { }
      for k,v in pairs(problems) do table.insert(choices, k) end
      print(parser:get_usage())
      print("\navailable problems:\n\t"..table.concat(choices, '\n\t'))
      return
   end

   local problem = args.ProblemClass()

   local function user_config_callback(runtime_cfg)
      for k,v in pairs(args) do
	 runtime_cfg[k:gsub('-','_')] = v
      end
   end

   if args.explain then
      problem:describe()
   else
      problem:update_model_parameters(args['model-parameters'])
      problem:run(user_config_callback, args['restart'])
   end
end

main()
