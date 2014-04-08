local argparse = require 'argparse'
local problems = require 'problems'
local m2app    = require 'm2app'
local logger   = require 'logger'


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
   :description 'astrophysical MHD code'
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
   parser:mutex(parser:flag '--relativistic',
		parser:flag '--newtonian')
   parser:mutex(parser:flag '--magnetized',
		parser:flag '--unmagnetized')
   parser:argument 'ProblemClass' :convert(problems)

   local args = parser:parse()
   local problem = args.ProblemClass()

   local function user_config_callback(runtime_cfg)
      for k,v in pairs(args) do
	 runtime_cfg[k:gsub('-','_')] = v
      end
      return args
   end

   if args.explain then
      problem:describe()
   else
      problem:update_model_parameters(args['model-parameters'])
      problem:run(user_config_callback)
   end
end

main()
