--------------------------------------------------------------------------------
--
--                            M2 Lua driver
--
--------------------------------------------------------------------------------

local colors  = require 'ansicolors'
local serpent = require 'serpent'
local lfs     = require 'lfs'
local struct  = require 'struct'
local class   = require 'class'
local m2lib   = require 'm2lib'
local hdf5    = require 'hdf5'
local array   = require 'array'
local logger  = require 'logger'


local function print_splash()
   print(([[
	    __  __       _                 
	    |  \/  |     | |___      _____  
	    | |\/| |_____| __\ \ /\ / / _ \ 
	    | |  | |_____| |_ \ V  V / (_) |
	    |_|  |_|      \__| \_/\_/ \___/ 

            astrophysical MHD code
            version: %s
            build date: %s
            author: Jonathan Zrake (2014)
            Stanford University, KIPAC
	 ]]):format(m2lib.M2_GIT_SHA, m2lib.M2_BUILD_DATE))
end
local function to_enum(s)
   local res = m2lib['M2_'..s:upper()]
   if not res then
      error(("m2 has no enum value '%s'"):format(s))
   else
      return res
   end
end


local m2Application = class.class 'm2Application'


------------------------------------------------
-- forward methods of m2Application to m2sim
------------------------------------------------
for _, method in pairs(struct.methods('m2sim')) do
   m2Application[method] = function(self, ...)
      return self._m2[method](self._m2, ...)
   end
end

------------------------------------------------
-- m2Application:[get,set]_parameter(value)
------------------------------------------------
for _, member in pairs(struct.members 'm2sim') do
   m2Application['set_'..member] = function(self, val)
      self._m2[member] = val
   end
end
for _, member in pairs(struct.members 'm2sim') do
   m2Application['get_'..member] = function(self)
      return self._m2[member]
   end
end

function m2Application:__init__(args)
   self._logger = logger.CommandLineLogger(class.classname(self))
   args = args or { }
   local resolution = args.resolution or {128,1,1}
   local lower = args.lower or {0,0,0}
   local upper = args.upper or {1,1,1}
   local scaling = {to_enum'linear', to_enum'linear', to_enum'linear'}
   local geometry = to_enum(args.geometry or 'cartesian')
   for i,v in ipairs(args['scaling'] or { }) do scaling[i] = to_enum(v) end
   self._m2 = m2lib.m2sim()
   self._m2.domain_resolution = m2lib.ivec4(0, unpack(resolution))
   self._m2.domain_extent_lower = m2lib.dvec4(0, unpack(lower))
   self._m2.domain_extent_upper = m2lib.dvec4(0, unpack(upper))
   self._m2.coordinate_scaling1 = scaling[1]
   self._m2.coordinate_scaling2 = scaling[2]
   self._m2.coordinate_scaling3 = scaling[3]
   self._m2.number_guard_zones0 = m2lib.ivec4(0,unpack(args.guard0 or {1,1,1}))
   self._m2.number_guard_zones1 = m2lib.ivec4(0,unpack(args.guard1 or {0,0,0}))
   self._m2.geometry = geometry
   self._m2.relativistic = args.relativistic and 1 or 0
   self._m2.magnetized = args.magnetized and 1 or 0
   self._m2:initialize()
   self.status = self._m2.status
end

function m2Application:set_problem(problem)
   self._problem = problem
end

function m2Application:print_splash()
   print_splash()
end

function m2Application:print_config()
   print()
   local items = { }
   for i,member in ipairs(struct.members(self._m2)) do
      items[member] = self._m2[member]
   end
   self:_log(serpent.block(items, {nocode=true}), 1)
   print()
end

function m2Application:print_status()
   print()
   local items = { }
   for i,member in ipairs(struct.members(self._m2.status)) do
      items[member] = self._m2.status[member]
   end
   self:_log(serpent.block(items, {nocode=true}), 1)
   print()
end

function m2Application:write_checkpoint_hdf5(fname)
   print()
   local fast_method = true
   local global_shape = self:global_shape()
   local h5file   = hdf5.File(fname, 'w')
   local h5prim   = hdf5.Group(h5file, 'prim')
   local h5status = hdf5.Group(h5file, 'status')
   local h5config = hdf5.Group(h5file, 'config')
   local wguard = { }
   local start = { }
   for i=1,#global_shape do
      start[i] = self._m2.number_guard_zones0[i]
      wguard[i] = global_shape[i] +
	 self._m2.number_guard_zones0[i] +
	 self._m2.number_guard_zones1[i]
   end

   if fast_method then
      for i,field in ipairs(struct.members('m2prim')) do
	 self:_log('writing '..fname..'/prim/'..field, 1)
	 local h5d = hdf5.DataSet(h5prim, field, 'w', {shape=global_shape})
	 local data = array.view(self:get_volume_data(field), 'double',
				 wguard, start, global_shape)
	 h5d[nil] = data
	 h5d:close()
      end
   else
      for _,field in ipairs(struct.members('m2prim')) do
	 self:_log('writing '..fname..'/prim/'..field, 1)
	 local h5d = hdf5.DataSet(h5prim, field, 'w', {shape=global_shape})
	 for n,v in self:volumes() do
	    if v.zone_type == m2lib.M2_ZONE_TYPE_FULL then
	       local I = {v.global_index[1],
			  v.global_index[2],
			  v.global_index[3]}
	       h5d[I] = v.prim[field]
	    end
	 end
	 h5d:close()
      end
   end

   for i,member in ipairs(struct.members(self._m2)) do
      h5config[member] = tostring(self._m2[member])
   end
   for i,member in ipairs(struct.members(self.status)) do
      h5status[member] = self.status[member]
   end
   if self._problem then
      self._problem:write_hdf5_problem_data(h5file)
   end
   h5file['version'] = 'm2: '..m2lib.M2_GIT_SHA
   h5file['build_date'] = m2lib.M2_BUILD_DATE
   h5file['time_stamp'] = os.date("%c")
   h5file:close()
   print()
end

function m2Application:global_shape()
   local N = self._m2.domain_resolution
   local S = { }
   if N[1] > 1 then S[#S+1] = N[1] end
   if N[2] > 1 then S[#S+1] = N[2] end
   if N[3] > 1 then S[#S+1] = N[3] end
   return S
end

function m2Application:_log(msg, indent)
   self._logger:log_message(debug.getinfo(2).name, msg, indent)
end


return {m2Application=m2Application,
	print_splash=print_splash,
	to_enum=to_enum}
