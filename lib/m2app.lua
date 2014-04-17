--------------------------------------------------------------------------------
--
--                            M2 Lua driver
--
--------------------------------------------------------------------------------

local colors   = require 'ansicolors'
local serpent  = require 'serpent'
local lfs      = require 'lfs'
local struct   = require 'struct'
local class    = require 'class'
local m2lib    = require 'm2lib'
local hdf5     = require 'hdf5'
local buffer   = require 'buffer'
local array    = require 'array'
local logger   = require 'logger'
local parallel = require 'parallel'

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
   args = args or { }

   local resolution = args.resolution or {128,1,1}
   local lower = args.lower or {0,0,0}
   local upper = args.upper or {1,1,1}
   local scaling = {to_enum'linear', to_enum'linear', to_enum'linear'}
   local geometry = to_enum(args.geometry or 'cartesian')
   local proc_sizes = {0, 0, 0}
   local periods = {0, 0, 0}
   local Ng0 = m2lib.ivec4()
   local Ng1 = m2lib.ivec4()
   for i,v in ipairs(args['scaling'] or { }) do scaling[i] = to_enum(v) end

   for n=1,3 do
      if (args.periods or { })[n] then periods[n] = 1 end
      if resolution[n] == 1 then proc_sizes[n] = 1 end
   end

   self._cart_comm = parallel.MPI_CartesianCommunicator(3, proc_sizes, periods)
   self._logger = logger.CommandLineLogger(class.classname(self))
   self._m2 = m2lib.m2sim()
   local start, size = self._cart_comm:partition(resolution)

   for n=1,3 do
      if resolution[n] > 1 then
	 Ng0[n] = 2
	 Ng1[n] = 2
	 if periods[n] == 0 then
	    local coords = self._cart_comm:get 'coords'
	    local dims = self._cart_comm:get 'dims'
	    if coords[n] == 0 then
	       Ng0[n] = 1 -- so that there's a layer of shells
	    end
	    if coords[n] == dims[n] - 1 then
	       Ng1[n] = 0
	    end
	 end
      end
      self._m2.local_grid_size [n] = size [n] + Ng0[n] + Ng1[n]
      self._m2.local_grid_start[n] = start[n] - Ng0[n]
      self._m2.periodic_dimension[n] = periods[n]
   end
   self._m2.cart_comm = self._cart_comm._comm
   self._m2.domain_resolution = m2lib.ivec4(0, unpack(resolution))
   self._m2.domain_extent_lower = m2lib.dvec4(0, unpack(lower))
   self._m2.domain_extent_upper = m2lib.dvec4(0, unpack(upper))
   self._m2.number_guard_zones0 = Ng0
   self._m2.number_guard_zones1 = Ng1
   self._m2.coordinate_scaling1 = scaling[1]
   self._m2.coordinate_scaling2 = scaling[2]
   self._m2.coordinate_scaling3 = scaling[3]
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
   self:_log(serpent.block(struct.items(self._m2), {nocode=true}), 1)
   print()
end

function m2Application:print_status()
   print()
   self:_log(serpent.block(struct.items(self._m2.status), {nocode=true}), 1)
   print()
end

function m2Application:memory_selections()
   -----------------------------------------------------------------------------
   -- Build descriptions of file and memory space selections
   -----------------------------------------------------------------------------
   local global_shape = self:global_shape()
   local coords = self._cart_comm:get 'coords'
   local period = self._cart_comm:get 'periods'
   local Ng0 = self._m2.number_guard_zones0
   local Ng1 = self._m2.number_guard_zones1
   local file = { exten={}, start={}, strid={}, count={}, block={} }
   local mems = { exten={}, start={}, strid={}, count={}, block={} }

   for n=1,#global_shape do
      file.exten[n] = self._m2.domain_resolution[n]
      file.start[n] = self._m2.local_grid_start[n] + Ng0[n]
      file.strid[n] = 1
      file.count[n] = self._m2.local_grid_size[n] - Ng0[n] - Ng1[n]
      file.block[n] = 1

      mems.exten[n] = self._m2.local_grid_size[n]
      mems.start[n] = Ng0[n]
      mems.strid[n] = 1
      mems.count[n] = self._m2.local_grid_size[n] - Ng0[n] - Ng1[n]
      mems.block[n] = 1

      -- If there is no periodicity along this axis then add another layer of
      -- data at the left boundary for face-centered data
      if period[n] == 0 then
	 if coords[n] == 0 then
	    file.exten[n] = file.exten[n] + 1
	    file.count[n] = file.count[n] + 1
	    mems.exten[n] = mems.exten[n] + 1
	    mems.start[n] = mems.start[n] - 1
	    mems.count[n] = mems.count[n] + 1
	 else
	    file.exten[n] = file.exten[n] + 1
	    file.start[n] = file.start[n] + 1
	 end
      end
   end

   for _,s in pairs {file, mems} do
      s.space = hdf5.DataSpace()
      s.space:set_extent(s.exten)
      s.space:select_hyperslab(s.start, s.strid, s.count, s.block)
   end

   return file, mems
end

function m2Application:write_checkpoint_hdf5(fname, extras)
   print()
   local file, mems = self:memory_selections()

   if self._cart_comm:rank() == 0 then
      local h5file = hdf5.File(fname, 'w')
      local h5prim = hdf5.Group(h5file, 'prim')
      local h5face = hdf5.Group(h5file, 'face_magnetic_flux')
      for _,field in ipairs(struct.members('m2prim')) do
	 local h5d = hdf5.DataSet(h5prim, field, 'w', {shape=file.exten})
	 h5d:close()
      end
      for field=1,3 do
	 local h5d = hdf5.DataSet(h5face, field, 'w', {shape=file.exten})
	 h5d:close()
      end
      h5file:close()
   end


   -----------------------------------------------------------------------------
   -- Write cell-centered primitive and face-centered magnetic flux to HDF5
   -----------------------------------------------------------------------------
   local function write()
      local h5file = hdf5.File(fname, 'r+')
      local h5prim = hdf5.Group(h5file, 'prim')
      local h5face = hdf5.Group(h5file, 'face_magnetic_flux')

      for _,field in ipairs(struct.members('m2prim')) do
	 self:_log('writing '..h5prim:path()..'/'..field, 1)
	 local h5d = hdf5.DataSet(h5prim, field, 'r+', {shape=file.exten})
	 local data = self:get_volume_data(field)
	 h5d:write(data, mems.space, file.space)
	 h5d:close()
      end
      for field=1,3 do
	 self:_log('writing '..h5face:path()..'/'..field, 1)
	 local h5d = hdf5.DataSet(h5face, field, 'r+', {shape=file.exten})
	 local data = self:get_face_data(field)
	 h5d:write(data, mems.space, file.space)
	 h5d:close()
      end

      h5file:close()
   end
   self._cart_comm:call_in_serial(write)


   -----------------------------------------------------------------------------
   -- Write m2 configuration, status, problem, and build meta-data to HDF5
   -----------------------------------------------------------------------------
   if self._cart_comm:rank() == 0 then
      local h5file = hdf5.File(fname, 'r+')
      local h5status = hdf5.Group(h5file, 'status')
      local h5config = hdf5.Group(h5file, 'config')

      for i,member in ipairs(struct.members(self._m2)) do
	 h5config[member] = tostring(self._m2[member])
      end
      for i,member in ipairs(struct.members(self.status)) do
	 h5status[member] = self.status[member]
      end
      if self._problem then
	 self._problem:write_hdf5_problem_data(h5file)
      end
      for k,v in pairs(extras or { }) do
	 h5file[k] = serpent.block(v)
      end
      h5file['version'] = 'm2: '..m2lib.M2_GIT_SHA
      h5file['build_date'] = m2lib.M2_BUILD_DATE
      h5file['time_stamp'] = os.date("%c")
      h5file:close()
   end
   print()
end

function m2Application:read_checkpoint_hdf5(fname)
   print()
   local file, mems = self:memory_selections()

   local h5file   = hdf5.File(fname, 'r')
   local h5prim   = hdf5.Group(h5file, 'prim')
   local h5face   = hdf5.Group(h5file, 'face_magnetic_flux')
   local h5status = hdf5.Group(h5file, 'status')

   local function prod(A)
      local p = 1
      for _,a in ipairs(A) do p = p * a end
      return p
   end


   -----------------------------------------------------------------------------
   -- Read cell-centered primitive and face-centered magnetic flux from HDF5
   -----------------------------------------------------------------------------
   for _,field in ipairs(struct.members('m2prim')) do
      self:_log('reading '..h5prim:path()..'/'..field, 1)
      local N = prod(mems.space:get_extent())
      local buf = buffer.new_buffer(N * buffer.sizeof(buffer.double))
      h5prim[field]:read(buf, mems.space, file.space)
      self:set_volume_data(field, buf)
   end
   for field=1,3 do
      self:_log('reading '..h5face:path()..'/'..field, 1)
      local N = prod(mems_face.space:get_extent())
      local buf = buffer.new_buffer(N * buffer.sizeof(buffer.double))
      h5face[field]:read(buf, mems.space, file.space)
      self:set_face_data(field, buf)
   end


   -----------------------------------------------------------------------------
   -- Read status and problem meta-data from HDF5
   -----------------------------------------------------------------------------
   for i,member in ipairs(struct.members(self.status)) do
      self.status[member] = h5status[member]:value()
   end
   if self._problem then
      self._problem:read_hdf5_problem_data(h5file)
   end

   h5file:close()
   self:print_status()
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
