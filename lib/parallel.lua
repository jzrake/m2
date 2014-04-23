local serpent = require 'serpent'
local buffer  = require 'buffer'
local class   = require 'class'
local MPI     = require 'MPI'


local function integer(size)
   return buffer.new_buffer(buffer.sizeof(buffer.int) * (size or 1))
end
local function ival(i, n)
   return buffer.get_typed(i, buffer.int, n or 0)
end
local function initialized()
   local a = integer(1)
   MPI.Initialized(a)
   return ival(a) == 1
end


-- =============================================================================
--
-- MPI_Communicator class
--
-- =============================================================================
local MPI_Communicator = class.class 'MPI_Communicator'
function MPI_Communicator:__init__()
   self._comm = MPI.COMM_WORLD
end
function MPI_Communicator:barrier()
   MPI.Barrier(self._comm)
end
function MPI_Communicator:rank()
   local rank = integer()
   MPI.Comm_rank(self._comm, rank)
   return ival(rank)
end
function MPI_Communicator:size()
   local size = integer()
   MPI.Comm_size(self._comm, size)
   return ival(size)
end
function MPI_Communicator:send(rank, message, tag)
   local serialized = serpent.block(message)
   local buf = buffer.new_buffer(serialized)
   MPI.Send(buf, #buf, MPI.BYTE, rank, tag or 0, self._comm)
end
function MPI_Communicator:recv(rank, tag)
   local status = MPI.Status()
   local recv_size = integer()
   MPI.Probe(rank, tag or 0, self._comm, status)
   MPI.Get_count(status, MPI.BYTE, recv_size)
   local buf = buffer.new_buffer(ival(recv_size))
   MPI.Recv(buf, #buf, MPI.BYTE, rank, tag or 0, self._comm, status)
   local err, val = serpent.load(tostring(buf))
   return val
end
function MPI_Communicator:call_in_serial(f, ...)
   for n=0,self:size()-1 do
      if n == self:rank() then
	 f(...)
      end
      self:barrier()
   end
end
function MPI_Communicator:print(rank, ...)
   rank = rank or 0
   if rank == 'all' then
      for n=0,self:size()-1 do
	 if n == self:rank() then
	    io.write(('[rank %d] '):format(n))
	    for _,s in ipairs {...} do
	       io.write(tostring(s)..'  ')
	    end
	    io.write('\n')
	 end
	 self:barrier()
      end
   else
      if self:rank() == rank then
	 io.write(...)
	 io.write('\n')
      end
   end
end



-- =============================================================================
--
-- MPI_CartesianCommunicator class
--
-- =============================================================================
local MPI_CartesianCommunicator = class.class('MPI_CartesianCommunicator',
					      MPI_Communicator)
function MPI_CartesianCommunicator:__init__(ndim, sizes, periods)
   local reorder = 1
   local wrapped_ax = integer(ndim)
   local proc_sizes = integer(ndim)
   local pcomm_size = integer()
   local parent = MPI.COMM_WORLD
   sizes = sizes or { }
   periods = periods or { }
   for n=0,ndim-1 do
      buffer.set_typed(proc_sizes, buffer.int, n, sizes[n+1] or 0)
      buffer.set_typed(wrapped_ax, buffer.int, n, periods[n+1] or 1)
   end
   self._ndim = ndim
   self._comm = MPI.Comm()
   MPI.Comm_size(parent, pcomm_size)
   MPI.Dims_create(ival(pcomm_size), ndim, proc_sizes)
   MPI.Cart_create(parent, ndim,
		   proc_sizes,
		   wrapped_ax, reorder, self._comm)
end
function MPI_CartesianCommunicator:__gc__()
   MPI.Comm_free(self._comm)
end
function MPI_CartesianCommunicator:cart_coords(rank)
   local proc_index = integer(self._ndim)
   local ret = { }
   MPI.Cart_coords(self._comm, rank, self._ndim, proc_index)
   for n=1,self._ndim do
      ret[n] = buffer.get_typed(proc_index, buffer.int, n-1)
   end
   return ret
end
function MPI_CartesianCommunicator:get(what)
   local dims = integer(self._ndim)
   local periods = integer(self._ndim)
   local coords = integer(self._ndim)
   local ret = { dims={ }, periods={ }, coords={ } }
   MPI.Cart_get(self._comm, self._ndim, dims, periods, coords)
   for n=1,self._ndim do
      ret.dims   [n] = buffer.get_typed(dims,    buffer.int, n-1)
      ret.periods[n] = buffer.get_typed(periods, buffer.int, n-1)
      ret.coords [n] = buffer.get_typed(coords,  buffer.int, n-1)
   end
   if what then
      return ret[what]
   else
      return ret
   end
end
function MPI_CartesianCommunicator:partition(shape)
   local dims = self:get 'dims'
   local coords = self:get 'coords'
   local start = { }
   local size = { }
   for n=1,#coords do
      local R = shape[n] % dims[n]
      local S = math.floor(shape[n] / dims[n])
      if coords[n] < R then
	 size[n] = S + 1
      else
	 size[n] = S
      end
      start[n] = 0
      for j=0,coords[n]-1 do
	 if j < R then
	    start[n] = start[n] + S + 1
	 else
	    start[n] = start[n] + S
	 end
      end
   end
   return start, size
end



-- =============================================================================
--
-- MPI_Communicator_S class
--
-- =============================================================================
local MPI_Communicator_S = class.class 'MPI_Communicator_S'
function MPI_Communicator_S:__init__() end
function MPI_Communicator_S:barrier() end
function MPI_Communicator_S:rank() return 0 end
function MPI_Communicator_S:size() return 1 end
function MPI_Communicator_S:send(rank, message, tag) end
function MPI_Communicator_S:recv(rank, tag) end
function MPI_Communicator_S:call_in_serial(f, ...) f(...) end
function MPI_Communicator_S:print(rank, ...)
   local n = rank
   io.write(('[rank %d] '):format(n))
   for _,s in ipairs {...} do
      io.write(tostring(s)..'  ')
   end
   io.write('\n')
end



-- =============================================================================
--
-- MPI_CartesianCommunicator_S class
--
-- =============================================================================
local MPI_CartesianCommunicator_S = class.class('MPI_CartesianCommunicator_S',
						MPI_Communicator_S)
function MPI_CartesianCommunicator_S:__init__(ndim, sizes, periods)
   self._ndim = ndim
   self._sizes = sizes
   self._periods = periods
end
function MPI_CartesianCommunicator_S:__gc__() end
function MPI_CartesianCommunicator_S:cart_coords(rank)
   local ret = { }
   for i=n,rank do ret[n] = 0 end
   return ret
end
function MPI_CartesianCommunicator_S:get(what)
   local dims = integer(self._ndim)
   local periods = integer(self._ndim)
   local coords = integer(self._ndim)
   local ret = { dims={ }, periods={ }, coords={ } }
   for n=1,self._ndim do
      ret.dims   [n] = 1
      ret.periods[n] = self._periods[n]
      ret.coords [n] = 0
   end
   if what then
      return ret[what]
   else
      return ret
   end
end
function MPI_CartesianCommunicator_S:partition(shape)
   local dims = self:get 'dims'
   local coords = self:get 'coords'
   local start = { }
   local size = { }
   return MPI_CartesianCommunicator.partition(self, shape)
end



local function test1()
   local comm = MPI_Communicator()
   if comm:size() ~= 2 then
      print "skipping test1, must be run on exactly 2 processes"
      return
   end
   if comm:rank() == 0 then
      comm:send(1, {"here is the message", { }})
   elseif comm:rank() == 1 then
      local msg = comm:recv(0)
      print(msg[1], msg[2])
   end
end



local function test2()
   local comm = MPI_CartesianCommunicator(3)
   for n=0, comm:size()-1 do
      if n == comm:rank() then
	 print(comm:rank(), unpack(comm:get 'coords'))
      end
      comm:barrier()
   end
end



if ... then -- if __name__ == "__main__"
   local parallel = { }
   if next(MPI) then
      if initialized() then
	 parallel.MPI_Communicator = MPI_Communicator
	 parallel.MPI_CartesianCommunicator = MPI_CartesianCommunicator
      end
   end
   if not next(parallel) then
      parallel.MPI_Communicator = MPI_Communicator_S
      parallel.MPI_CartesianCommunicator = MPI_CartesianCommunicator_S
   end
   parallel.initialized = initialized
   return parallel
else
   test1()
   test2()
   print(debug.getinfo(1).source, ": All tests passed")
end
