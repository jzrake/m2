local serpent = require 'serpent'
local buffer = require 'buffer'
local class = require 'class'
local MPI = require 'MPI'


local function integer(size)
   return buffer.new_buffer(buffer.sizeof(buffer.int) * (size or 1))
end
local function ival(i, n)
   return buffer.get_typed(i, buffer.int, n or 0)
end


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
   for n=0, self:size()-1 do
      if n == self:rank() then
	 f(...)
      end
      self:barrier()
   end
end

function MPI_Communicator:print(rank, ...)
   rank = rank or 0
   if rank == 'all' then
      for n=0, self:size()-1 do
	 if n == self:rank() then
	    io.write(('[rank %d] '):format(n))
	    io.write(...)
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

local MPI_CartesianCommunicator = class.class('MPI_CartesianCommunicator',
					      MPI_Communicator)

function MPI_CartesianCommunicator:__init__(ndim)
   local reorder = 1
   local wrapped_ax = integer(ndim)
   local proc_sizes = integer(ndim)
   local pcomm_size = integer()
   local parent = MPI.COMM_WORLD

   for n=0,ndim-1 do
      buffer.set_typed(wrapped_ax, buffer.int, n, 1)
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
   for n=1, self._ndim do
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
   for n=1, self._ndim do
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
   for n=1, #coords do
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

local function test1()
   local comm = MPI_Communicator()
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
   parallel.MPI_Communicator = MPI_Communicator
   parallel.MPI_CartesianCommunicator = MPI_CartesianCommunicator
   return parallel
else
   MPI.Init()
   test2()
   collectgarbage()
   MPI.Finalize()
   print(debug.getinfo(1).source, ": All tests passed")
end
