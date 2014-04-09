local serpent = require 'serpent'
local struct = require 'struct'
local colors = require 'ansicolors'
local m2lib = require 'm2lib'

local function test1()
   local f = m2lib.dvec4(1.1, 2.2, 3.3, 4.4)
   assert(f[0] == 1.1)
   assert(f[1] == 2.2)
   assert(f[2] == 3.3)
   assert(f[3] == 4.4)
end

local function test2()
   local m2sim = m2lib.m2sim()
   local I = m2lib.ivec4()
   local J = m2sim.domain_resolution
   local K = m2sim.domain_resolution
   local M = m2sim.domain_resolution
   I[0] = 130
   m2sim.domain_resolution = I
   assert(m2sim.domain_resolution[0] == 130)
end

local function test3()
   local m2sim = m2lib.m2sim()
   local S = m2sim.status
   local instance_table, obj = struct.debug(S)
   assert(instance_table[obj].__parent__ == m2sim)
   collectgarbage() -- leave S, gc, make sure it's not dead
   assert(instance_table[obj] ~= nil)
   S = nil -- delete S, gc, make sure it's dead
   collectgarbage()
   local instance_table, obj = struct.debug(m2sim)
   assert(instance_table[obj] ~= nil)
   m2sim = nil
   collectgarbage()
   assert(instance_table[obj] == nil)
end

local function test4()
   local m2sim = m2lib.m2sim()
   m2sim.domain_resolution = m2lib.ivec4(0, 16, 16, 16)
   m2sim:initialize()
   for i,v in m2sim:volumes() do
      assert(i == v.global_index[0])
   end
end

local function test5()
   local m2sim = m2lib.m2sim()
   for k,v in pairs(m2sim) do
      print(k, v)
   end
end

test1()
test2()
test3()
test4()
test5()

--m2lib.self_test()

