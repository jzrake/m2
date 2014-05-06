local serpent = require 'serpent'
local class   = require 'class'
local logger  = require 'logger'
local m2lib   = require 'm2lib'
local m2app   = require 'm2app'
local hdf5    = require 'hdf5'



--------------------------------------------------------------------------------
-- TestProblem base class
--------------------------------------------------------------------------------
local TestProblem = class.class 'TestProblem'
function TestProblem:write_hdf5_problem_data(h5file)
   h5file['problem_name'] = class.classname(self)
   local h5mp = hdf5.Group(h5file, 'model_parameters', 'w')
   for k,v in pairs(self.model_parameters or { }) do
      h5mp[k] = tostring(v)
   end
end
function TestProblem:read_hdf5_problem_data(h5file)
   local h5mp = hdf5.Group(h5file, 'model_parameters', 'r+')
   for k,v in pairs(self.model_parameters or { }) do
      self.model_parameters[k] = h5mp[k]:value()
   end
end
function TestProblem:set_runtime_defaults(cfg) end
function TestProblem:update_model_parameters(user_mp)
   local problem_mp = self.model_parameters or { }
   for k,v in pairs(user_mp or { }) do
      if problem_mp[k] ~= nil then
         problem_mp[k] = v
      else
         error(("problem class '%s' does not accept model parameter '%s'")
         :format(class.classname(self), k))
      end
   end
end
function TestProblem:describe()
   if self.explanation then
      print('\n'..self.explanation..'\n')
   end
   local mp = { }
   for n in pairs(self.model_parameters or { }) do table.insert(mp, n) end
   table.sort(mp)
   print("model parameters for problem '"..class.classname(self).."':\n")
   for _,k in ipairs(mp) do
      local v = self.model_parameters[k]
      local d = self.model_parameters_doc[k]
      print(("%-40s %s"):format(("%16s .......... %s"):format(k, v), d))
   end
   print()
end
function TestProblem:reconfigure_after_failure(m2, attempt)
   return false
end
function TestProblem:run(user_config_callback, restart_file)
   -----------------------------------------------------------------------------
   -- Build runtime_cfg table from with precedence given first to the
   -- user_config_callback (probably feeding in command line arguments), then
   -- the restart file if it exists, and finally the problem's
   -- set_runtime_defaults method.
   -----------------------------------------------------------------------------
   local runtime_cfg = {
      tmax = 0.4,
      output_path = 'data',
      msg_cadence = 1
   }
   self:set_runtime_defaults(runtime_cfg)
   if restart_file then
      local h5f = hdf5.File(restart_file, 'r')
      local err, restored_runtime_cfg = serpent.load(h5f['runtime_cfg']:value())
      for k,v in pairs(restored_runtime_cfg) do
         runtime_cfg[k] = v
      end
      h5f:close()
   end
   if user_config_callback then
      user_config_callback(runtime_cfg)
   end

   local log = logger.CommandLineLogger(class.classname(self))
   local m2 = self:build_m2(runtime_cfg)
   local dims = m2._cart_comm:get 'dims'

   m2:print_splash()
   m2:print_config()
   log:log_message('run', serpent.block(runtime_cfg), 1)
   log:log_message('run', ("domain decomposition: [%d %d %d]")
   :format(dims[1], dims[2], dims[3]), 1)

   if restart_file then
      m2:read_checkpoint_hdf5(restart_file)
      m2:run_initial_data()
   else
      m2:run_initial_data(self.initial_data_cell,
      self.initial_data_face,
      self.initial_data_edge)
   end
   self:describe()


   if runtime_cfg.vis then
      m2:visualize()
   end


   while m2.status.time_simulation < runtime_cfg.tmax do
      local cad = m2:get_cadence_checkpoint_hdf5()
      local now = m2.status.time_simulation
      local las = m2.status.time_last_checkpoint_hdf5
      local num = m2.status.checkpoint_number_hdf5

      if cad > 0.0 and now - las > cad then
         -- Print and then update the status
         num = num + 1
         m2:print_status()
         m2.status.checkpoint_number_hdf5 = num
         m2.status.time_last_checkpoint_hdf5 = las + cad

         -- Write data, config, status, problem data, and serialized runtime_cfg
         -- into the checkpoint file
         m2:write_checkpoint_hdf5(
         ('%s/chkpt.%04d.h5'):format(runtime_cfg.output_path, num),
         {runtime_cfg=runtime_cfg})
      end

      if m2.status.iteration_number > 0 and
         m2.status.iteration_number % runtime_cfg.msg_cadence == 0 then
         log:log_message('run', m2.status:get_message())
      end

      local keep_trying = true
      local cached_config = m2:get_config{
         -- only these values should be modified by the failure handler
         'cfl_parameter',
         'plm_parameter',
         'quartic_solver',
         'riemann_solver',
         'rk_order',
         'suppress_extrapolation_at_unhealthy_zones'
      }

      m2.status.drive_attempt = 0
      while keep_trying do
         m2:drive()
         if m2.status.error_code ~= 0 then
            m2.status.drive_attempt = m2.status.drive_attempt + 1
            keep_trying = self:reconfigure_after_failure(
            m2,
            m2.status.drive_attempt)
            print(m2.status.num_failures)
         else
            keep_trying = false
         end
      end

      if m2.status.error_code ~= 0 then
         log:log_message('run', 'exiting due to failures')
         break
      end

      if m2.status.drive_attempt > 0 then m2:update_config(cached_config) end

      collectgarbage()
   end

   m2:write_checkpoint_hdf5(
   ('%s/chkpt.final.h5'):format(runtime_cfg.output_path),
   {runtime_cfg=runtime_cfg})
end



local function outflow_bc_flux(axis)
   local flux = 'flux'..axis
   local nhat = m2lib.dvec4()
   nhat[axis] = 1.0
   local function bc(V0)
      if V0.global_index[axis] == -1 then
         local V1 = V0:neighbor(axis, 1)
         if V1 == nil then return end
         V0[flux] = V1.aux:fluxes(nhat)
      else
         V0[flux] = V0.aux:fluxes(nhat)
      end
  end
  return bc
end



local DensityWave = class.class('DensityWave', TestProblem)
DensityWave.explanation = [[
--------------------------------------------------------------------------------
-- Smooth sinusoidal density and/or pressure fluctuation
--
--------------------------------------------------------------------------------]]
function DensityWave:__init__()
   local mps = { }
   local doc = { }
   self.model_parameters = mps
   self.model_parameters_doc = doc

   -- ==========================================================================
   mps.dims = 1; doc.dims = 'domain dimensionality'
   mps.B0 = 0.0; doc.B0 = 'magnetic field strength'
   -- ==========================================================================

   local pi = math.pi
   local L = 1
   self.initial_data_cell = function(x)
      local dims = mps.dims
      local v1 = dims >= 1 and 1.0 or 0.0
      local v2 = dims >= 2 and 1.0 or 0.0
      local v3 = dims >= 3 and 1.0 or 0.0
      local r1 = dims >= 1 and x[1] or 0.0
      local r2 = dims >= 2 and x[2] or 0.0
      local r3 = dims >= 3 and x[3] or 0.0
      local r = r1 + r2 + r3
      return { 1.0 + 0.5 * math.sin(4*pi*r/L), 1.0, v1, v2, v3}
   end
   self.initial_data_face = function(x, n)
      return { mps.B0 * (n[1] + n[2] + n[3]) }
   end
end
function DensityWave:build_m2(runtime_cfg)
   local mps = self.model_parameters
   local dims = mps.dims
   local N1 = dims >= 1 and (runtime_cfg.resolution or 64) or 1
   local N2 = dims >= 2 and (runtime_cfg.resolution or 64) or 1
   local N3 = dims >= 3 and (runtime_cfg.resolution or 64) or 1
   local build_args = {lower={0.0, 0.0, 0.0},
   upper={1.0, 1.0, 1.0},
   resolution={N1,N2,N3},
   periods={true,true,true},
   scaling={'linear'},
   relativistic=false,
   magnetized=false,
   geometry='cartesian'}
   if runtime_cfg.relativistic then build_args.relativistic = true end
   if runtime_cfg.unmagnetized then build_args.magnetized = false end
   local m2 = m2app.m2Application(build_args)
   m2:set_cadence_checkpoint_hdf5(runtime_cfg.hdf5_cadence or 0.0)
   m2:set_cadence_checkpoint_tpl(runtime_cfg.tpl_cadence or 0.0)
   m2:set_gamma_law_index(5./3)
   m2:set_rk_order(runtime_cfg.rkorder or 2)
   m2:set_cfl_parameter(0.3)
   m2:set_plm_parameter(2.0)
   m2:set_interpolation_fields(m2lib.M2_PRIMITIVE)
   m2:set_riemann_solver(runtime_cfg.riemann_solver or m2lib.M2_RIEMANN_HLLE)
   m2:set_problem(self)
   return m2
end



--------------------------------------------------------------------------------
-- base class for shocktube probems
--------------------------------------------------------------------------------
local Shocktube = class.class('Shocktube', TestProblem)
function Shocktube:build_m2(runtime_cfg)
   local build_args = {lower={0.0, 0.0, 0.0},
   upper={1.0, 1.0, 1.0},
   resolution={512,1,1},
   periods={false,false,false},
   scaling={'linear'},
   relativistic=false,
   magnetized=true,
   geometry='cartesian'}
   if runtime_cfg.relativistic then build_args.relativistic = true end
   if runtime_cfg.unmagnetized then build_args.magnetized = false end
   if runtime_cfg.resolution then
      build_args.resolution[1] = runtime_cfg.resolution
   end
   local m2 = m2app.m2Application(build_args)
   m2:set_cadence_checkpoint_hdf5(runtime_cfg.hdf5_cadence or 0.0)
   m2:set_cadence_checkpoint_tpl(runtime_cfg.tpl_cadence or 0.0)
   m2:set_gamma_law_index(self.gamma_law_index or 5./3)
   m2:set_rk_order(runtime_cfg.rkorder or 2)
   m2:set_cfl_parameter(0.8)
   m2:set_plm_parameter(2.0)
   m2:set_interpolation_fields(m2lib.M2_PRIMITIVE)
   m2:set_riemann_solver(runtime_cfg.riemann_solver or m2lib.M2_RIEMANN_HLLE)
   m2:set_boundary_conditions_flux1(outflow_bc_flux(1))
   m2:set_problem(self)
   return m2
end



local BrioWu = class.class('BrioWu', Shocktube)
BrioWu.explanation = [[
--------------------------------------------------------------------------------
-- Brio-Wu test
--
-- http://www.astro.princeton.edu/~jstone/Athena/tests/brio-wu/Brio-Wu.html
--
-- From section 5 of: Brio, M. & C.C. Wu, "An Upwind Differencing Scheme for the
-- Equations of Ideal Magnetohydrodynamics", Journal of Computational Physics,
-- 75, 400-422 (1988).
--------------------------------------------------------------------------------]]
function BrioWu:__init__()
   self.gamma_law_index = 2.0
   self.initial_data_cell = function(x)
      if x[1] < 0.5 then
         return {1.000, 1.0, 0.0, 0.0, 0.0}
      else
         return {0.125, 0.1, 0.0, 0.0, 0.0}
      end
   end
   self.initial_data_face = function(x, n)
      if x[1] < 0.5 then
         return {0.75*n[1] + n[2]}
      else
         return {0.75*n[1] - n[2]}
      end
   end
end
function BrioWu:set_runtime_defaults(runtime_cfg)
   runtime_cfg.tmax = 0.15
end



local RyuJones = class.class('RyuJones', Shocktube)
RyuJones.explanation = [[
--------------------------------------------------------------------------------
-- Ryu and Jones Test 2A (Newtonian only)
--
-- http://www.astro.princeton.edu/~jstone/Athena/tests/rj2a/RJ2a.html
--
-- From section 4 of: Ryu, D. & Jones, T.W., "Numerical Magnetohydrodynamics in
-- Astrophysics: Algorithm and Tests for One-Dimensional Flow", Astro. J., 442,
-- 228-258 (1995).
--------------------------------------------------------------------------------]]
function RyuJones:set_runtime_defaults(runtime_cfg)
   runtime_cfg.tmax = 4.0
end
function RyuJones:__init__()
   self.initial_data_cell = function(x)
      if x[1] < 0.5 then
         return {1.08, 0.95, 1.2, 0.01, 0.5}
      else
         return {1.00, 1.00, 0.0, 0.00, 0.0}
      end
   end
   self.initial_data_face = function(x, n)
      if x[1] < 0.5 then
         return {2.0*n[1] + 4.0*n[2] + 2.0*n[3]}
      else
         return {2.0*n[1] + 3.6*n[2] + 2.0*n[3]}
      end
   end
end



local ContactWave = class.class('ContactWave', Shocktube)
ContactWave.explanation = [[
--------------------------------------------------------------------------------
-- The contact wave problem is meant for testing the HLLC family of Riemann
-- solvers. When HLLC solvers are working, solutions which vary only in the
-- density and tangential velocity will be numerically time-independent.
--------------------------------------------------------------------------------]]
function ContactWave:__init__()
   local mps = { }
   local doc = { }
   self.model_parameters = mps
   self.model_parameters_doc = doc

   -- ==========================================================================
   mps.vt = 0.0; doc.vt = 'transverse velocity magnitude'
   mps.Bt = 0.0; doc.Bt = 'transverse magnetic field strength'
   -- ==========================================================================

   self.initial_data_cell = function(x)
      local vt = mps.vt
      if x[1] < 0.5 then
         return {1.00, 1.00, 0.0, vt, vt}
      else
         return {0.10, 1.00, 0.0,-vt,-vt}
      end
   end
   self.initial_data_face = function(x, n)
      local Bt = mps.Bt
      if x[1] < 0.5 then
         return {1.0*n[1] + Bt*n[2] + Bt*n[3]}
      else
         return {1.0*n[1] - Bt*n[2] - Bt*n[3]}
      end
   end
end



local BlastMHD = class.class('BlastMHD', TestProblem)
BlastMHD.explanation = [[
--------------------------------------------------------------------------------
-- MHD blast test
--
-- http://www.astro.princeton.edu/~jstone/Athena/tests/blast/blast.html
--------------------------------------------------------------------------------]]
function BlastMHD:__init__()
   local mps = { }
   local doc = { }
   self.model_parameters = mps
   self.model_parameters_doc = doc

   -- ==========================================================================
   mps.three_d = false; doc.three_d = 'run in three dimensions'
   -- ==========================================================================

   self.initial_data_cell = function(x)
      local r = (x[1]^2 + x[2]^2 + x[3]^2)^0.5
      if r < 0.1 then
         return {1.0, 10.0, 0.0, 0.0, 0.0}
      else
         return {1.0, 0.10, 0.0, 0.0, 0.0}
      end
   end
   self.initial_data_face = function(x, n)
      local B = { 2^-0.5, 2^-0.5, 0.0 }
      return {B[1]*n[1] + B[2]*n[2] + B[3]*n[3]}
   end
end
function BlastMHD:set_runtime_defaults(runtime_cfg)
   runtime_cfg.tmax = 0.4
end
function BlastMHD:build_m2(runtime_cfg)
   local build_args = {lower={-0.5,-0.5,-0.5},
   upper={ 0.5, 0.5, 0.5},
   resolution={64,64,1},
   periods={true,true,true},
   scaling={'linear', 'linear', 'linear'},
   relativistic=false,
   magnetized=true,
   geometry='cartesian'}
   if runtime_cfg.relativistic then build_args.relativistic = true end
   if runtime_cfg.unmagnetized then build_args.magnetized = false end
   if runtime_cfg.resolution then
      build_args.resolution[1] = runtime_cfg.resolution
      build_args.resolution[2] = runtime_cfg.resolution
   end
   if self.model_parameters.three_d then
      build_args.resolution[3] = runtime_cfg.resolution or 64
   end
   local m2 = m2app.m2Application(build_args)
   m2:set_cadence_checkpoint_hdf5(runtime_cfg.hdf5_cadence or 0.1)
   m2:set_cadence_checkpoint_tpl(runtime_cfg.tpl_cadence or 0.0)
   m2:set_gamma_law_index(5./3)
   m2:set_rk_order(runtime_cfg.rkorder or 2)
   m2:set_cfl_parameter(0.4)
   m2:set_plm_parameter(2.0)
   m2:set_interpolation_fields(m2lib.M2_PRIMITIVE)
   m2:set_riemann_solver(runtime_cfg.riemann_solver or m2lib.M2_RIEMANN_HLLE)
   --m2:set_quartic_solver(m2lib.M2_QUARTIC_NONE)
   m2:set_problem(self)
   return m2
end



local Jet = class.class('Jet', TestProblem)
Jet.explanation = [[
--------------------------------------------------------------------------------
-- Stability of relativistic magnetized jets
--------------------------------------------------------------------------------]]
function Jet:__init__()
   local m2jet = require 'm2jet'
   local mps = m2jet.jet_model_parameters()
   local doc = { }
   self.model_parameters = mps
   self.model_parameters_doc = doc

   self.initial_data_cell = function(x)
      return {mps.ambient_density,
      mps.mabient_pressure,
      0.0, 0.0, 0.0}
   end
   self.initial_data_edge = function(x, n)
      local nu = mps.nu_parameter
      local r = x[1]
      local t = x[2]
      local R = r * math.sin(t)
      local P = 0.25 * r^nu * (1.0 - math.abs(math.cos(t)))
      local Af = P / (R + 0.01) -- A_phi
      return {-Af * n[3]}
   end
end
function Jet:set_runtime_defaults(runtime_cfg)
   runtime_cfg.tmax = 128.0
end
function Jet:build_m2(runtime_cfg)
   local r0 = 1.0
   local r1 = 100.0
   local N = runtime_cfg.resolution or 64
   local build_args = {lower={r0, 0.0, 0.0},
   upper={r1, math.pi/2, 2*math.pi},
   resolution={ },
   periods={false,false,true},
   scaling={'logarithmic', 'linear', 'linear'},
   relativistic=false,
   magnetized=true,
   geometry='spherical'}
   build_args.resolution[1] = N / 2 * math.floor(math.log10(r1/r0))
   build_args.resolution[2] = N / 2
   build_args.resolution[3] = 1
   local m2 = m2app.m2Application(build_args)
   m2:set_cadence_checkpoint_hdf5(runtime_cfg.hdf5_cadence or 4.0)
   m2:set_cadence_checkpoint_tpl(runtime_cfg.tpl_cadence or 0.0)
   m2:set_gamma_law_index(5./3)
   m2:set_rk_order(runtime_cfg.rkorder or 2)
   m2:set_cfl_parameter(0.4)
   m2:set_plm_parameter(1.8)
   m2:set_interpolation_fields(m2lib.M2_PRIMITIVE)
   m2:set_riemann_solver(m2lib.M2_RIEMANN_HLLE)
   --m2:set_add_physical_source_terms(m2lib.jet_add_physical_source_terms)
   m2:set_boundary_conditions_cell(m2lib.jet_boundary_conditions_cell)
   m2:set_problem(self)
   return m2
end



local MagnetarWind = class.class('MagnetarWind', TestProblem)
MagnetarWind.explanation = [[
--------------------------------------------------------------------------------
-- Relativistic toroidally magnetized wind injected through a spherical inner
-- boundary
--------------------------------------------------------------------------------]]
function MagnetarWind:__init__()
   local mps = { }
   local doc = { }
   self.model_parameters = mps
   self.model_parameters_doc = doc

   -- ==========================================================================
   mps.three_d = false; doc.three_d='run in three dimensions'
   mps.r_outer=100; doc.r_outer='outer radius (inner is 1.0)'
   mps.r_cavity=10; doc.r_cavity='cavity radius'
   mps.d_star=100;  doc.d_star='star density (wind density is 1.0)'
   mps.B_wind=8.00; doc.B_wind='wind toroidal field value'
   mps.g_wind=8.00; doc.g_wind='wind Lorentz factor'
   mps.g_pert=0.00; doc.g_pert='fractional azimuthal variation in g_wind'
   -- ==========================================================================

   self.initial_data_cell = function(x)
      local d0
      if x[1] < mps.r_cavity  then
         d0 = 0.1
      else
         d0 = mps.d_star
      end
      return {d0, 0.01, 0.0, 0.0, 0.0}
   end
end
function MagnetarWind:set_runtime_defaults(runtime_cfg)
   runtime_cfg.tmax = 2.0
end
function MagnetarWind:build_m2(runtime_cfg)
   local build_args = {lower={  1.0, 0.0, 0.0},
   upper={100.0, math.pi, 2*math.pi},
   resolution={128,64,1},
   periods={false,false,true},
   scaling={'logarithmic', 'linear', 'linear'},
   relativistic=true,
   magnetized=true,
   geometry='spherical'}
   if runtime_cfg.newtonian then build_args.relativistic = false end
   if runtime_cfg.unmagnetized then build_args.magnetized = false end
   local mps = self.model_parameters
   local N = runtime_cfg.resolution or 64
   local r0 = 1.0
   local r1 = mps.r_outer
   build_args.upper[1] = r1
   build_args.resolution[1] = N * math.floor(math.log10(r1/r0))
   build_args.resolution[2] = N
   build_args.resolution[3] = N * 2
   if not mps.three_d then
      build_args.resolution[3] = 1
   end
   local function wind_inner_boundary_flux(V0)
      local nhat = m2lib.dvec4(0,1,0,0)
      local t = 0.5 * (V0.x0[2] + V0.x1[2])
      local p = 0.5 * (V0.x0[3] + V0.x1[3])
      if V0.global_index[1] == -1 then
         local d = 1.0
         local B = mps.B_wind
         local g = mps.g_wind
         local h = mps.g_pert
         g = g * (1.0 + h * math.sin(6*p)^2)
         V0.prim.v1 =(1.0 - g^-2)^0.5
         V0.prim.v2 = 0.0
         V0.prim.v3 = 0.0
         V0.prim.B1 = 0.0
         V0.prim.B2 = 0.0
         V0.prim.B3 = B * math.sin(t)
         V0.prim.d = d
         V0.prim.p = d * 0.01
         V0:from_primitive()
         V0.flux1 = V0.aux:fluxes(nhat)
      else
         V0.flux1 = V0.aux:fluxes(nhat)
      end
   end
   local m2 = m2app.m2Application(build_args)
   m2:set_cadence_checkpoint_hdf5(runtime_cfg.hdf5_cadence or 4.0)
   m2:set_cadence_checkpoint_tpl(runtime_cfg.tpl_cadence or 0.0)
   m2:set_gamma_law_index(4./3)
   m2:set_rk_order(runtime_cfg.rkorder or 2)
   m2:set_cfl_parameter(runtime_cfg.cfl_parameter or 0.3)
   m2:set_plm_parameter(runtime_cfg.plm_parameter or 1.5)
   m2:set_interpolation_fields(m2lib.M2_PRIMITIVE_AND_FOUR_VELOCITY)
   m2:set_riemann_solver(m2lib.M2_RIEMANN_HLLE)
   m2:set_boundary_conditions_flux1(wind_inner_boundary_flux)
   m2:set_boundary_conditions_flux2(outflow_bc_flux(2))
   m2:set_quartic_solver(runtime_cfg.quartic or m2lib.M2_QUARTIC_FULL_COMPLEX)
   m2:set_suppress_extrapolation_at_unhealthy_zones(1)
   m2:set_pressure_floor(runtime_cfg.pressure_floor or 1e-8)
   m2:set_problem(self)
   return m2
end
function MagnetarWind:reconfigure_after_failure(m2, attempt)
   local responses = {
     function()
         print('reduce cfl -> 0.1', attempt)
         print('reduce plm -> 1.0', attempt)
         m2:set_cfl_parameter(0.1)
         m2:set_plm_parameter(1.0)
         m2:set_quartic_solver(m2lib.M2_QUARTIC_NONE)
      end,
      function()
         print('reduce plm -> 0.0', attempt)
         m2:set_plm_parameter(0.0)
      end,
   }
   if responses[attempt] then
      responses[attempt]()
      return true
   else
      return false
   end
end


return {
   DensityWave  = DensityWave,
   RyuJones     = RyuJones,
   BrioWu       = BrioWu,
   ContactWave  = ContactWave,
   BlastMHD     = BlastMHD,
   Jet          = Jet,
   MagnetarWind = MagnetarWind
}
