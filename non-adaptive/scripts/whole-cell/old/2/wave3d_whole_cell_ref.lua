-- load pre-implemented lua functions
ug_load_script("ug_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("plugins/Limex/limex_util.lua")

--AssertPluginsLoaded({"neuro_collection", "Limex", "Parmetis"})
--AssertPluginsLoaded({"neuro_collection", "Limex"})

InitUG(3, AlgebraType("CPU", 1))

--EnableLUA2C(true)  -- speed up evaluation of lua functions by c program

-------------------------------------
-- parse command line parameters  ---
-------------------------------------

-- num ERM refs
numERMRefs = util.GetParamNumber("-numERMRefs", 0)

-- choice of grid name
gridName = util.GetParam("-grid", "aDendrite.ugx")

-- vtk folder
vtkFolder = util.GetParam("-vtkFolder", "vtk/")


-- refinements (global and at ERM)
--numGlobRefs = util.GetParamNumber("-numGlobRefs", 0)
numGlobRefs = util.GetParamNumber("-numPostRefs", 0)


-- which ER mechanisms are to be activated?
setting = util.GetParam("-setting", "all")
setting = string.lower(setting)
validSettings = {}
validSettings["all"] = 0;
validSettings["none"] = 0;
validSettings["ip3r"] = 0;
validSettings["ryr"] = 0;
if (validSettings[setting] == nil) then
    error("Unknown setting " .. setting)
end

-- densities
ryrDensity = util.GetParamNumber("-ryrDensity", 0.86) -- was 3.0 or 0.86?
pmcaDensity = util.GetParamNumber("-pmcaDensity", 500)
ncxDensity = util.GetParamNumber("-ncxDensity", 15)

-- buffer
totalBuffer = util.GetParamNumber("-totBuf", 4*40.0e-6)


-- choice of algebra
-- useBlockAlgebra = util.HasParamOption("-block")


-- choice of solver setup
solverID = util.GetParam("-solver", "GMG")
solverID = string.upper(solverID)
validSolverIDs = {}
validSolverIDs["GMG"] = 0
validSolverIDs["GS"] = 0
validSolverIDs["ILU"] = 0
if (validSolverIDs[solverID] == nil) then
    error("Unknown solver identifier " .. solverID)
end

-- error tolerance for Limex iteration
tol = util.GetParamNumber("-tol", 0.01)

-- specify "-verbose" to output linear solver convergence
verbose = util.HasParamOption("-verbose")

-- parameters for instationary simulation
dt = util.GetParamNumber("-dt", 1e-2)
endTime = util.GetParamNumber("-endTime", 1.0)

-- choose outfile directory
outDir = util.GetParam("-outName", "wave3d")
outDir = outDir .. "/"

-- specify -vtk to generate vtk output
generateVTKoutput = util.HasParamOption("-vtk")
pstep = util.GetParamNumber("-pstep", dt, "plotting interval")

-------------------------
--  problem constants  --
-------------------------
-- setting-dependent variables
withIP3R = true
withRyR = true
withSERCAandLeak = true

if setting == "none" then 
  withIP3R = false
  withRyR = false
  withSERCAandLeak = false
end

if setting == "ip3r" then
  withRyR = false
end

if setting == "ryr" then
  withIP3R = false
end


-- total cytosolic calbindin concentration
-- (four times the real value in order to simulate four binding sites in one)
totalClb = totalBuffer

-- diffusion coefficients
D_cac = 220.0
D_cae = 220.0
D_ip3 = 280.0
D_clb = 20.0

-- calbindin binding rates
k_bind_clb = 27.0e06
k_unbind_clb = 19

-- initial concentrations
ca_cyt_init = 5.0e-08
ca_er_init = 2.5e-4
ip3_init = 4.0e-8
clb_init = totalClb / (k_bind_clb/k_unbind_clb*ca_cyt_init + 1)




-- IP3 constants
reactionRateIP3 = 0.11
equilibriumIP3 = 4.0e-08
reactionTermIP3 = -reactionRateIP3 * equilibriumIP3

-- ER densities
IP3Rdensity = 17.3 


local v_s = 6.5e-27  -- V_S param of SERCA pump
local k_s = 1.8e-7   -- K_S param of SERCA pump
local j_ip3r = 3.7606194166520605e-23   -- single channel IP3R flux (mol/s) - to be determined via gdb
local j_ryr = 1.1201015633466695e-21    -- single channel RyR flux (mol/s) - to be determined via gdb
                      -- ryr1: 1.1204582669024472e-21 
---[[
-- equilibration using SERCA
leakERconstant = 3.8e-17
local j_leak = ca_er_init-ca_cyt_init -- leak proportionality factor


SERCAfluxDensity = leakERconstant * j_leak
if withIP3R then 
  SERCAfluxDensity = SERCAfluxDensity + IP3Rdensity * j_ip3r
end
if withRyR then
  SERCAfluxDensity = SERCAfluxDensity + ryrDensity * j_ryr
end
SERCAdensity = SERCAfluxDensity / (v_s/(k_s/ca_cyt_init+1.0)/ca_er_init)
if (SERCAdensity < 0) then error("SERCA flux density is outward for these density settings!") end
--]]
--[[
-- equilibration using leakage
SERCAdensity = 1973.0
SERCAflux = v_s / (k_s / ca_cyt_init + 1.0) / ca_er_init

netEquilFlux = SERCAdensity*SERCAflux
if withIP3R then 
  netEquilFlux = netEquilFlux - IP3Rdensity * j_ip3r
end
if withRyR then
  netEquilFlux = netEquilFlux - ryrDensity * j_ryr
end

leakERconstant = netEquilFlux / (ca_er_init - ca_cyt_init)
if (leakERconstant < 0) then
  error("ER leakage flux density is outward for these density settings!")
end
--]]

-- PM densities

vdccDensity = 0.0  -- 1.0
leakPMconstant =  pmcaDensity * 6.9672131147540994e-24  -- single pump PMCA flux (mol/s)
        + ncxDensity *  6.7567567567567566e-23  -- single pump NCX flux (mol/s)
        --+ vdccDensity * (-1.5752042094823713e-25)    -- single channel VGCC flux (mol/s)
        -- *1.5 // * 0.5 for L-type // T-type
if (leakPMconstant < 0) then error("PM leak flux is outward for these density settings!") end


-- activation pattern
caEntryDuration = 0.01 -- was 0.001
--caEntryDelay = 0.5
function synCurrentDensityCa(x, y, z, t, si)  
  -- single spike (~1200 ions)
  local influx
  if t <= caEntryDuration --and t >= caEntryDelay
  then influx = 2.5e-2 * (1.0 - t/caEntryDuration) * 10 -- 10x stronger!
  else influx = 0.0
  end
    return influx
end

ip3EntryDelay = 0.000
ip3EntryDuration = 0.2
function synCurrentDensityIP3(x, y, z, t, si)
  local influx
  if t > ip3EntryDelay and t <= ip3EntryDelay+ip3EntryDuration
  then influx = 7.5e-5 * (1.0 - t/ip3EntryDuration) * 10 -- 10x stronger!
  else influx = 0.0
  end
  
    return influx
end


-------------------------------
-- setup approximation space --
-------------------------------
-- required subsets
reqSubsets = 
{
  "cyt",
  "er",
  "pm",
  "erm",
  "meas1", "meas1ER",
  "meas2", "meas2ER",
  "meas2r", "meas2rER",
  "meas3", "meas3ER",
  "meas4", "meas4ER",
  "meas5", "meas5ER",
}

-- pre refinements
numRefs = util.GetParamNumber("-numPreRefs", 0)

-- load domain
dom = util.CreateDomain(gridName, numRefs, reqSubsets, not integrityCheck)

-- element quality statistics
-- print("Quality...")
-- ElementQualityStatistics(dom:grid(), 3)
-- exit()

-- compose domain parts
cytVol = "cyt, meas1, meas2, meas2r, meas3, meas4, meas5" -- outer domain
erVol = "er, meas1ER, meas2ER, meas2rER, meas3ER, meas4ER, meas5ER" -- inner domain
plMem = "pm" -- plasma membrane
plMem_vec = {"pm"} 
erMem = "erm" -- endoplasmic reticulum membrane
erMemVec = {"erm"}
bdry = "meas2, meas2r" -- non-internal boundaries only (outer domain!)

-- compose domain parts for function spaces
outerDomain = cytVol .. ", " .. plMem .. ", " .. erMem .. ", ".. bdry
innerDomain = erVol .. ", " .. erMem

-- create approximation space
approxSpace = ApproximationSpace(dom)
approxSpace:add_fct("ca_cyt", "Lagrange", 1, outerDomain) 
approxSpace:add_fct("ca_er", "Lagrange", 1, innerDomain)
approxSpace:add_fct("clb", "Lagrange", 1, outerDomain)
approxSpace:add_fct("ip3", "Lagrange", 1, outerDomain)
approxSpace:add_fct("o2", "Lagrange", 1, "erm")
approxSpace:add_fct("c1", "Lagrange", 1, "erm")
approxSpace:add_fct("c2", "Lagrange", 1, "erm")
approxSpace:init_levels()
approxSpace:init_surfaces()
approxSpace:init_top_surface()
approxSpace:print_layout_statistic()
approxSpace:print_statistic()

-- use block algebra or not
if useBlockAlgebra then OrderCuthillMcKee(approxSpace, true) end

-- in parallel environments: domain distribution
-- balancer.partitioner = "parmetis"
balancer.partitioner = "dynBisection"
-- balancer.partitioner = ""
balancer.staticProcHierarchy = true
balancer.firstDistLvl = -1
balancer.redistSteps = 0
-- balancer.parallelElementThreshold = 4
balancer.ParseParameters()
balancer.PrintParameters()
loadBalancer = balancer.CreateLoadBalancer(dom)

if loadBalancer ~= nil then
  loadBalancer:enable_vertical_interface_creation(solverID == "GMG")
  if balancer.partitioner == "parmetis" then
    mu = ManifoldUnificator(dom)
    mu:add_protectable_subsets("erm")
    cdgm = ClusteredDualGraphManager()
    cdgm:add_unificator(SiblingUnificator())
    cdgm:add_unificator(mu)
    balancer.defaultPartitioner:set_dual_graph_manager(cdgm)
  end
  balancer.Rebalance(dom, loadBalancer)
  loadBalancer:estimate_distribution_quality()
  loadBalancer:print_quality_records()
  if balancer.partitioner == "parmetis" then
    print("Edge cut on base level: "..balancer.defaultPartitioner:edge_cut_on_lvl(0))
  end
end

strat = SurfaceMarking(dom)
strat:add_surface("erm", "cyt")
refiner = HangingNodeDomainRefiner(dom)
for i = 1, numERMRefs do
  strat:mark_without_error(refiner, approxSpace)
  refiner:refine()
end


-- more global refs (post refs after domain distribution)
if numGlobRefs > 0 then 
  local refiner = GlobalDomainRefiner(dom)  
  for i = 1, numGlobRefs do
    refiner:refine()
  end
end

-- print some statistics on parallel layout
print(dom:domain_info():to_string())
SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), outDir .. "grid/refined_grid_hierarchy_p" .. ProcRank() .. ".ugx", 1.0)
SaveParallelGridLayout(dom:grid(), outDir .. "grid/parallel_grid_layout_p"..ProcRank()..".ugx", 1.0)

--------------------------
-- setup discretization --
--------------------------

-- diffusion --
diffCaCyt = ConvectionDiffusion("ca_cyt", cytVol, "fv1")
diffCaCyt:set_diffusion(D_cac)

diffCaER = ConvectionDiffusion("ca_er", "er", "fv1")
diffCaER:set_diffusion(D_cae)

diffClb = ConvectionDiffusion("clb", cytVol, "fv1")
diffClb:set_diffusion(D_clb)


diffIP3 = ConvectionDiffusion("ip3", cytVol, "fv1")
diffIP3:set_diffusion(D_ip3)
diffIP3:set_reaction_rate(reactionRateIP3)
diffIP3:set_reaction(reactionTermIP3)



-- buffering --
discBuffer = BufferFV1(cytVol) -- where buffering occurs
discBuffer:add_reaction(
  "clb",                     -- the buffering substance
  "ca_cyt",                  -- the buffered substance
  totalClb,                  -- total amount of buffer
  k_bind_clb,    -- binding rate constant
  k_unbind_clb  -- unbinding rate constant
)


-- er membrane transport systems
ip3r = IP3R({"ca_cyt", "ca_er", "ip3"})
ip3r:set_scale_inputs({1e3,1e3,1e3})
ip3r:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ryr = RyRImplicit({"ca_cyt", "ca_er", "o2", "c1", "c2"}, {"erm"})
ryr:set_scale_inputs({1e3, 1e3, 1.0, 1.0, 1.0})
ryr:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

serca = SERCA({"ca_cyt", "ca_er"})
serca:set_scale_inputs({1e3,1e3})
serca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakER = Leak({"ca_er", "ca_cyt"})
leakER:set_scale_inputs({1e3,1e3})
leakER:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)


discIP3R = MembraneTransportFV1(erMem, ip3r)
discIP3R:set_density_function(IP3Rdensity)

discRyR = MembraneTransportFV1(erMem, ryr)
discRyR:set_density_function(ryrDensity)

discSERCA = MembraneTransportFV1(erMem, serca)
discSERCA:set_density_function(SERCAdensity)


--discERLeak = MembraneTransportFv1(erMem, leakER)
discERLeak = MembraneTransportFV1(erMem, leakER)
discERLeak:set_density_function(1e12*leakERconstant/(1e3)) -- from mol/(um^2 s M) to m/s


-- plasma membrane transport systems
pmca = PMCA({"ca_cyt", ""})
pmca:set_constant(1, 1.0)
pmca:set_scale_inputs({1e3,1.0})
pmca:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

ncx = NCX({"ca_cyt", ""})
ncx:set_constant(1, 1.0)
ncx:set_scale_inputs({1e3,1.0})
ncx:set_scale_fluxes({1e15}) -- from mol/(um^2 s) to (mol um)/(dm^3 s)

leakPM = Leak({"", "ca_cyt"})
leakPM:set_constant(0, 1.0)
leakPM:set_scale_inputs({1.0,1e3})
leakPM:set_scale_fluxes({1e3}) -- from mol/(m^2 s) to (mol um)/(dm^3 s)


discPMCA = MembraneTransportFV1("pm", pmca)
discPMCA:set_density_function(pmcaDensity)

discNCX = MembraneTransportFV1("pm", ncx)
discNCX:set_density_function(ncxDensity)

discPMLeak = MembraneTransportFV1("pm", leakPM)
discPMLeak:set_density_function(1e12*leakPMconstant / (1.0-1e3*ca_cyt_init))



-- synaptic activity
synapseInfluxCa = UserFluxBoundaryFV1("ca_cyt", "meas2")
synapseInfluxCa:set_flux_function("synCurrentDensityCa")
synapseInfluxIP3 = UserFluxBoundaryFV1("ip3", "meas2")
synapseInfluxIP3:set_flux_function("synCurrentDensityIP3")


-- domain discretization --
domDisc = DomainDiscretization(approxSpace)

domDisc:add(diffCaCyt)
domDisc:add(diffCaER)
domDisc:add(diffClb)
domDisc:add(diffIP3)
domDisc:add(discBuffer)


if withIP3R then
  domDisc:add(discIP3R)
end
if withRyR then
 domDisc:add(discRyR)
 domDisc:add(ryr) -- also add ryr as elem disc (for state variables)
end
if withSERCAandLeak then
  domDisc:add(discSERCA)
  domDisc:add(discERLeak)
end

domDisc:add(discPMCA)
domDisc:add(discNCX)
domDisc:add(discPMLeak)

domDisc:add(synapseInfluxCa)
if withIP3R then
  domDisc:add(synapseInfluxIP3)
end


-- setup time discretization --
timeDisc = ThetaTimeStep(domDisc)
timeDisc:set_theta(1.0) -- 1.0 is implicit Euler

-- create operator from discretization
op = AssembledOperator()
op:set_discretization(timeDisc)
op:init()



------------------
-- solver setup --
------------------
-- debug writer
dbgWriter = GridFunctionDebugWriter(approxSpace)
dbgWriter:set_base_dir(outDir)
dbgWriter:set_vtk_output(false)

-- biCGstab --
convCheck = ConvCheck()
convCheck:set_minimum_defect(1e-50)
convCheck:set_reduction(1e-8)
convCheck:set_verbose(verbose)


if (solverID == "ILU") then
    bcgs_steps = 1000
    ilu = ILU()
    ilu:set_sort(true)
    bcgs_precond = ilu
elseif (solverID == "GS") then
    bcgs_steps = 1000
    bcgs_precond = GaussSeidel()
  --  bcgs_precond:set_debug(GridFunctionDebugWriter(approxSpace))
else -- (solverID == "GMG")
  gmg = GeometricMultiGrid(approxSpace)
  gmg:set_discretization(timeDisc)
  gmg:set_base_level(0)
  gmg:set_gathered_base_solver_if_ambiguous(true)
  
  -- treat SuperLU problems with Dirichlet constraints by using constrained version
  gmg:set_base_solver(SuperLU())
 --  gmg:set_base_solver(LU())
  
  smoother = GaussSeidel()
  gmg:set_smoother(smoother)
  gmg:set_smooth_on_surface_rim(true)
  gmg:set_cycle_type(1)
  gmg:set_num_presmooth(3)
  gmg:set_num_postsmooth(3)
  --gmg:set_rap(true)
 -- gmg:set_debug(GridFunctionDebugWriter(approxSpace))
  
    bcgs_steps = 1000
  bcgs_precond = gmg
end

convCheck:set_maximum_steps(bcgs_steps)

bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(bcgs_precond)
bicgstabSolver:set_convergence_check(convCheck)
--bicgstabSolver:set_debug(dbgWriter)

--- non-linear solver ---
-- convergence check
newtonConvCheck = CompositeConvCheck(approxSpace, 10, 1e-17, 1e-08)
--newtonConvCheck:set_component_check("ca_cyt, ca_er, clb, ip3", 1e-18, 1e-10)
newtonConvCheck:set_verbose(true)
newtonConvCheck:set_time_measurement(true)
--newtonConvCheck:set_adaptive(true)


-- Newton solver
newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(bicgstabSolver)
--newtonSolver:set_convergence_check(newtonConvCheck)
--newtonSolver:set_debug(dbgWriter)
--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))

newtonSolver:init(op)

if (generateVTKoutput) then out = VTKOutput() end

-- get grid function
u = GridFunction(approxSpace)

-- set initial value
InterpolateInner(ca_cyt_init, u, "ca_cyt", 0.0)
InterpolateInner(ca_er_init, u, "ca_er", 0.0)
InterpolateInner(clb_init, u, "clb", 0.0)
InterpolateInner(ip3_init, u, "ip3", 0.0)

if withRyR then
  ryr:calculate_steady_state(u)
end

-- begin time loop
time = 0.0
uOld = u:clone()
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

numSteps = util.GetParamNumber("-numTimeSteps", 10000)
writeInterval = util.GetParamNumber("-writeInterval", 100)
for step=0, numSteps do
  print("Time step: " .. step)
  timeDisc:prepare_step(solTimeSeries, dt)
  newtonSolver:apply(u)

  time = solTimeSeries:time(0) + dt
 
  name = vtkFolder .. "/solution3d" 
  if (generateVTKoutput and (step % writeInterval == 0)) then out:print(name, u, step, time) end

  oldestSol = solTimeSeries:oldest()
  VecScaleAssign(oldestSol, 1.0, u)
  solTimeSeries:push_discard_oldest(oldestSol, time)
  outPath = vtkFolder .. "/meas/meas"

  -- at source
  take_measurement(u, time, "meas2", "ca_cyt", outPath)
  take_measurement(u, time, "meas2ER", "ca_er", outPath)
  -- 30 micron away from source downstream (after first BP)
  take_measurement(u, time, "meas1", "ca_cyt", outPath)
  take_measurement(u, time, "meas1ER", "ca_er", outPath)
  -- 50 micron away from source downstream (before first BP)
  take_measurement(u, time, "meas3", "ca_cyt", outPath)
  take_measurement(u, time, "meas3ER", "ca_er", outPath)
  -- 70 micron away from source downstream (before first BP)
  take_measurement(u, time, "meas4", "ca_cyt", outPath)
  take_measurement(u, time, "meas4ER", "ca_er", outPath)
  -- 80 micron away from source downstream (before first BP)
  take_measurement(u, time, "meas5", "ca_cyt", outPath)
  take_measurement(u, time, "meas5ER", "ca_er", outPath)
  -- at soma
  take_measurement(u, time, "meas2r", "ca_cyt", outPath)
  take_measurement(u, time, "meas2rER", "ca_er", outPath)

end

if (generateVTKoutput) then out:write_time_pvd(name, u) end
