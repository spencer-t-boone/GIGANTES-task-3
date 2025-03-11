
# Generate a LPO with multiple shooting
from godot.core import num, tempo, astro, events, ipfwrap
from godot.core import autodif as ad
from godot.core.autodif import bridge as br
import godot.core.util as c_util
from godot.model import interface, frames, common, prop
from godot import cosmos
from godot.cosmos import util
import numpy as np
import matplotlib.pyplot as plt # Plotting like maltab
from ruamel import yaml
import time, os, copy
import pygmo as pg
import pygmo_plugins_nonfree as ppnf
import midas

from aux_traj_fun import add_ctr, add_man, add_match, create_trajectory_config, add_end_point

# Avoid propagator log
c_util.suppressLogger()

# Setup WORHP optimiser
# worhp_path = os.environ.get('WORHP_BINARY')
# worhp = ppnf.worhp()
# # Setup optimiser options
# worhp.set_integer_option('MaxIter', 100)
# worhp.set_integer_option('MaxCalls', 5000)
# worhp.set_numeric_option('TolOpti', 1e-3)
# worhp.set_integer_option('Algorithm', 1)
# worhp.set_bool_option('FeasibleInit', False)


# Create universe
config_uni = util.load_yaml('universe_1.yml')
universe = cosmos.Universe( config_uni )
XMUE=universe.constants.getMu('Saturn') 
XMUM=universe.constants.getMu('Enceladus') 
R_earth=universe.constants.getRadius('Saturn') 
R_moon=universe.constants.getRadius('Enceladus') 
LU_cr3bp = 238000
MU = XMUM/(XMUM + XMUE)
TU_cr3bp = np.sqrt(LU_cr3bp**3/(XMUM+XMUE))


# Select specific orbit to analyse 
orbit_types = ['L1 nrho', 'L2 nrho', 'butterfly', 'period-3']
orbit_type = orbit_types[3]
MU = 1.9011497893288988e-7

if orbit_type == 'L1 nrho':
    x0_cr3bp = np.array([0.9999391020169318, 0, -0.001183305369127517, 0, -0.01686127728398247, 0])
    period = 2.282108846114326*TU_cr3bp/86400
    n_orb = 14
    n_points = 4
    scales_vect = [25e0, 25e0, 25e0, 0.2, 0.2, 0.2]
    delta_vect = np.array([25e0, 25e0, 25e0, 0.2, 0.2, 0.2])

elif orbit_type == 'L2 nrho':
    x0_cr3bp = np.array([1.000064575406844, 0, -0.001183305369127517, 9.288141002844474e-13, 0.01685176088135196, 1.251253788081552e-11])
    period = 2.287781993612279*TU_cr3bp/86400
    n_orb = 15
    n_points = 4
    scales_vect = [30e0, 30e0, 30e0, 0.2, 0.2, 0.2]
    delta_vect = np.array([30e0, 30e0, 30e0, 0.2, 0.2, 0.2])

elif orbit_type == 'butterfly':
    x0_cr3bp = np.array([0.9979394855172609, 0, 0.003723544484837956, 0, -0.0008522149920766561, 0]) # Apoapse
    period = 3.593492349524928*TU_cr3bp/86400
    n_orb = 20
    n_points = 4
    scales_vect = [50e0, 50e0, 50e0, 0.2, 0.2, 0.2]
    delta_vect = np.array([50e0, 50e0, 50e0, 0.2, 0.2, 0.2])
    

elif orbit_type == 'period-3':
    x0_cr3bp = np.array([0.9973675664626385, 0, 0.004768373628104226, 0, 0.00570366409390525, 0]) # Apoapse
    period = 6.856114820633476*TU_cr3bp/86400
    n_orb = 6
    n_points = 12
    scales_vect = [25e0, 25e0, 25e0, 0.2, 0.2, 0.2]
    delta_vect = np.array([25e0, 25e0, 25e0, 0.2, 0.2, 0.2])




# CR3BP orbit
CRTBP = midas.astro.crtbp.AdimensionalCRTBP(MU)
IntegratorCRTBP = midas.astro.crtbp.AdimensionalIntegratorCRTBP(MU)

#################################################################################
# Create multiple shooting points
E0 = tempo.Epoch('2040-05-01T00:00:00.0 TDB')
# get points
dt = period/n_points
mjd_vect_ref = E0.mjd() + np.arange(0, n_orb*period + dt, dt) 
x0_vect_ref = []
for m_ in mjd_vect_ref:
    # Get CR3BP state
    dt_cr3bp = (m_ - E0.mjd()) % period * 86400/TU_cr3bp
    [t_cr3bp, x_cr3bp] = IntegratorCRTBP.integrate(x0_cr3bp, dt_cr3bp)
    x0_rot = np.copy(x_cr3bp[-1])
    # Moon-centred dimensional state
    lu = universe.frames.distance('Saturn','Enceladus', tempo.Epoch(str(m_) + ' TDB'))
    tu = np.sqrt(lu**3/(XMUM+XMUE))
    x0_rot[0] -= 1 - MU
    x0_rot[:3] *= lu
    x0_rot[3:] *= lu/tu
    x0_vect_ref.append(x0_rot)
    

#################################################################################
# Prepare timeline
config_traj = create_trajectory_config()
TT = []
i_ = 0
x_point_ref = []
# Add points
for m_, x_ in zip(mjd_vect_ref, x0_vect_ref):
    add_ctr(TT, i_, m_ + dt, x_)
    if i_ < n_orb*n_points:
        add_man(TT, i_)
        add_match(TT, i_, dt)
    i_ += 1
    x_point_ref.append(x_)
# add final point
add_end_point(TT, i_-1, 1.0)

# Assemble trajectory
config_traj['timeline'] = TT
traj = cosmos.Trajectory(universe, config_traj)

#################################################################################
# Assemble problem
config_prob = {
    'objective': {'type': 'minimise',
                  'point': 'end_point',
                  'value': 'SC_dv',
                  'scale': 1e3},
    'parameters': {'free': []}}
FF = []
SS = {}
BB = {}
cart_vect = ['pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z']

units_vect = [' km', ' km', ' km', ' km/s', ' km/s', ' km/s']

for i_ in range(n_points*n_orb):
    if not (i_==0): 
        FF.append('ctr' + str(i_)+'_dt')
        FF.append('ctr' + str(i_)+'_SC_mass')
        FF.append('ctr' + str(i_)+'_SC_dv')
        SS['ctr' + str(i_)+'_dt'] = '1 day'
        SS['ctr' + str(i_)+'_SC_mass'] = '100 kg'
        SS['ctr' + str(i_)+'_SC_dv'] = '1 mm/s'
        BB['ctr' + str(i_)+'_dt'] = ['-0.5 day', '0.5 day']
        BB['ctr' + str(i_)+'_SC_dv'] = ['0 m/s', '1000 m/s']
        BB['ctr' + str(i_)+'_SC_mass'] = ['0 kg', '1000 kg']
    for c_, s_, u_, b_, d_ in zip(cart_vect, scales_vect, units_vect, x_point_ref[i_],delta_vect):
        FF.append('ctr'+str(i_)+'_SC_center_'+c_)
        SS['ctr'+str(i_)+'_SC_center_'+c_] = str(s_)+u_
        BB['ctr'+str(i_)+'_SC_center_'+c_] = [str(b_-d_)+u_,str(b_+d_)+u_]

    # Add manoeuvre
    FF.append('man'+str(i_)+'_dv1')
    FF.append('man'+str(i_)+'_dv2')
    FF.append('man'+str(i_)+'_ras')
    FF.append('man'+str(i_)+'_dec')
    SS['man'+str(i_)+'_dv1'] = '1 mm/s'
    SS['man'+str(i_)+'_dv2'] = '1 mm/s'
    BB['man'+str(i_)+'_dv1'] = ['0 m/s', '2 m/s']
    BB['man'+str(i_)+'_dv2'] = ['0 m/s', '2 m/s']
    # Add match point
    FF.append('match'+str(i_)+'_right_dt')
    SS['match'+str(i_)+'_right_dt'] = '1 day'
    BB['match'+str(i_)+'_right_dt'] = ['1 s', str(period/n_points+1)+' day']
# Add last point
i_ = n_points*n_orb 
FF.append('ctr' + str(i_)+'_dt')
FF.append('ctr' + str(i_)+'_SC_mass')
FF.append('ctr' + str(i_)+'_SC_dv')
SS['ctr' + str(i_)+'_dt'] = '1 day'
SS['ctr' + str(i_)+'_SC_mass'] = '100 kg'
SS['ctr' + str(i_)+'_SC_dv'] = '1 mm/s'
BB['ctr' + str(i_)+'_dt'] = ['-0.5 day', '0.5 day']
BB['ctr' + str(i_)+'_SC_dv'] = ['0 m/s', '100 m/s']
BB['ctr' + str(i_)+'_SC_mass'] = ['0 kg', '1000 kg']
for c_, s_, u_, b_, d_ in zip(cart_vect, scales_vect, units_vect, x_point_ref[0],delta_vect):
    FF.append('ctr'+str(i_)+'_SC_center_'+c_)
    SS['ctr'+str(i_)+'_SC_center_'+c_] = str(s_)+u_
    BB['ctr'+str(i_)+'_SC_center_'+c_] = [str(b_-d_)+u_,str(b_+d_)+u_]

FF.pop(0) # remove initial x coordinate
FF.pop(0) # remove initial y coordinate
FF.pop(-5) # remove final y coordinate
# FF.pop(-5) # remove final x coordinate


# Setup bounds on left and right points
config_prob['parameters']['free'] = FF
config_prob['parameters']['scales'] = SS
config_prob['parameters']['bounds'] = BB
util.save_yaml(config_prob, "prob_out.yml")
prob = cosmos.Problem(universe, [traj], config_prob, useGradient=True)


from aux_bisection_fun import DummyConstraint
dummy_con = DummyConstraint()
prob.addConstraint('dummy', dummy_con)

# plot initial guess
fig,ax = plt.subplots(1,3)
ax[0].grid()
ax[0].set_aspect('equal')
ax[0].set_xlabel('$X_{rot}$ (km)')
ax[0].set_ylabel('$Y_{rot}$ (km)')
ax[1].grid()
ax[1].set_aspect('equal')
ax[1].set_xlabel('$X_{rot}$ (km)')
ax[1].set_ylabel('$z_{rot}$ (km)')
ax[2].grid()
ax[2].set_aspect('equal')
ax[2].set_xlabel('$Y_{rot}$ (km)')
ax[2].set_ylabel('$z_{rot}$ (km)')

traj.compute(partials=False)
XROT=[]
# E_grid = tempo.EpochRange(traj.point('ctr0'),traj.point('end_point')).createGrid(3600)
# for E_ in E_grid:
#     x_rot = universe.frames.vector6('Enceladus','SC_center', 'SEROT', E_)
#     XROT.append(np.concatenate(([E_.mjd()],x_rot)))
# ax[0].plot([x_[1] for x_ in XROT],[x_[2] for x_ in XROT],'k--')
# ax[1].plot([x_[1] for x_ in XROT],[x_[3] for x_ in XROT],'k--')


# Prepare Pygmo problem
problem = pg.problem(prob)
tol_con = 1e-5
problem.c_tol = [tol_con] * problem.get_nc()


# Prepare population with initial guess
x0 = prob.get_x()
pop = pg.population(problem, 0)
pop.push_back(x0)

# Solve!
#algo = pg.algorithm(worhp)
ip = pg.ipopt()
ip.set_numeric_option("tol",1e-3)
algo = pg.algorithm(ip)
algo.set_verbosity(1)
pop = algo.evolve(pop)

#################################################################################
# Solution
traj.compute(partials=False)
traj_up = traj.applyParameterChanges()
conf_update = util.deep_update(config_traj, traj_up)
#util.save_yaml(conf_update, "trajectories/traj_out_multiple_shooting.yml")

dv_tot = universe.evaluables.get('SC_dv').eval(traj.point('end_point'))
print(f'Total dv: {dv_tot*1e3:.2e} m/s')
t_tot = traj.point('end_point') - traj.point('ctr0')
print(f'Total time: {t_tot/86400:.2f} days')

E_grid = tempo.EpochRange(traj.point('ctr0'),traj.point('end_point')).createGrid(60)
XROT=[]
for E_ in E_grid:
    x_rot = universe.frames.vector6('Enceladus','SC_center', 'SEROT', E_)
    XROT.append(np.concatenate(([E_.mjd()],x_rot)))
ax[0].plot([x_[1] for x_ in XROT],[x_[2] for x_ in XROT])
ax[1].plot([x_[1] for x_ in XROT],[x_[3] for x_ in XROT])
ax[2].plot([x_[2] for x_ in XROT],[x_[3] for x_ in XROT])
#fig.savefig('figures/multiple_shooting_fig.png')
plt.show(block = False)

