# Firstly, specify location of netcdf input files
from configure import Scenario, FileExplorer
fpath = FileExplorer(node='local', model_name='cmems', key="APSO")
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic, reader_global_landmask
import datetime

# Simulation settings (time, releases, initialization scenario)
date_init = datetime.datetime(2017, 1, 1, 0, 0)  # beginning of first simulation;
duration_days = datetime.timedelta(days=60)  # simulation duration in days;
release_step = datetime.timedelta(days=2)  # days between releases of particles;
release_n = 3   # total number of releases for simulation
date_limit = date_init + (release_step*(release_n-1)) + duration_days  # final date that will be accessed
time_step_hours = datetime.timedelta(hours=6)  # simulation time step (negative time is backwards stepping of model)
save_time_step_hours = datetime.timedelta(hours=6)  # save time step
phys_states = fpath.get_phys_states(date_init, date_limit)  # input netcdf file with physics (u, v ,T ...)

# load readers before the main loop of the simulation:
reader_phys_states = reader_netCDF_CF_generic.Reader(phys_states)  # read input variables
reader_landmask = reader_global_landmask.Reader()  # high resolution coast for particle beaching etc.

# loop over releases:
for r_i in range(0, release_n):
    date_release = date_init + (release_step*r_i)
    breakpoint()
    print('Beginning simulation: year = ' + str() + ', release number ' + str(r_i + 1))

    print(reader_phys_states)


    # initialize OpenDrift object (model type);
    o = OceanDrift(loglevel=20)  # log_level= 0 for full diagnostics, 50 for none
    o.add_reader([reader_landmask, reader_phys_states])

    # Parameterization of scenario object for initialization and running
    scenario = Scenario(fpath, reader_phys_states, duration_days,
                        time_step_hours, save_time_step_hours, release_step,
                        y_start, r_i, d_start, m_start, m_end)
    #todo: need to change arguments to scenario;

    # Initialization of simulation
    o.disable_vertical_motion()
    o.seed_elements(lon=scenario.site_lon_init,
                    lat=scenario.site_lat_init,
                    time=scenario.t_init,
                    number=scenario.n_part,
                    radius=scenario.radius)

    # Running of simulation
    o.run(duration=scenario.duration,
          time_step=scenario.time_step,
          time_step_output=scenario.save_time_step,
          outfile=scenario.trajectory_file,
          export_variables=scenario.export_variables)














