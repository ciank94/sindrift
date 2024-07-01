# Firstly, specify location of netcdf input files
from configure import Scenario, FileExplorer
fpath = FileExplorer(node='local', model_name='cmems', key="SOIN")
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic, reader_global_landmask
import datetime

# Simulation settings (time, releases, initialization scenario)
date_init = datetime.datetime(2016, 11, 1, 0, 0)  # beginning of first simulation;
duration_days = datetime.timedelta(days=3)  # simulation duration in days;
release_step = datetime.timedelta(days=3)  # days between releases of particles;
time_step = datetime.timedelta(hours=1)  # simulation time step (negative time is backwards stepping of model)
save_step = datetime.timedelta(hours=6)  # save time step
release_n = 1   # total number of releases for simulation
date_limit = date_init + (release_step*(release_n-1)) + duration_days  # final date that will be accessed
phys_states = fpath.get_phys_states(date_init, date_limit)  # input netcdf file with physics (u, v ,T ...)

# load readers before the main loop of the simulation:
reader_phys_states = reader_netCDF_CF_generic.Reader(phys_states)  # read input variables
print(reader_phys_states)
reader_landmask = reader_global_landmask.Reader()  # high resolution coast for particle beaching etc.
print(reader_landmask)

# loop over releases:
for release_i in range(0, release_n):
    date_release = date_init + (release_step*release_i)
    print('===========================')
    print('Beginning simulation for year: ' + str(date_release) + ', release number ' + str(release_i + 1))
    print('Simulation duration: ' + str(duration_days.days) + ' days')
    print('Simulation end: ' + str(date_release + duration_days))
    print('===========================')

    # initialize OpenDrift object (model type);
    o = OceanDrift(loglevel=20)  # log_level= 0 for full diagnostics, 50 for none
    o.add_reader([reader_landmask, reader_phys_states])  # add readers for input

    # Parameterization of scenario object for initialization and running
    scenario = Scenario(fpath, date_release, duration_days, time_step, save_step, release_step, release_i
                        , reader_phys_states)

    # Initialization of simulation
    o.disable_vertical_motion()
    o.set_config('drift:advection_scheme', 'runge-kutta')
    o.seed_elements(lon=scenario.site_lon_init,
                    lat=scenario.site_lat_init,
                    time=scenario.date_init,
                    number=scenario.n_part,
                    radius=scenario.radius)

    # Running of simulation
    o.run(duration=scenario.duration,
          time_step=scenario.time_step,
          time_step_output=scenario.save_time_step,
          outfile=scenario.trajectory_file,
          export_variables=scenario.export_variables)

    # clear reference to class after the model run:
    del o

    print('===========================')
    print('Finished simulation for year: ' + str(date_release) + ', release number ' + str(release_i + 1))
    print('Simulation duration: ' + str(duration_days.days) + ' days')
    print('Simulation end: ' + str(date_release + duration_days))
    print('===========================')
















