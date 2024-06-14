# Firstly, specify location of netcdf input files
from configure import Scenario, FileExplorer
fpath = FileExplorer(node='local', model_name='cmems', key="APSO")
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic, reader_global_landmask

# Simulation settings (time, releases, initialization scenario)
y_start = 2017  # first year of simulation #todo: define start and end year instead of loop as below
m_start = 1  # start month
d_start = 1  # start day of month - otherwise
m_end = 12  # end month
time_step_hours = 1  # simulation time step (negative time is backwards stepping of model)
save_time_step_hours = 6  # save time step
duration_days = 365  # simulation duration in days;
release_end = 1   # total number of releases for simulation
release_n_days = 1  # number of days between releases (time=start_time + i*time_step)
release_step = 24*release_n_days  # number of hours between releases
phys_states = fpath.check_phys_states(y_start, m_start, m_end)  # input netcdf file with physics (u, v ,T ...)


for r_i in range(0, release_end):
    print('Beginning simulation: year = ' + str(y_start) + ', release number ' + str(r_i + 1))
    reader_phys_states = reader_netCDF_CF_generic.Reader(phys_states)  # read input variables
    print(reader_phys_states)


    # initialize OpenDrift object (model type);
    o = OceanDrift(loglevel=20)  # log_level= 0 for full diagnostics, 50 for none
    reader_landmask = reader_global_landmask.Reader()  # high resolution coast for particle beaching etc.
    o.add_reader([reader_landmask, reader_phys_states])

    # Parameterization of scenario object for initialization and running
    scenario = Scenario(fpath, reader_phys_states, duration_days,
                        time_step_hours, save_time_step_hours, release_step,
                        y_start, r_i, d_start, m_start, m_end)

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














