# Firstly, specify location of netcdf input files
from configure import Scenario, FileExplorer
fpath = FileExplorer(node='local', model_name='sinmod')
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic, reader_global_landmask

# Simulation settings (time, releases, initialization scenario)
y_start = 2020  # first year of simulation
y_end = 2021  # final year of simulation (note: only used if release in [y_end - 1] extends into [y_end])
time_step_hours = 6  # simulation time step (negative time is backwards stepping of model)
save_time_step_hours = 6  # save time step
duration_days = 28  # simulation duration in days;
release_end = 1   # total number of releases for simulation
release_n_days = 1  # number of days between releases (time=start_time + i*time_step)
release_step = 24*release_n_days  # number of hours between releases
init_keys = ["SG8H"]  # key names for initialization scenario: defines lat-long start points, number of particles etc.

for y_i in range(y_start, y_end):
    for r_i in range(0, release_end):
        for key in init_keys:
            # Input netcdf file with physics (u, v ,T ...):
            if fpath.model_name == "sinmod":
                phys_states = fpath.phys_states_path + fpath.phys_states_file_prefix + str(y_i) + '01' + '.nc' #todo: generalise to multiple months
            else:
                phys_states = fpath.phys_states_path + fpath.phys_states_file_prefix + str(y_i) + '.nc'
            print('Beginning simulation: year = ' + str(y_i) + ', release number ' + str(r_i+1))

            # Initialize OpenDrift object (model type)
            o = OceanDrift(loglevel=20)  # log_level= 0 for full diagnostics, 50 for none
            reader_phys_states = reader_netCDF_CF_generic.Reader(phys_states)  # read input variables
            reader_landmask = reader_global_landmask.Reader()  # high resolution coast for particle beaching etc.
            o.add_reader([reader_landmask, reader_phys_states])  # add readers to model instance

            # Parameterization of scenario object for initialization and running
            scenario = Scenario(fpath, reader_phys_states, duration_days,
                                time_step_hours, save_time_step_hours, release_step, y_i, r_i, key)

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














