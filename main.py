# Firstly, specify location of netcdf input files
node = 'local'
if node == 'local':
    import sys
    sys.path.insert(0, 'C:/Users/ciank/PycharmProjects/sinmod/opendrift')  # add opendrift local path
    infile_path = 'A:/Cian_sinmod/copernicus_client/results/'
    outfile_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/results/'
else:  # assume job is on server
    infile_path = '/cluster/projects/nn9828k/Cian_sinmod/copernicus_client/results/'
    outfile_path = infile_path

from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic, reader_global_landmask
from configure import Scenario
from post_process import Process

# Simulation settings (time, releases, initialization scenario)
file_prefix = 'CMEMS_GLPHYS_D_full_'  # File identifier
y_start = 2000  # first year of simulation
y_end = 2004  # final year of simulation (note: only used if release in [y_end - 1] extends into [y_end])
time_step_hours = 1/4  # simulation time step (negative time is backwards stepping of model)
save_time_step_hours = 1  # save time step
duration_days = 200  # simulation duration in days;
release_end = 20   # total number of releases for simulation
release_n_days = 5  # number of days between releases (time=start_time + i*time_step)
release_step = 24*release_n_days  # number of hours between releases
init_keys = ["APSO"]  # key names for initialization scenario: defines lat-long start points, number of particles etc.

# test processing
test = False

for y_i in range(y_start, y_end):
    for r_i in range(0, release_end):
        for key in init_keys:

            # Input netcdf file with physics (u, v ,T ...):
            phys_states = infile_path + file_prefix + str(y_i) + '.nc'
            print('Beginning simulation: year = ' + str(y_i) + ', release number ' + str(r_i))

            # Initialize OpenDrift object (model type)
            o = OceanDrift(loglevel=20)  # log_level= 0 for full diagnostics, 50 for none
            reader_samples = reader_netCDF_CF_generic.Reader(phys_states)  # read input variables
            reader_landmask = reader_global_landmask.Reader()  # high resolution coast for particle beaching etc.
            o.add_reader([reader_landmask, reader_samples])  # add readers to model instance

            # Parameterization of scenario object for initialization and running
            scenario = Scenario(infile_path, outfile_path, reader_samples, duration_days,
                                time_step_hours, save_time_step_hours, release_step, y_i, r_i, key)

            # Initialization of simulation
            o.disable_vertical_motion()
            o.seed_elements(lon=scenario.lon_init,
                            lat=scenario.lat_init,
                            time=scenario.t_init,
                            number=scenario.n_part,
                            radius=scenario.radius)

            # Running of simulation
            o.run(duration=scenario.duration,
                  time_step=scenario.time_step,
                  time_step_output=scenario.save_time_step,
                  outfile=scenario.trajectory_file,
                  export_variables=scenario.export_variables)

            # Post-process simulation file, saving intermediate data (unique particle visits, transit times ...)
            Process(scenario.trajectory_file, outfile_path, y_i, r_i, key, test)










