from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_netCDF_CF_generic, reader_global_landmask
from configure import Case

# Firstly, specify location of netcdf input files
node = 'local'
if node == 'local':
    infile_path = 'A:/Cian_sinmod/copernicus_client/results/'
    outfile_path = 'C:/Users/ciank/PycharmProjects/sinmod/Krill_data/SINdrift/results/'
else:  # assume job is on server
    infile_path = '/cluster/projects/nn9828k/Cian_sinmod/copernicus_client/results/'
    outfile_path = infile_path

# Simulation settings (time, releases, initialization scenario)
file_prefix = 'CMEMS_GLPHYS_D_full_'  # File identifier
y_start = 2000  # first year of simulation
y_end = 2011  # final year of simulation (note: only used if release in [y_end - 1] extends into [y_end])
time_step_hours = 6  # negative time is backwards stepping of model
duration_days = 200  # look into (time=start_time + i*time_step) for setting the start and end time of simulations;
release_end = 10   # total number of releases for simulation
release_n_days = 5  # number of days between releases
release_step = 24*release_n_days  # number of hours between releases
init_keys = ["APSO"]  # key names for initialization scenario: defines lat-long start points, number of particles etc.

for y_i in range(y_start, y_end):
    for r_i in range(0, release_end):
        for key in init_keys:

            # Specify netcdf name and print information about simulation:
            phys_states = infile_path + file_prefix + str(y_i) + '.nc'
            print('Beginning simulation: year = ' + str(y_i) + ', release number ' + str(r_i))

            # Initialize OpenDrift object
            o = OceanDrift(loglevel=20)  # log_level= 0 for full diagnostics, 50 for none
            reader_samples = reader_netCDF_CF_generic.Reader(phys_states)  # read input variables
            reader_landmask = reader_global_landmask.Reader()  # high resolution coast for particle beaching etc.
            o.add_reader([reader_landmask, reader_samples])  # add readers to model instance

            # builds parameterization for the specific case with case object
            case = Case(infile_path, outfile_path, reader_samples, duration_days, time_step_hours, release_step, y_i,
                        r_i, key)

            # Initialization of simulation
            o.disable_vertical_motion()
            o.seed_elements(lon=case.lon_init,
                            lat=case.lat_init,
                            time=case.t_init,
                            number=case.n_part,
                            radius=case.radius)

            # Running of simulation
            o.run(duration=case.duration,
                  time_step=case.time_step,
                  outfile=case.trajectory_file,
                  export_variables=case.export_variables)






#print(reader_samples)






