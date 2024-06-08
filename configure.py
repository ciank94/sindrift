from datetime import timedelta
import numpy as np
import netCDF4 as nc


class Scenario:
    def __init__(self, infile_path, outfile_path, reader_phys_states, duration_days, time_step_hours,
                 save_time_step_hours, release_step, y_i, r_i, key):
        # Simulation settings for scenario
        self.key = key  # key used to identify initialization scenario
        self.year = y_i  # simulation year
        self.release_n = r_i + 1  # release number starts at one
        self.z = 50  # release depth
        self.duration = timedelta(hours=24 * duration_days)  # simulation duration in hours as datetime object
        self.time_step = timedelta(hours=time_step_hours)  # simulation time step in hours as datetime object
        self.save_time_step = timedelta(hours=save_time_step_hours)  # how often the file is saved
        self.release_step = release_step  # number of hours between releases
        self.export_variables = ['lon', 'lat']  # choose variables to export from simulation to nc file
        self.infile_path = infile_path  # path to server with nc files
        self.outfile_path = outfile_path  # path to save the output file

        # Information from the phys_states_file:
        self.phys_states_timestep = reader_phys_states.time_step.seconds
        self.reader_t_init = reader_phys_states.start_time
        self.t_init = self.reader_t_init + r_i * timedelta(hours=self.release_step)
        self.reader_end_time = reader_phys_states.end_time

        # Initialize scenario parameters and save them as attributes to scenario file
        self.scenario_initialization()  # furnish initialization scenario with class attributes at beginning
        self.scenario_file = (self.outfile_path + self.key + '_' + str(self.year) + '_R'
                         + str(self.release_n) + '_trajectory_scenario.nc')
        self.init_scenario_netcdf()
        return

    def init_scenario_netcdf(self):
        self.outfile = nc.Dataset(self.scenario_file, 'w')

        # main simulation settings as attributes:
        self.outfile.title = 'OpenDrift trajectory analysis'
        self.outfile.infile_path = self.infile_path  # path to server with nc files
        self.outfile.outfile_path = self.outfile_path  # path to save the output file
        self.outfile.trajectory_file = self.trajectory_file
        self.outfile.scenario_key = self.key  # key name given to scenario
        self.outfile.scenario_description = self.description  # description of scenario

        # information about release:
        self.outfile.n_part = self.n_part
        self.outfile.radius = self.radius
        self.outfile.release_number = self.release_n  # release number starts at one
        self.outfile.release_depth = self.z  # release depth
        self.outfile.release_step_hours = self.release_step  # number of hours between releases


        # simulation time information
        self.outfile.sim_start_day = self.t_init.day # simulation start time
        self.outfile.sim_start_month = self.t_init.month  # simulation start time
        self.outfile.sim_start_month = self.t_init.month  # simulation start time
        self.outfile.sim_start_year = self.year  # simulation year
        self.outfile.sim_duration_days = self.duration.days # simulation duration
        self.outfile.sim_time_step_seconds = self.time_step.seconds # simulation time step
        self.outfile.sim_save_time_step_seconds = self.save_time_step.seconds # simulation save time step

        # phys_states reader information
        self.outfile.reader_time_step_seconds = self.phys_states_timestep  # phys_states time step
        #self.outfile.reader_start_time = self.reader_t_init  # reader start time
        #self.outfile.reader_end_time = self.reader_end_time  # reader end time
        self.outfile.export_variables = self.export_variables  # variables exported from simulation

        # domain and site limits:
        self.outfile.site_lon_min = self.site_lon_min
        self.outfile.site_lon_max = self.site_lon_max
        self.outfile.site_lat_min = self.site_lat_min
        self.outfile.site_lat_max = self.site_lat_max
        self.outfile.domain_lon_min = self.domain_lon_min
        self.outfile.domain_lon_max = self.domain_lon_max
        self.outfile.domain_lat_min = self.domain_lat_min
        self.outfile.domain_lat_max = self.domain_lat_max
        self.outfile.bin_resolution = self.bin_res


        # set_dimensions of file:
        self.outfile.createDimension('n_lon_bins', self.shp_lon_bins)
        self.outfile.createDimension('n_lat_bins', self.shp_lat_bins)
        self.outfile.createVariable('lon_bin_vals', 'f4', ('n_lon_bins', ))
        self.outfile.createVariable('lat_bin_vals', 'f4', ('n_lat_bins', ))
        self.outfile['lon_bin_vals'][:] = self.lon_bin_vals
        self.outfile['lat_bin_vals'][:] = self.lat_bin_vals

        print('Closing: ' + self.scenario_file)
        self.outfile.close()

    def scenario_initialization(self):
        # Parameterization for sites in regular grid
        self.n_part = 3600   # number of particles & sites initialized
        self.radius = 0  # radius of initialized particles in metres (zero in case of regular grid)
        self.bin_res = 0.06  # resolution of bins (lat and lon) for analysis
        if self.key == "SG800":
            self.get_SG_bounds()
        else:
            print('WARNING: missing key configuration in get_scenario')

        # Define trajectory output file
        self.trajectory_file = (self.outfile_path + self.key + '_' + str(self.year) + '_R'
                                    + str(self.release_n) + '_trajectory.nc')  # Trajectory file path and name

        # Parameterize a grid for analysis:
        self.lat_bin_vals = np.arange(self.domain_lat_min, self.domain_lat_max, self.bin_res)
        self.lon_bin_vals = np.arange(self.domain_lon_min, self.domain_lon_max, self.bin_res)
        self.shp_lon_bins = np.shape(self.lon_bin_vals)[0]
        self.shp_lat_bins = np.shape(self.lat_bin_vals)[0]

        # Set up site positions based on scenario parameterization
        step_lon = (-1 * (self.site_lon_min - self.site_lon_max) / np.sqrt(self.n_part))
        step_lat = (-1 * (self.site_lat_min - self.site_lat_max) / np.sqrt(self.n_part))
        lons = np.arange(self.site_lon_min, self.site_lon_max, step_lon)
        lats = np.arange(self.site_lat_min, self.site_lat_max, step_lat)
        self.site_lat_init, self.site_lon_init = np.meshgrid(lats, lons)  # Initialize lonlat coordinates for sites
        return


    def get_SG_bounds(self):
        self.description = " Initialize with regular grid in South Georgia 800m domain"
        self.site_lon_min = -37
        self.site_lon_max = -36
        self.site_lat_min = -56
        self.site_lat_max = -55.2
        self.domain_lon_min = -42.48
        self.domain_lon_max = -30.64
        self.domain_lat_min = -57.19
        self.domain_lat_max = -50.92
        return
