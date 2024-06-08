from datetime import timedelta
import numpy as np
import netCDF4 as nc


class Scenario:
    def __init__(self, infile_path, outfile_path, reader, duration_days, time_step_hours, save_time_step_hours,
                 release_step, y_i, r_i, key):
        self.lon_init = None  # Initialization longitude(s)
        self.lat_init = None  # Initialization latitude(s)
        self.n_part = None  # number of particles initialized
        self.radius = None  # radius of initialized particles in metres (zero in case of regular grid)
        self.description = None  # Case description (todo: write to nc file as metadata)
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
        if time_step_hours < 0:  # for backward simulation start with final time on nc file
            self.t_init = reader.end_time
        else:
            self.t_init = reader.start_time
            self.t_init = self.t_init + r_i * timedelta(hours=self.release_step)

        #Initialize scenario parameters and save them as attributes to scenario file
        self.scenario_initialization()  # furnish initialization scenario with class attributes at beginning
        self.scenario_file = (self.outfile_path + self.key + '_' + str(self.year) + '_R'
                         + str(self.release_n) + '_trajectory_scenario.nc')
        return

    def init_scenario_netcdf(self):
        self.outfile = nc.Dataset(self.scenario_file, 'w')
        self.outfile.createDimension('lon_bins', self.shp_lon_bins)
        self.outfile.createDimension('lat_bins', self.shp_lat_bins)
        self.outfile.createVariable('dom_paths', 'i4', ('lon_bins', 'lat_bins'))
        self.outfile.variables['dom_paths'][:] = 0

    def scenario_initialization(self):
        # Parameterization for sites in regular grid
        self.n_part = 3600
        self.radius = 0
        self.bin_res = 0.06


        if self.key == "APSO":
            self.description = " Initialize with regular grid covering the Antarctic Peninsula and South Orkney Islands"
            self.site_lon_min = -70
            self.site_lon_max = -40
            self.site_lat_min = -68
            self.site_lat_max = -57

        elif self.key == "SG800":
            self.description = " Initialize with regular grid in South Georgia 800m domain"
            self.site_lon_min = -37
            self.site_lon_max = -36
            self.site_lat_min = -56
            self.site_lat_max = -55.2
            #todo: make function get SG site grid and analysis grid- use preformed values- based on matdisp boundaries;

        else:
            print('WARNING: missing key configuration in get_scenario')

        # Define trajectory output file
        self.trajectory_file = (self.outfile_path + self.key + '_' + str(self.year) + '_R'
                                    + str(self.release_n) + '_trajectory.nc')  # Trajectory file path and name

        # Parameterize a grid for analysis:
        self.lat_range = np.arange(self.site_lat_min - 20, self.site_lat_max + 15, self.bin_res)
        self.lon_range = np.arange(self.min_lon - 20, self.max_lon + 15, self.bin_res)
        self.shp_lon_bins = np.shape(self.lon_range)[0]
        self.shp_lat_bins = np.shape(self.lat_range)[0]

        # Set up site positions based on scenario parameterization
        step_lon = (-1 * (self.site_lon_min - self.site_lon_max) / np.sqrt(self.n_part))
        step_lat = (-1 * (self.site_lat_min - self.site_lat_max) / np.sqrt(self.n_part))
        lons = np.arange(self.site_lon_min, self.site_lon_max, step_lon)
        lats = np.arange(self.site_lat_min, self.site_lat_max, step_lat)
        self.site_lat_init, self.site_lon_init = np.meshgrid(lats, lons)
        return
