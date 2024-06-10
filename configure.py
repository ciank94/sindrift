from datetime import timedelta
import numpy as np
import netCDF4 as nc
import sys

class FileExplorer:
    def __init__(self, node, model_name):
        self.node = node
        self.model_name = model_name
        self.local_drift_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/'
        self.mounted_remote_drift_path = 'A:/Cian_sinmod/opendrift/'
        self.remote_drift_path = '/cluster/projects/nn9828k/Cian_sinmod/opendrift/'
        self.init_paths()
        return

    def init_paths(self):
        if self.node == 'local':
            sys.path.insert(0, 'C:/Users/ciank/PycharmProjects/sinmod/opendrift')  # add opendrift local path
            self.phys_states_path = self.mounted_remote_drift_path + 'phys_states/'
            self.trajectory_path = self.local_drift_path + 'trajectory/'
            self.analysis_path = self.local_drift_path + 'analysis/'
            self.figures_path = self.local_drift_path + 'figures/'
        elif self.node == 'remote':
            self.phys_states_path = self.remote_drift_path + 'phys_states/'
            self.trajectory_path = self.remote_drift_path + 'trajectory/'
            self.analysis_path = self.remote_drift_path + 'analysis/'
            self.figures_path = self.remote_drift_path + 'figures/'
        else:
            sys.exit('Specify the correct node in FileExplorer, node = ' + self.node + ' is not an option')

        if self.model_name == "sinmod":
            self.phys_states_file_prefix = "samplesNSEW_" # File identifier
        elif self.model_name == "cmems":
            self.phys_states_file_prefix = 'CMEMS_GLPHYS_D_full_'
        else:
            sys.exit('Specify the correct model_name in FileExplorer, model_name = ' + self.model_name + ' is not an option')

        return

    def search_path(self):
        import os
        path_list = os.listdir(self.analysis_path)
        self.file_list = []
        for file_name in path_list:
            if file_name[0:4] == "SG8H":
                self.file_list.append(file_name)
        return


class Scenario:
    def __init__(self, fpath, reader_phys_states, duration_days, time_step_hours,
                 save_time_step_hours, release_step, y_i, r_i, key):
        # Simulation settings for scenario
        self.key = key  # key used to identify initialization scenario
        self.year = y_i  # simulation year
        self.release_n = r_i + 1  # release number starts at one
        self.z = 15  # release depth
        self.duration = timedelta(hours=24 * duration_days)  # simulation duration in hours as datetime object
        self.time_step = timedelta(hours=time_step_hours)  # simulation time step in hours as datetime object
        self.save_time_step = timedelta(hours=save_time_step_hours)  # how often the file is saved
        self.release_step = release_step  # number of hours between releases
        self.export_variables = ['lon', 'lat']  # choose variables to export from simulation to nc file

        # Information from the phys_states_file:
        self.phys_states_file = reader_phys_states.name
        self.phys_states_timestep = reader_phys_states.time_step.seconds
        self.reader_t_init = reader_phys_states.start_time
        self.t_init = self.reader_t_init + r_i * timedelta(hours=self.release_step)
        self.reader_end_time = reader_phys_states.end_time

        # Initialize scenario parameters and save them as attributes to scenario file
        self.scenario_initialization()  # furnish initialization scenario with class attributes at beginning
        self.trajectory_file = (fpath.trajectory_path + self.key + '_' + str(self.year) + '_R'
                                + str(self.release_n) + '_trajectory.nc')  # Trajectory output file name
        self.analysis_file = (fpath.analysis_path + self.key + '_' + str(self.year) + '_R'
                         + str(self.release_n) + '_trajectory_analysis.nc')  # Analysis output file path name
        self.init_scenario_netcdf(fpath)
        return

    def init_scenario_netcdf(self, fpath):
        self.outfile = nc.Dataset(self.analysis_file, 'w')

        # main simulation settings as attributes:
        self.outfile.title = 'OpenDrift trajectory analysis'
        self.outfile.ocean_model = fpath.model_name
        self.outfile.server = fpath.node
        self.outfile.trajectory_file_prefix = self.key + '_' + str(self.year) + '_R' + str(self.release_n) + '_'
        self.outfile.phys_states_file_prefix = fpath.phys_states_file_prefix
        self.outfile.phys_states_path = fpath.phys_states_path
        self.outfile.phys_states_file = self.phys_states_file
        self.outfile.trajectory_path = fpath.trajectory_path
        self.outfile.trajectory_file = self.trajectory_file
        self.outfile.analysis_path = fpath.analysis_path
        self.outfile.analysis_file = self.analysis_file
        self.outfile.figures_path = fpath.figures_path
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

        print('Closing: ' + self.analysis_file)
        self.outfile.close()
        return

    def scenario_initialization(self):
        # Parameterization for sites in regular grid
        self.n_part = 10000   # number of particles & sites initialized
        self.radius = 0  # radius of initialized particles in metres (zero in case of regular grid)
        self.bin_res = 0.06  # resolution of bins (lat and lon) for analysis
        if self.key == "SG8H":
            self.get_SG_bounds()
        else:
            print('WARNING: missing key configuration in get_scenario')
            exit()

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
        self.site_lon_min = -38
        self.site_lon_max = -35
        self.site_lat_min = -57
        self.site_lat_max = -55.2
        self.domain_lon_min = -42.48
        self.domain_lon_max = -30.64
        self.domain_lat_min = -57.19
        self.domain_lat_max = -50.92
        return




