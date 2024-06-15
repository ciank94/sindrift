import os
import os.path
import numpy as np
import netCDF4 as nc
import xarray as xr
import sys


class FileExplorer:
    def __init__(self, node, model_name, key):
        self.node = node
        self.model_name = model_name
        self.local_drift_path = 'C:/Users/ciank/PycharmProjects/sinmod/sindrift/'
        self.mounted_remote_drift_path = 'A:/Cian_sinmod/opendrift/'
        self.remote_drift_path = '/cluster/projects/nn9828k/Cian_sinmod/opendrift/'
        self.key = key
        self.init_paths()
        return

    def init_paths(self):
        valid_key_list = ["SG8H", "APSO"]
        if self.key not in valid_key_list:
            sys.exit('Not a valid key pointing to an initialisation scenario')
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

    def mounted_paths(self):
        self.analysis_path = self.mounted_remote_drift_path + 'analysis/'
        self.trajectory_path = self.mounted_remote_drift_path + 'trajectory/'
        self.figures_path = self.local_drift_path + 'figures/'
        return


    def get_phys_states(self, date_init, date_limit):
        if self.model_name == "sinmod":
            if not date_init.month == date_limit.month:
                phys_states_list = []
                for kMonth in range(date_init.month, date_limit.month+1):
                    phys_states_i = (self.phys_states_path + self.phys_states_file_prefix + str(date_init.year) +
                                   "{:02d}".format(kMonth) + '.nc')
                    phys_states_list.append(phys_states_i)

                phys_states = xr.open_mfdataset(phys_states_list)
            else:
                phys_states = (self.phys_states_path + self.phys_states_file_prefix + str(date_init.year) +
                               "{:02d}".format(date_init.month) + '.nc')
        elif self.model_name == "cmems":
            if not date_init.year == date_limit.year:
                phys_states_list = []
                for y in range(date_init.year, date_limit.year + 1):
                    phys_states_i = (self.phys_states_path + self.phys_states_file_prefix + str(y) + '.nc')
                    phys_states_list.append(phys_states_i)
                phys_states = xr.open_mfdataset(phys_states_list)
            else:
                phys_states = self.phys_states_path + self.phys_states_file_prefix + str(date_init.year) + '.nc'
        else:
            sys.exit('incorrect model name in check_phys_states')
        return phys_states

    def search_path(self, year):
        path_list = os.listdir(self.analysis_path)
        self.file_list = []
        print('===========================')
        print('File list for analysis: ')
        print('---------------------------')
        for file_name in path_list:
            if file_name[0:4] == self.key:
                if file_name[5:9] == str(year):
                    self.file_list.append(file_name)
                    print(file_name)
        print('===========================')
        return

class Scenario:
    def __init__(self, fpath, date_release, duration_days, time_step, save_step, release_step, release_i,
                 reader_phys_states):
        # Simulation settings for scenario
        self.key = fpath.key  # key used to identify initialization scenario
        self.year = date_release.year  # simulation year
        self.release_n = release_i + 1  # release number starts at one
        self.duration = duration_days  # simulation duration in hours as datetime object
        self.time_step = time_step  # simulation time step in hours as datetime object
        self.save_time_step = save_step  # how often the file is saved
        self.release_step = release_step  # number of hours between releases
        self.export_variables = ['lon', 'lat']  # choose variables to export from simulation to nc file

        # Information from the phys_states_file:
        self.phys_states_file = reader_phys_states.name
        self.phys_states_timestep = reader_phys_states.time_step.seconds
        self.reader_date_init = reader_phys_states.start_time
        self.reader_date_end = reader_phys_states.end_time
        self.date_init = date_release
        self.date_end = date_release + duration_days

        # Initialize scenario parameters and save them as attributes to scenario file
        self.scenario_initialization()  # furnish initialization scenario with class attributes at beginning
        self.trajectory_file_name = self.key + '_' + str(self.year) + '_R' + str(self.release_n) + '_trajectory.nc'
        self.trajectory_file = (fpath.trajectory_path + self.trajectory_file_name)  # Trajectory output file name
        self.analysis_file_name = self.key + '_' + str(self.year) + '_R' + str(self.release_n) + '_trajectory_analysis.nc'
        self.analysis_file = (fpath.analysis_path + self.analysis_file_name)  # Analysis output file path name
        self.init_scenario_netcdf(fpath)
        return


    def scenario_initialization(self):
        # Parameterization for sites in regular grid
        self.n_part = 10000   # number of particles & sites initialized
        self.radius = 0  # radius of initialized particles in metres (zero in case of regular grid)
        self.bin_res = 0.2  # resolution of bins (lat and lon) for analysis
        if self.key == "SG8H":
            self.get_SG_bounds()
        elif self.key == "APSO":
            self.get_AP_bounds()
        else:
            sys.exit('WARNING: missing key configuration in get_scenario')

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


    def get_AP_bounds(self):
        self.description = " Initialize with regular grid in Antarctic Peninsula nemo domain"
        self.site_lon_min = -67
        self.site_lon_max = -40
        self.site_lat_min = -69
        self.site_lat_max = -57
        # self.site_lon_max = -50
        # self.site_lat_min = -60
        # self.site_lat_max = -58
        self.domain_lon_min = -69
        self.domain_lon_max = -30
        self.domain_lat_min = -65
        self.domain_lat_max = -48
        self.z = 50  # release depth
        self.APSO_square_polyons()
        return

    def SSMUs(self, fpath):
        import geopandas as gpd
        shape_p = gpd.read_file(fpath.figures_path + "ssmusPolygon.shp")
        g1 = shape_p[0].geometry
        return


    def get_SG_bounds(self):
        self.description = " Initialize with regular grid in South Georgia 800m domain"
        self.site_lon_min = -42
        self.site_lon_max = -37.5
        self.site_lat_min = -57
        self.site_lat_max = -54
        self.domain_lon_min = -42.48
        self.domain_lon_max = -30.64
        self.domain_lat_min = -57.19
        self.domain_lat_max = -50.92
        self.z = 15  # release depth
        self.SG_square_polygons()
        return

    def APSO_square_polyons(self):
        self.n_polys = 2
        self.lon_lims_poly = np.zeros([self.n_polys, 2])
        self.lat_lims_poly = np.zeros([self.n_polys, 2])

        # Define area:
        lon_lims = [-39.5, -35]
        lat_lims = [-54, -53]
        poly_desc = 'polygon number 1: SG northern part of island'
        self.assign_polygons(lon_lims, lat_lims, poly_desc, i=0)

        # Define area:
        lon_lims = [-42.48, -30.64]
        lat_lims = [-57.19, -50.92]
        poly_desc = 'polygon number 2: SG domain boundaries (800m model)'
        self.assign_polygons(lon_lims, lat_lims, poly_desc, i=1)
        return



    def SG_square_polygons(self):
        import matplotlib.pyplot as plt
        self.n_polys = 3
        self.lon_lims_poly = np.zeros([self.n_polys, 2])
        self.lat_lims_poly = np.zeros([self.n_polys, 2])

        # Define area:
        lon_lims = [-39.25, -37.85]
        lat_lims = [-53.85, -53.55]
        poly_desc = 'polygon number 1: SG hotspot on the north western part of the island'
        self.assign_polygons(lon_lims, lat_lims, poly_desc, i=0)

        # area 2
        lon_lims = [-36.1, -35.1]
        lat_lims = [-54.85, -53.4]
        poly_desc = 'polygon number 2: SG hotspot on the north eastern part of the island'
        self.assign_polygons(lon_lims, lat_lims, poly_desc, i=1)

        # test case
        lon_lims = [-34, -32]
        lat_lims = [-56, -54]
        poly_desc = 'polygon number 3: test case for SG recruitment and retention'
        self.assign_polygons(lon_lims, lat_lims, poly_desc, i=2)
        return

    def assign_polygons(self, lon_lims, lat_lims, poly_desc, i):
        nline_indent = '\n\t\t\t\t\t\t  '
        self.lon_lims_poly[i, 0:2] = lon_lims
        self.lat_lims_poly[i, 0:2] = lat_lims
        if i == 0:
            self.poly_desc = poly_desc
        else:
            self.poly_desc = self.poly_desc + nline_indent + poly_desc
        return


    def init_scenario_netcdf(self, fpath):
        self.outfile = nc.Dataset(self.analysis_file, 'w')

        # main simulation settings as attributes:
        self.outfile.title = 'OpenDrift trajectory analysis'
        self.outfile.ocean_model = fpath.model_name
        self.outfile.server = fpath.node
        self.outfile.trajectory_file_prefix = self.key + '_' + str(self.year) + '_R' + str(self.release_n) + '_'
        self.outfile.phys_states_file_prefix = fpath.phys_states_file_prefix
        self.outfile.trajectory_path = fpath.trajectory_path
        self.outfile.trajectory_file = self.trajectory_file_name
        self.outfile.analysis_path = fpath.analysis_path
        self.outfile.analysis_file = self.analysis_file_name
        self.outfile.figures_path = fpath.figures_path
        self.outfile.scenario_key = self.key  # key name given to scenario
        self.outfile.scenario_description = self.description  # description of scenario

        # information about release:
        self.outfile.n_part = self.n_part
        self.outfile.radius = self.radius
        self.outfile.release_number = self.release_n  # release number starts at one
        self.outfile.release_depth = self.z  # release depth
        self.outfile.release_step_days = self.release_step.days  # number of hours between releases

        # simulation time information
        self.outfile.sim_start_day = self.date_init.day # simulation start time
        self.outfile.sim_start_month = self.date_init.month  # simulation start time
        self.outfile.sim_start_year = self.date_init.year  # simulation year
        self.outfile.sim_end_day = self.date_end.day
        self.outfile.sim_end_month = self.date_end.month  # simulation end time
        self.outfile.sim_end_year = self.date_end.year

        self.outfile.sim_duration_days = self.duration.days # simulation duration
        self.outfile.sim_time_step_seconds = self.time_step.seconds # simulation time step
        self.outfile.sim_save_time_step_seconds = self.save_time_step.seconds # simulation save time step

        # phys_states reader information
        self.outfile.reader_time_step_seconds = self.phys_states_timestep  # phys_states time step
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

        if self.n_polys > 0:
            self.outfile.createDimension('n_polygons', self.n_polys)  # n_polys is defined for the key
            self.outfile.createDimension('points', 4)
            self.outfile.createVariable('square_polygons', 'f4', ('n_polygons', 'points'))
            self.outfile['square_polygons'].coords_order = ('order of coords [0:4]: (lon_1, lon_2, lat_1, lat_2), '
                                                            'use coords = ((lon_1, lat_1), '
                                                            '(lon_1, lat_2), (lon_2, lat_2),'
                                                            '(lon_2, lat_1), (lon_1, lat_1))')
            for n_poly in range(0, self.n_polys):
                self.outfile['square_polygons'][n_poly, 0:2] = self.lon_lims_poly[n_poly, :]
                self.outfile['square_polygons'][n_poly, 2:4] = self.lat_lims_poly[n_poly, :]
                self.outfile['square_polygons'].polygon_descriptions = self.poly_desc

        print('Closing: ' + self.analysis_file)
        self.outfile.close()
        return







