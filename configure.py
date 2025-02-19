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
        self.local_drift_path = "C:/Users/ciank/PycharmProjects/sinmod/sindrift/"
        self.mounted_remote_drift_path = "A:/Cian_sinmod/opendrift/"
        self.remote_drift_path = "/cluster/projects/nn9828k/Cian_sinmod/opendrift/"
        self.key = key
        self.init_paths()
        return

    def init_paths(self):
        valid_key_list = ["SG8H", "BSSI", "SOIN", "SGCM"]
        if self.key not in valid_key_list:
            sys.exit("Not a valid key pointing to an initialisation scenario")
        # todo: use dictionary with folder names
        if self.node == "local":
            sys.path.insert(
                0, "C:/Users/ciank/PycharmProjects/sinmod/opendrift"
            )  # add opendrift local path
            self.phys_states_path = self.mounted_remote_drift_path + "phys_states/"
            self.trajectory_path = self.local_drift_path + "trajectory/"
            self.analysis_path = self.local_drift_path + "analysis/"
            self.compile_path = self.local_drift_path + "compile/"
            self.figures_path = self.local_drift_path + "figures/"
        elif self.node == "remote":
            self.phys_states_path = self.remote_drift_path + "phys_states/"
            self.trajectory_path = self.remote_drift_path + "trajectory/"
            self.analysis_path = self.remote_drift_path + "analysis/"
            self.compile_path = self.remote_drift_path + "compile/"
            self.figures_path = self.remote_drift_path + "figures/"
        else:
            sys.exit(
                "Specify the correct node in FileExplorer, node = "
                + self.node
                + " is not an option"
            )

        if self.model_name == "sinmod":
            self.phys_states_file_prefix = "samplesNSEW_"  # File identifier
        elif self.model_name == "cmems":
            self.phys_states_file_prefix = "CMEMS_GLPHYS_D_full_"
        else:
            sys.exit(
                "Specify the correct model_name in FileExplorer, model_name = "
                + self.model_name
                + " is not an option"
            )
        return

    def mounted_paths(self):
        self.analysis_path = self.mounted_remote_drift_path + "analysis/"
        self.trajectory_path = self.mounted_remote_drift_path + "trajectory/"
        self.figures_path = self.local_drift_path + "figures/"
        return

    def local_phys_states(self):
        self.phys_states_path = self.local_drift_path + "phys_states/"
        return

    def get_phys_states(self, date_init, date_limit):
        if self.model_name == "sinmod":
            if not date_init.month == date_limit.month:
                phys_states_list = []
                for kMonth in range(date_init.month, date_limit.month + 1):
                    phys_states_i = (
                        self.phys_states_path
                        + self.phys_states_file_prefix
                        + str(date_init.year)
                        + "{:02d}".format(kMonth)
                        + ".nc"
                    )
                    phys_states_list.append(phys_states_i)

                phys_states = xr.open_mfdataset(phys_states_list)
            else:
                phys_states = (
                    self.phys_states_path
                    + self.phys_states_file_prefix
                    + str(date_init.year)
                    + "{:02d}".format(date_init.month)
                    + ".nc"
                )
        elif self.model_name == "cmems":
            if not date_init.year == date_limit.year:
                phys_states_list = []
                for y in range(date_init.year, date_limit.year + 1):
                    phys_states_i = (
                        self.phys_states_path
                        + self.phys_states_file_prefix
                        + str(y)
                        + ".nc"
                    )
                    phys_states_list.append(phys_states_i)
                phys_states = xr.open_mfdataset(phys_states_list)
            else:
                phys_states = (
                    self.phys_states_path
                    + self.phys_states_file_prefix
                    + str(date_init.year)
                    + ".nc"
                )
        else:
            sys.exit("incorrect model name in check_phys_states")
        return phys_states

    def search_path(self, year, release_start, release_end):
        path_list = os.listdir(self.analysis_path)
        self.file_list = []
        release_range = np.arange(release_start, release_end + 1, 1)
        print("===========================")
        print("File list for analysis: ")
        print("---------------------------")
        for file_name in path_list:
            if file_name[0:4] == self.key:
                if file_name[5:9] == str(year):
                    if file_name[12:13] == "_":
                        if int(file_name[11:12]) in release_range:
                            self.file_list.append(file_name)
                            print(file_name)
                    else:
                        if int(file_name[11:13]) in release_range:
                            self.file_list.append(file_name)
                            print(file_name)

        print("===========================")
        return


class Scenario:
    def __init__(
        self,
        fpath,
        date_init,
        date_release,
        duration_days,
        time_step,
        save_step,
        release_step,
        release_i,
        reader_phys_states,
    ):
        # Simulation settings for scenario
        self.key = fpath.key  # key used to identify initialization scenario
        self.year = date_release.year  # simulation year
        self.release_n = release_i + 1  # release number starts at one
        self.duration = duration_days  # simulation duration in hours as datetime object
        self.time_step = time_step  # simulation time step in hours as datetime object
        self.save_time_step = save_step  # how often the file is saved
        self.release_step = release_step  # number of hours between releases
        self.export_variables = [
            "lon",
            "lat",
            "z",
        ]  # choose variables to export from simulation to nc file

        # Information from the phys_states_file:
        self.phys_states_file = reader_phys_states.name
        self.phys_states_timestep = reader_phys_states.time_step.seconds
        self.reader_date_init = reader_phys_states.start_time
        self.reader_date_end = reader_phys_states.end_time
        self.date_init = date_release
        self.date_end = date_release + duration_days

        # Initialize scenario parameters and save them as attributes to scenario file
        self.scenario_initialization()  # furnish initialization scenario with class attributes at beginning
        self.trajectory_file_name = (
            self.key
            + "_"
            + str(date_init.year)
            + "_R"
            + str(self.release_n)
            + "_trajectory.nc"
        )
        self.trajectory_file = (
            fpath.trajectory_path + self.trajectory_file_name
        )  # Trajectory output file name
        self.analysis_file_name = (
            self.key
            + "_"
            + str(self.year)
            + "_R"
            + str(self.release_n)
            + "_trajectory_analysis.nc"
        )
        self.analysis_file = (
            fpath.analysis_path + self.analysis_file_name
        )  # Analysis output file path name
        self.init_scenario_netcdf(fpath)
        return

    def scenario_initialization(self):
        # Parameterization for sites in regular grid
        self.n_part = 10000  # number of particles & sites initialized
        self.radius = 0  # radius of initialized particles in metres (zero in case of regular grid)
        self.bin_res = 0.2  # resolution of bins (lat and lon) for analysis
        if self.key == "SG8H":
            self.get_NEMO_bounds()
            self.get_SG8H_bounds()
        elif self.key == "SOIN":
            self.get_NEMO_bounds()
            self.get_SOIN_bounds()
        elif self.key == "BSSI":
            self.get_NEMO_bounds()
            self.get_BSSI_bounds()
        elif self.key == "SGCM":
            self.get_SGret_bounds()
            self.get_SGCM_bounds()
        else:
            sys.exit("WARNING: missing key configuration in get_scenario")

        # Parameterize a grid for analysis:
        self.lat_bin_vals = np.arange(
            self.domain_lat_min, self.domain_lat_max, self.bin_res
        )
        self.lon_bin_vals = np.arange(
            self.domain_lon_min, self.domain_lon_max, self.bin_res
        )
        self.shp_lon_bins = np.shape(self.lon_bin_vals)[0]
        self.shp_lat_bins = np.shape(self.lat_bin_vals)[0]

        # Set up site positions based on scenario parameterization
        step_lon = -1 * (self.site_lon_min - self.site_lon_max) / np.sqrt(self.n_part)
        step_lat = -1 * (self.site_lat_min - self.site_lat_max) / np.sqrt(self.n_part)
        lons = np.arange(self.site_lon_min, self.site_lon_max, step_lon)
        lats = np.arange(self.site_lat_min, self.site_lat_max, step_lat)
        self.site_lat_init, self.site_lon_init = np.meshgrid(
            lats, lons
        )  # Initialize lonlat coordinates for sites
        return

    def get_NEMO_bounds(self):
        self.domain_lon_min = -69
        self.domain_lon_max = -30
        self.domain_lat_min = -65
        self.domain_lat_max = -48
        self.bound_name = "NEMO"
        return

    def get_SGret_bounds(self):
        self.domain_lon_min = -69
        self.domain_lon_max = -30
        self.domain_lat_min = -65
        self.domain_lat_max = -48
        self.bound_name = "SGret"
        return

    def get_BSSI_bounds(self):
        self.description = (
            " Initialize with regular grid in Antarctic Peninsula nemo domain"
        )
        self.site_lon_min = -63
        self.site_lon_max = -54
        self.site_lat_min = -65
        self.site_lat_max = -60.5
        self.z = 75  # release depth
        return

    def get_SOIN_bounds(self):
        self.description = (
            " Initialize with regular grid in South Orkney islands nemo domain"
        )
        self.site_lon_min = -48
        self.site_lon_max = -44.35
        self.site_lat_min = -61
        self.site_lat_max = -59.2
        self.z = 75  # release depth
        return

    def get_SG8H_bounds(self):
        self.description = " Initialize with regular grid in South Georgia 800m domain"
        self.site_lon_min = -42
        self.site_lon_max = -37.5
        self.site_lat_min = -57
        self.site_lat_max = -54
        self.z = 75  # release depth
        return

    def get_SGCM_bounds(self):
        self.description = " Initialize with regular grid in South Georgia nemo domain"
        self.site_lon_min = -39.4
        self.site_lon_max = -35
        self.site_lat_min = -55
        self.site_lat_max = -53.2
        self.z = 75  # release depth
        return

    def init_scenario_netcdf(self, fpath):
        self.outfile = nc.Dataset(self.analysis_file, "w")

        # main simulation settings as attributes:
        self.outfile.title = "OpenDrift trajectory analysis"
        self.outfile.ocean_model = fpath.model_name
        self.outfile.server = fpath.node
        self.outfile.trajectory_file_prefix = (
            self.key + "_" + str(self.year) + "_R" + str(self.release_n) + "_"
        )
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
        self.outfile.release_step_days = (
            self.release_step.days
        )  # number of hours between releases

        # simulation time information
        self.outfile.sim_start_day = self.date_init.day  # simulation start time
        self.outfile.sim_start_month = self.date_init.month  # simulation start time
        self.outfile.sim_start_year = self.date_init.year  # simulation year
        self.outfile.sim_end_day = self.date_end.day
        self.outfile.sim_end_month = self.date_end.month  # simulation end time
        self.outfile.sim_end_year = self.date_end.year
        self.outfile.sim_duration_days = self.duration.days  # simulation duration
        self.outfile.sim_time_step_seconds = (
            self.time_step.seconds
        )  # simulation time step
        self.outfile.sim_save_time_step_seconds = (
            self.save_time_step.seconds
        )  # simulation save time step

        # phys_states reader information
        self.outfile.reader_time_step_seconds = (
            self.phys_states_timestep
        )  # phys_states time step
        self.outfile.export_variables = (
            self.export_variables
        )  # variables exported from simulation

        # domain and site limits:
        self.outfile.site_lon_min = self.site_lon_min
        self.outfile.site_lon_max = self.site_lon_max
        self.outfile.site_lat_min = self.site_lat_min
        self.outfile.site_lat_max = self.site_lat_max
        self.outfile.bounds = self.bound_name
        self.outfile.domain_lon_min = self.domain_lon_min
        self.outfile.domain_lon_max = self.domain_lon_max
        self.outfile.domain_lat_min = self.domain_lat_min
        self.outfile.domain_lat_max = self.domain_lat_max
        self.outfile.bin_resolution = self.bin_res

        # set_dimensions of file:
        self.outfile.createDimension("n_lon_bins", self.shp_lon_bins)
        self.outfile.createDimension("n_lat_bins", self.shp_lat_bins)
        self.outfile.createVariable("lon_bin_vals", "f4", ("n_lon_bins",))
        self.outfile.createVariable("lat_bin_vals", "f4", ("n_lat_bins",))
        self.outfile["lon_bin_vals"][:] = self.lon_bin_vals
        self.outfile["lat_bin_vals"][:] = self.lat_bin_vals

        print("Closing: " + self.analysis_file)
        self.outfile.close()
        return
