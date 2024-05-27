from datetime import timedelta
import numpy as np


class Scenario:
    def __init__(self, infile_path, outfile_path, reader, duration_days, time_step_hours, save_time_step_hours,
                 release_step, y_i, r_i, key):
        self.trajectory_file = None  # Trajectory file path and name
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
        self.scenario_initialization()  # furnish initialization scenario with class attributes at beginning
        return

    def scenario_initialization(self):
        # Define trajectory output file
        self.trajectory_file = (self.outfile_path + self.key + '_' + str(self.year) + '_R'
                                + str(self.release_n) + '_trajectory.nc')
        if self.key == "APSO":
            self.description = " Initialize with regular grid covering the Antarctic Peninsula and South Orkney Islands"
            self.n_part = 10000
            self.radius = 0
            lon_min = -70
            lon_max = -40
            lat_min = -68
            lat_max = -57
            step_lon = (-1 * (lon_min - lon_max) / np.sqrt(self.n_part))
            step_lat = (-1 * (lat_min - lat_max) / np.sqrt(self.n_part))
            lons = np.arange(lon_min, lon_max, step_lon)
            lats = np.arange(lat_min, lat_max, step_lat)
            self.lat_init, self.lon_init = np.meshgrid(lats, lons)
        else:
            print('WARNING: missing key configuration in get_scenario')
        return
