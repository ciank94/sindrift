from datetime import timedelta
import numpy as np


class Case:
    def __init__(self, reader, path, time_step_hours, days, y_i):
        self.trajectory_file = None
        self.description = None
        self.lon_init = None
        self.lat_init = None
        self.name = None
        self.year = y_i
        self.n_part = 10000
        self.radius = 10000
        self.z = 50
        self.duration = timedelta(hours=24*days)
        self.time_step = timedelta(hours=time_step_hours)
        self.export_variables = ['lon', 'lat']
        self.path = path
        if time_step_hours < 0:
            self.t_init = reader.end_time
        else:
            self.t_init = reader.start_time
        return

    def get_scenarios(self, key):
        self.name = key
        self.trajectory_file = self.path + key + '_' + str(self.year) + '_trajectory.nc'
        if self.name == "APSO":
            lon_min = -70
            lon_max = -40
            step_lon = (-1 * (lon_min - lon_max) / np.sqrt(self.n_part))
            lat_min = -68
            lat_max = -57
            step_lat = (-1 * (lat_min - lat_max) / np.sqrt(self.n_part))
            lons = np.arange(lon_min, lon_max, step_lon)
            lats = np.arange(lat_min, lat_max, step_lat)
            self.lat_init, self.lon_init = np.meshgrid(lats, lons)
            self.radius = 0
        elif self.name == "SG_NE":
            self.description = "Important fishing ground at NE"
            self.lat_init = -53.8
            self.lon_init = -36

        elif self.name == "SG_NW":
            self.description = "Important fishing ground at NW"
            self.lat_init = -53.75
            self.lon_init = -38.5
        else:
            print('missing key configuration in get_config_params')
        return
