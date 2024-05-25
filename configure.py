from datetime import timedelta


class Case:
    def __init__(self, reader, path, time_step_hours, days):
        self.trajectory_file = None
        self.description = None
        self.lon_init = None
        self.lat_init = None
        self.name = None
        self.n_part = 10000
        self.radius = 10000
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
        self.trajectory_file = self.path + key + '_case_trajectory.nc'
        if self.name == "SG_NE":
            self.description = "Important fishing ground at NE"
            self.lat_init = -53.8
            self.lon_init = -36

        if self.name == "SG_NW":
            self.description = "Important fishing ground at NW"
            self.lat_init = -53.75
            self.lon_init = -38.5
        else:
            print('missing key configuration in get_config_params')
        return
