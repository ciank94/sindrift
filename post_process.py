import numpy as np
import netCDF4 as nc


class Process:

    def __init__(self, trajectory_file, outfile_path,  y_i, r_i, key, test):
        # Organise file data
        self.lon_id = None
        self.lat_id = None
        self.lat_t = None
        self.lon_t = None
        self.bin_res = 0.06
        self.outfile_path = outfile_path
        self.key = key
        self.year = y_i
        self.release_n = r_i + 1
        self.nc_file = nc.Dataset(trajectory_file)

        if test:
            self.p_stride = 10
            self.t_stride = 10
        else:
            self.p_stride = 1
            self.t_stride = 1

        # extract data from trajectory file
        self.lat = self.nc_file['lat']
        self.lon = self.nc_file['lon']
        self.shp_p = np.shape(self.lat)[0]
        self.shp_t = np.shape(self.lat)[1]
        #self.z = self.nc_file['z']

        # Initialize settings for dominant pathway calculation:
        self.min_lon = -75
        self.max_lon = -30
        self.min_lat = -70
        self.max_lat = -50
        self.lat_range = np.arange(self.min_lat - 20, self.max_lat + 15, self.bin_res)
        self.lon_range = np.arange(self.min_lon - 20, self.max_lon + 15, self.bin_res)
        self.shp_lon_bins = np.shape(self.lon_range)[0]
        self.shp_lat_bins = np.shape(self.lat_range)[0]
        self.dom_matrix = np.zeros([self.shp_lon_bins, self.shp_lat_bins])
        self.dom_file = (self.outfile_path + self.key + '_' + str(self.year) + '_R'
                         + str(self.release_n) + '_dominant_paths.npy')
        self.tp_file = (self.outfile_path + self.key + '_' + str(self.year) + '_R'
                         + str(self.release_n) + '_trajectory_tp.nc')
        self.dominant_paths()
        return

    def dominant_paths(self):
        #self.init_ncfile()
        for p_i in range(0, self.shp_p, self.p_stride):
            lon_p = self.lon[p_i, :]
            lat_p = self.lat[p_i, :]
            #self.write_ncfile(lon_p, lat_p, p_i)
            temp_points = np.ones([self.shp_t, 2])*-1
            self.print_pstep(descriptor="Dominant path calculation", p=p_i)

            for t in range(0, self.shp_t, self.t_stride):
                self.lon_t = lon_p[t]
                self.lat_t = lat_p[t]
                self.get_closest_point()
                temp_points[t, 0] = self.lon_id
                temp_points[t, 1] = self.lat_id

            # Then, look at the unique lat-long ids for particle p_i
            unique_rows = np.unique(temp_points, axis=0).astype(int)
            unique_rows = unique_rows[unique_rows[:, 0] > -1]
            if np.shape(unique_rows)[0] > 0:
                self.dom_matrix[unique_rows[:, 0], unique_rows[:, 1]] = self.dom_matrix[
                                                                              unique_rows[:, 0], unique_rows[:, 1]] + 1
        #print('Closing transposed file')
        #self.outfile.close()

        print('Saving: ' + self.dom_file)
        np.save(self.dom_file, self.dom_matrix)
        return

    def write_ncfile(self, lon_p, lat_p, p_i):
        self.outfile.variables['lat'][:, p_i] = lat_p.T
        self.outfile.variables['lon'][:, p_i] = lon_p.T
        return

    def print_pstep(self, descriptor, p):
        if p == 0:
            print("First step in calculation")
        else:
            print(descriptor + " : " + str(p + 1) + " of " + str(self.shp_p))
        if p == self.shp_p:
            print("Final particle calculation")
        return

    def get_closest_point(self):
        self.lon_id = np.argmin(np.sqrt((self.lon_t - self.lon_range[:]) ** 2))
        self.lat_id = np.argmin(np.sqrt((self.lat_t - self.lat_range[:]) ** 2))
        return

    def init_ncfile(self):
        self.outfile = nc.Dataset(self.tp_file, 'w')
        self.outfile.createDimension('time', self.shp_t)
        self.outfile.createDimension('particles', self.shp_p)
        self.outfile.createVariable('lat', 'f4', ('time', 'particles'))
        self.outfile.createVariable('lon', 'f4', ('time', 'particles'))
        return



