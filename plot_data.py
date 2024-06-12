import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.patches import Polygon
from shapely.geometry.polygon import LinearRing, Polygon
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import netCDF4 as nc
import os

class PlotData:

    def __init__(self, fpath, analysis_file):
        # files for analysis
        self.analysis_file = analysis_file
        self.analysis_df = nc.Dataset(analysis_file, mode='r')
        self.trajectory_df = nc.Dataset(self.analysis_df.trajectory_file, mode='r')
        self.fpath = fpath

        # extract data from trajectory file
        self.load_trajectory_df()

        # extract data from analysis file
        self.load_analysis_df()
        self.save_prefix = (self.analysis_df.scenario_key + '_' + str(self.analysis_df.sim_start_year) + '_R'
                            + str(self.analysis_df.release_number) + '_')

        # For creating the bathymetry file
        self.load_bathymetry()

        # parameters for plotting:
        self.load_plot_params()
        return

    def load_plot_params(self):

        # Dominant pathways color scaling:
        self.d_scale = 40
        self.max_scale = 0.8

        # plotting parameters
        self.bath_contours = np.linspace(0, 3000, 10)
        self.bath_cmap = plt.get_cmap('Blues')
        #self.dom_cmap = plt.get_cmap('Reds')
        self.dom_cmap = plt.get_cmap('OrRd')
        self.depth_colors = np.arange(0, 4500, 200)

        # offsets for figure boundaries:
        self.lon_offset = 1
        self.lat_offset = 3
        return

    def load_analysis_df(self):
        self.lon_bin_vals = self.analysis_df.variables['lon_bin_vals'][:]
        self.lat_bin_vals = self.analysis_df.variables['lat_bin_vals'][:]
        self.bin_res = self.analysis_df.bin_resolution
        self.n_lat_bins = np.shape(self.lat_bin_vals)[0]
        self.n_lon_bins = np.shape(self.lat_bin_vals)[0]
        self.key = self.analysis_df.scenario_key
        self.min_lon = self.analysis_df.domain_lon_min
        self.min_lat = self.analysis_df.domain_lat_min
        self.max_lon = self.analysis_df.domain_lon_max
        self.max_lat = self.analysis_df.domain_lat_max
        self.polygons = self.analysis_df['square_polygons']
        self.recruit = self.analysis_df['recruit'][:]
        self.CG = self.analysis_df['CG'][:]
        return

    def load_trajectory_df(self):
        self.lat = self.trajectory_df['lat']
        self.lon = self.trajectory_df['lon']
        self.shp_p = np.shape(self.lat)[0]
        self.shp_t = np.shape(self.lat)[1]
        # self.z = self.nc_file['z']
        return

    def load_bathymetry(self):
        self.bath_res = 0.04  # bathymetry resolution
        self.bath_file = self.fpath.figures_path + 'bath.npy'
        self.bath_file_lon = self.fpath.figures_path + 'bath_lon.npy'
        self.bath_file_lat = self.fpath.figures_path + 'bath_lat.npy'
        if not os.path.exists(self.bath_file):
            self.get_bathfile()
            print('Creating ' + self.bath_file)
        self.bath = np.load(self.bath_file)
        self.bath_lon = np.load(self.bath_file_lon)
        self.bath_lat = np.load(self.bath_file_lat)
        return

    def get_bathfile(self):
        bathfile_gebco = self.fpath.figures_path + 'gebco_2023.nc'

        # take in the gebco bathymetry file
        bath_f = nc.Dataset(bathfile_gebco)
        e = np.array(bath_f['elevation'][:]).astype(float)
        lat_e = np.array(bath_f['lat'][:])
        lon_e = np.array(bath_f['lon'][:])

        e[e > 0] = np.nan
        e *= -1

        new_lat = np.arange(np.min(lat_e), np.max(lat_e), self.bath_res)
        new_lon = np.arange(np.min(lon_e), np.max(lon_e), self.bath_res)

        shp_lat = np.shape(new_lat)[0]
        shp_lon = np.shape(new_lon)[0]
        e_new = np.zeros([shp_lat, shp_lon])
        for i in range(0, shp_lat):
            for j in range(0, shp_lon):
                lat_id = np.argmin(np.sqrt((lat_e[:] - new_lat[i]) ** 2))
                lon_id = np.argmin(np.sqrt((lon_e[:] - new_lon[j]) ** 2))
                e_new[i, j] = e[lat_id, lon_id]

        np.save(self.bath_file, e_new)
        np.save(self.bath_file_lat, new_lat)
        np.save(self.bath_file_lon, new_lon)
        bath_f.close()
        return


    def plot_dom_paths(self, dom_paths):
        self.init_plot()
        self.plot_background()
        n_levels = np.arange(np.min(dom_paths), np.max(dom_paths), np.max(dom_paths)/20)
        self.plot1 = plt.contourf(self.lon_bin_vals, self.lat_bin_vals, dom_paths.T, levels=n_levels, cmap=self.dom_cmap, transform=ccrs.PlateCarree(), extend='both')
        self.c_max = self.max_scale * np.max(dom_paths)
        self.add_cbar(c_max= np.max(dom_paths), caxis_title='unique_particles')
        plt_name = self.save_prefix + "_dom_paths"
        self.save_plot(plt_name)
        return

    def plot_recruits(self):
        n_poly = 0
        recruit = self.recruit[:,:,n_poly]
        recruit_t = recruit[:,0]
        recruit_i = recruit[:,1]
        id1 = np.where(recruit_t>0)
        self.init_plot()
        self.plot_background()
        i_vals = recruit_i[id1[0][:]].astype(int)
        particle_n = id1[0][:].astype(int)
        x = np.array(self.lon[particle_n, :])
        y = np.array(self.lat[particle_n, :])
        self.ax.scatter(x[:, 0], y[:, 0], alpha=0.7, c='r')
        [self.ax.plot(x[i, 0:i_vals[i]+3], y[i, 0:i_vals[i]+3], alpha=0.1, linewidth=1, c='k') for i in range(0, np.shape(id1)[1])]
        poly = self.get_polygon(0)
        x,y = poly.exterior.xy
        self.ax.plot(x, y, color='#6699cc', alpha=0.7,
                linewidth=3, solid_capstyle='round', zorder=2)

        #self.ax.plot(polygon1)
        #self.ax.add_patch(polygon1, facecolor='r', alpha=0.4)
        plt.show()
        breakpoint()



    def get_polygon(self, poly_n):
        self.lon_1 = self.analysis_df.variables['square_polygons'][poly_n, 0]
        self.lon_2 = self.analysis_df.variables['square_polygons'][poly_n, 1]
        self.lat_1 = self.analysis_df.variables['square_polygons'][poly_n, 2]
        self.lat_2 = self.analysis_df.variables['square_polygons'][poly_n, 3]
        coords = ((self.lon_1, self.lat_1), (self.lon_1, self.lat_2), (self.lon_2, self.lat_2),
                  (self.lon_2, self.lat_1), (self.lon_1, self.lat_1))  # tuple item;
        poly = Polygon(coords)
        return poly

    def plot_trajectory(self):
        self.init_plot()
        self.plot_background()
        self.plot_depth()
        self.add_cbar(caxis_title='depth (m)', c_max=4500)

        step_v = np.floor(np.shape(self.lon)[0] / 1000).astype(int)
        step_v2 = np.floor(np.shape(self.lon)[1] / (np.shape(self.lon)[1] * 0.99)).astype(int)

        lon_1 = self.lon[0:-1:step_v, 0:-1:step_v2]
        lat_1 = self.lat[0:-1:step_v, 0:-1:step_v2]
        c_vals = np.arange(0, np.shape(lat_1)[1])
        c_vals = c_vals * np.ones([np.shape(lat_1)[0], np.shape(lat_1)[1]])
        #

        plt.scatter(lon_1, lat_1, c=c_vals, s=1.3, cmap='YlOrRd', alpha=1, linewidth=1.3, linestyle='-', marker='_')
        plt.scatter(lon_1[:, 0], lat_1[:, 0], s=18, facecolor='yellow', edgecolors='k', alpha=0.9, linewidth=0.6)
        plt.scatter(lon_1[:, -1], lat_1[:, -1], s=18, facecolor='red', edgecolors='k', alpha=0.9, linewidth=0.6)

        plt_name = self.save_prefix + "worms"
        self.save_plot(plt_name)
        return

    def plot_background(self):
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor='face',
                                                facecolor='lightgrey')
        self.ax.add_feature(land_10m)
        self.ax.coastlines(resolution='10m', linewidth=0.7)
        plt.contour(self.bath_lon, self.bath_lat, self.bath, self.bath_contours, colors='k', alpha=0.2, linewidths=0.7,
                    transform=ccrs.PlateCarree())

        # set extent and grid lines;
        gl = self.ax.gridlines(draw_labels=True, alpha=0.4)
        gl.top_labels = False
        gl.right_labels = False
        self.ax.set_extent([self.min_lon-self.lon_offset, self.max_lon+self.lon_offset, self.min_lat-self.lat_offset,
                            self.max_lat+self.lat_offset])
        return

    def plot_depth(self):
        # levels= self.depth_colors
        self.plot1 = plt.contourf(self.bath_lon, self.bath_lat, self.bath, levels=self.depth_colors,
                                  cmap=self.bath_cmap,
                                  transform=ccrs.PlateCarree(), extend='both')
        return

    def init_plot(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(projection=ccrs.PlateCarree())
        return

    def add_cbar(self, c_max, caxis_title):
        cbar = plt.colorbar(self.plot1, extend='both')
        cbar.ax.set_ylabel(caxis_title, loc='center', size=9, weight='bold')
        cbar.ax.tick_params(labelsize=10, rotation=0)
        plt.clim(0, c_max)
        return

    def save_plot(self, plt_name):
        savefile = self.fpath.figures_path + plt_name + '.png'
        print('Saving file: ' + savefile)
        plt.savefig(savefile, dpi=400)
        plt.close()
        return




class FuseData:
    def __init__(self, infile_path, outfile_path, year_ids, release_ids, key_ids):
        self.dom_path_i = None
        self.outfile_path = outfile_path
        self.infile_path = infile_path
        self.year_id = year_ids
        self.release_ids = release_ids
        self.key_id = key_ids
        self.n_parts = 10000

        # Fusion and statistics of datasets for plotting:
        self.fuse_dominant_paths()

    def fuse_dominant_paths(self):
        c_i = 0
        year = self.year_id
        key = self.key_id
        for release in self.release_ids:
            c_i = c_i + 1
            file = self.infile_path + key + '_' + str(year) + '_R' + str(release) + '_trajectory_analysis.nc'
            f1 = nc.Dataset(file)
            if c_i == 1:
                self.dom_path_i = f1['dom_paths'][:]
                f1.close()
            else:
                self.dom_path_i = self.dom_path_i + f1['dom_paths'][:]
                f1.close()

        #self.dom_path_i = self.dom_path_i/(c_i*self.n_parts)
        #breakpoint()
        self.dom_path_i = self.dom_path_i / np.nanmax(self.dom_path_i)
        if self.dom_path_i[0,0] > 3*np.mean(np.unique(self.dom_path_i)):
            self.dom_path_i[0, 0] = 0

        self.dom_path_i = self.dom_path_i / np.nanmax(self.dom_path_i)





