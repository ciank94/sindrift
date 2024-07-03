import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.patches import Polygon
from shapely.geometry.polygon import LinearRing, Polygon
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import netCDF4 as nc
import os
from netCDF4 import num2date


class ReadAnalysis:

    def __init__(self, fpath, analysis_file):
        self.analysis_file = fpath.analysis_path + analysis_file
        self.analysis_df = nc.Dataset(self.analysis_file, 'r')
        self.tr_file = fpath.trajectory_path + self.analysis_df.trajectory_file
        self.fpath = fpath
        return

    def read_trajectory_df(self):
        self.trajectory_df = nc.Dataset(self.tr_file, 'r')
        # extract data from trajectory file
        self.lat = self.trajectory_df.variables['lat']
        self.lon = self.trajectory_df.variables['lon']
        self.shp_p = np.shape(self.lat)[0]
        self.shp_t = np.shape(self.lat)[1]
        time = self.trajectory_df.variables['time']
        self.time = num2date(time, time.units)
        return

    def read_analysis_df(self):
        # subset variables from file in dictionary
        self.analysis_vardict = dict()
        for key in self.analysis_df.variables:
            self.analysis_vardict[key] = self.analysis_df.variables[key]

        # subset dimension size from file
        self.n_lat_bins = self.analysis_df.dimensions['n_lat_bins'].size
        self.n_lon_bins = self.analysis_df.dimensions['n_lon_bins'].size
        self.shp_t = self.analysis_df.dimensions['obs'].size

        # subset attributes from file
        self.bin_res = self.analysis_df.bin_resolution
        self.key = self.analysis_df.scenario_key
        self.min_lon = self.analysis_df.domain_lon_min
        self.min_lat = self.analysis_df.domain_lat_min
        self.max_lon = self.analysis_df.domain_lon_max
        self.max_lat = self.analysis_df.domain_lat_max
        self.p_shp = self.analysis_df.n_part
        self.sim_start_day = self.analysis_df.sim_start_day
        self.sim_start_month = self.analysis_df.sim_start_month
        self.sim_start_year = self.analysis_df.sim_start_year
        self.duration = self.analysis_df.sim_duration_days
        self.release_number = self.analysis_df.release_number
        return

class StoreRelease:
    def __init__(self, fpath):
        self.shp_r = np.shape(fpath.file_list)[0]
        self.counter_r = -1
        self.fpath = fpath

    def init_variables(self, rd):
        # initialise time series variables:
        self.CG_lon = np.zeros([rd.shp_t, self.shp_r])
        self.CG_lat = np.zeros([rd.shp_t, self.shp_r])
        self.o2_exp = np.zeros([rd.shp_t, self.shp_r])
        self.chl_exp = np.zeros([rd.shp_t, self.shp_r])
        self.temp_exp = np.zeros([rd.shp_t, self.shp_r])

        # initialise for 1D data
        self.recruit_number = np.zeros([self.shp_r])
        self.recruit_time = np.zeros([self.shp_r])
        self.year_r = np.zeros([self.shp_r])
        self.month_r = np.zeros([self.shp_r])
        self.day_r = np.zeros([self.shp_r])
        self.release_number_r = np.zeros([self.shp_r])
        self.p_shp = rd.p_shp
        self.site_recruits = np.zeros([rd.p_shp])
        self.dom_paths = np.zeros(np.shape(rd.dom_paths))
        self.rd = rd
        self.key = rd.scenario_key
        self.sim_start_year = rd.sim_start_year
        self.lon_bin_vals = rd.analysis_vardict['lon_bin_vals'][:]
        self.lat_bin_vals = rd.analysis_vardict['lat_bin_vals'][:]
        return

    def store_variables(self, rd):
        self.CG_lon[:, self.counter_r] = rd.analysis_vardict['CG_lon'][:]
        self.CG_lat[:, self.counter_r] = rd.analysis_vardict['CG_lat'][:]
        self.o2_exp[:, self.counter_r] = rd.analysis_vardict['o2_exp'][:]
        self.chl_exp[:, self.counter_r] = rd.analysis_vardict['chl_exp'][:]
        self.temp_exp[:, self.counter_r] = rd.analysis_vardict['temp_exp'][:]
        self.count_recruits(rd)
        self.dom_paths = self.dom_paths + rd.dom_paths
        return

    def store_trajectory_variables(self, t_df):
        self.lon_init = np.array(t_df.lon[:, 0])
        self.lat_init = np.array(t_df.lat[:, 0])
        return

    def count_recruits(self, rd):
        recruit = rd.recruit[:, :]
        recruit_t = recruit[:, 0]
        id1 = np.where(recruit_t > 0)
        site_vals = recruit_t>0
        self.site_recruits = self.site_recruits + site_vals*1
        self.recruit_number[self.counter_r] = np.shape(id1)[1]
        self.recruit_time[self.counter_r] = np.mean(recruit_t[recruit_t>0])
        self.year_r[self.counter_r] = self.sim_start_year
        self.month_r[self.counter_r] = rd.sim_start_month
        self.day_r[self.counter_r] = rd.sim_start_day
        self.release_number_r[self.counter_r] = rd.release_number
        return

    def write_data(self):
        import pandas as pd
        from datetime import datetime
        self.save_prefix = (self.key + '_' + str(self.sim_start_year) + '_')
        day_r = self.day_r
        month_r = self.month_r
        year_r = self.year_r
        shp_r = np.shape(self.day_r)[0]
        dates = ["{:02d}".format(day_r[i].astype(int)) + '/' + "{:02d}".format(month_r[i].astype(int)) #+ '/' + str(year_r[i].astype(int))
                    for i in range(0, shp_r)]
        df_table = pd.DataFrame([self.release_number_r, dates, self.recruit_number.astype(int), self.recruit_time]).T
        df_table = df_table.sort_values(by=[0])
        table_name = self.fpath.figures_path + self.save_prefix + 'recruit_table.csv'
        df_table.to_csv(table_name, header=['release_number', 'date', 'recruit_number', 'recruit_time'])
        recruit_table = pd.read_csv(table_name)
        self.recruit_table = recruit_table.iloc[:,1:np.shape(recruit_table)[1]] # adds extra column for some reason
        return

    def counter(self):
        self.counter_r = self.counter_r + 1
        print("Analysing release number: " + str(self.counter_r + 1))
        return


class PlotData:

    def __init__(self, df):
        self.fpath = df.fpath
        self.save_prefix = (df.key + '_' + str(df.sim_start_year) + '_')

        # polygon for analysis
        self.n_poly = 0
        self.poly = self.get_polygon(self.n_poly)

        # For creating the bathymetry file
        self.load_bathymetry()

        # parameters for plotting:
        self.load_plot_params()
        #self.lon_lat_extent()
        return


    def load_plot_params(self):

        # Dominant pathways color scaling:
        self.d_scale = 40
        self.max_scale = 0.4

        # plotting parameters
        self.bath_contours = np.linspace(0, 3000, 10)
        self.bath_cmap = plt.get_cmap('Blues')
        #self.dom_cmap = plt.get_cmap('Reds')
        self.dom_cmap = plt.get_cmap('OrRd')
        self.depth_colors = np.arange(0, 4500, 200)
        self.site_recruit_cmap = plt.get_cmap('jet')

        # offsets for figure boundaries:
        self.lon_offset = 1
        self.lat_offset = 3
        return



    def plot_CG_paths(self, df):
        self.init_plot()
        self.plot_background()
        [self.ax.plot(df.CG[:, 0, i], df.CG[:, 1, i], alpha=0.9, linewidth=0.5, c='k') for i in range(0, np.shape(df.CG)[2])]
        plt_name = self.save_prefix + "CG_plot"
        self.save_plot(plt_name)
        return

    def plot_site_recruits(self, df):
        self.init_plot()
        self.plot_background(background="SOI")
        c_vals = df.site_recruits[:]/(df.counter_r + 1) * 100
        self.ax.scatter(df.lon_init, df.lat_init, s=10, facecolor='none', edgecolors='gray', alpha=0.3, linewidth=0.2)
        self.plot1 = plt.scatter(df.lon_init, df.lat_init, c=c_vals,s=10, edgecolors='gray', vmin=np.nanmean(c_vals)/2,
                               vmax=np.nanmean(c_vals)*2, linewidth=0.2, cmap=self.site_recruit_cmap)
        self.add_cbar(c_max=np.nanmax(c_vals), caxis_title='Percentage (%)')
        self.save_plot(plt_name=self.save_prefix + 'site_recruits')
        return






    def plot_init(self):
        self.init_plot()
        self.plot_background()
        plt.scatter(self.lon[:,0], self.lat[:,0])
        plt_name = self.save_prefix + "init_plot"
        self.save_plot(plt_name)
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
        # todo: use digitize instead of argmin;
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


    def plot_dom_paths(self, df):
        self.init_plot()
        self.plot_background(background='n')
        dom_paths = df.dom_paths.astype(float)
        dom_paths[dom_paths==0] = np.nan
        dom_paths = (dom_paths/((df.counter_r +1)*10000)) * 100
        max_vals = np.nanmax(dom_paths)*self.max_scale
        n_levels = np.arange(np.nanmin(dom_paths), max_vals, max_vals/25)
        self.plot1 = plt.contourf(df.lon_bin_vals, df.lat_bin_vals, dom_paths.T, levels=n_levels, cmap=self.dom_cmap, transform=ccrs.PlateCarree(), extend='both')
        self.c_max = max_vals
        self.add_cbar(c_max= self.c_max, caxis_title='unique_particles')
        plt_name = self.save_prefix + "dom_paths"
        self.save_plot(plt_name)
        return

    def plot_time_recruits(self, df):
        dates = df.recruit_table.date
        recruit_n = df.recruit_table.recruit_number
        recruit_t = df.recruit_table.recruit_time/24

        fig, ax1 = plt.subplots()
        plt.xticks(rotation=45, fontsize=10)
        ax2 = ax1.twinx()
        ax1.set_ylabel('recruit time (days)', color='b', fontsize=11)
        ax2.set_ylabel('recruit number', color='r', fontsize=11)
        ax1.plot(dates, recruit_t, c = 'b')
        ax2.plot(dates, recruit_n, c = 'r')
        #ax1.set_xlabel('date', fontsize=8)
        plt.grid(alpha=0.45)  # nice and clean grid
        plt_name = self.save_prefix + "recruit_times"
        self.save_plot(plt_name)



    def plot_recruits(self):
        recruit = self.recruit[:,:]
        recruit_t = recruit[:,0]
        recruit_i = recruit[:,1]
        id1 = np.where(recruit_t>0)
        self.init_plot()
        self.plot_background()
        i_vals = recruit_i[id1[0][:]].astype(int)
        particle_n = id1[0][:].astype(int)
        x = np.array(self.lon[particle_n, :])
        y = np.array(self.lat[particle_n, :])
        [self.ax.plot(x[i, 0:i_vals[i]], y[i, 0:i_vals[i]], alpha=0.1, linewidth=0.5, c='r') for i in range(0, np.shape(id1)[1])]
        self.ax.scatter(x[:, 0], y[:, 0], alpha=0.7, c='r', s=1)
        x,y = self.poly.exterior.xy
        self.ax.plot(x, y, color='y', alpha=0.9,
                linewidth=2, solid_capstyle='round', zorder=2)
        plt_name = self.save_prefix + "recruit"
        self.save_plot(plt_name)

        n_recruits = np.shape(id1)[1]


        #self.ax.plot(polygon1)
        #self.ax.add_patch(polygon1, facecolor='r', alpha=0.4)
        #
        #breakpoint()
        return n_recruits



    def get_polygon(self, poly_n):
        lon_lims = [-39.5, -35]
        lat_lims = [-54, -53]
        self.lon_1 = lon_lims[0]
        self.lon_2 = lon_lims[1]
        self.lat_1 = lat_lims[0]
        self.lat_2 = lat_lims[1]
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
        plt.scatter(lon_1[:, -1], lat_1[:, -1], s=18, facecolor='red', edgecolors='k', alpha=0.9, linewidth=0.6)
        plt.scatter(lon_1[:, 0], lat_1[:, 0], s=5, facecolor='yellow', edgecolors='g', alpha=0.3, linewidth=0.5)

        plt_name = self.save_prefix + "worms"
        self.save_plot(plt_name)
        return

    def plot_background(self, background):
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
        #self.ax.set_extent(
            #[self.min_lon, self.max_lon, self.min_lat, self.max_lat])
        if background == "AP":
            self.AP_lon_lat_extent()
            self.ax.set_extent(
                [self.min_lon, self.max_lon, self.min_lat,
                 self.max_lat])
            #self.ax.set_extent(
             #   [self.min_lon - self.lon_offset, self.max_lon + self.lon_offset, self.min_lat - self.lat_offset,
              #   self.max_lat + self.lat_offset])
        elif background == "SOI":
            self.SOI_lon_lat_extent()
            self.ax.set_extent(
                [self.min_lon, self.max_lon, self.min_lat,
                 self.max_lat])
        else:
            self.gen_lon_lat_extent()
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
        cbar = plt.colorbar(self.plot1, extend='both', pad=0.01)
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

    def AP_lon_lat_extent(self):
        self.min_lon = -65.3
        self.max_lon = -51
        self.min_lat = -69
        self.max_lat = -56
        return

    def SOI_lon_lat_extent(self):
        #todo: put all site_lims in same function and add ifelse statements
        self.min_lon = -49
        self.max_lon = -41
        self.min_lat = -64.5
        self.max_lat = -58
        return

    def gen_lon_lat_extent(self):
        # SG extent;
        self.min_lon = -45
        self.max_lon = -33
        self.min_lat = -57
        self.max_lat = -52
        return






