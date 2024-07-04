import numpy as np
import netCDF4 as nc
import csv
import os
import pandas as pd
from netCDF4 import num2date

class StoreReleases:

    def __init__(self, file_explorer, release_number):
        self.f_path = file_explorer
        self.release_number = release_number
        return

    def read_analysis_df(self, fpath, analysis_file, counter):
        self.analysis_file = fpath.analysis_path + analysis_file
        self.analysis_df = nc.Dataset(self.analysis_file, 'r')
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
        self.shp_p = self.analysis_df.n_part
        self.shp_t = np.shape(self.analysis_vardict['chl_exp'])[0]
        self.sim_start_day = self.analysis_df.sim_start_day
        self.sim_start_month = self.analysis_df.sim_start_month
        self.sim_start_year = self.analysis_df.sim_start_year
        self.duration = self.analysis_df.sim_duration_days
        self.release_number = self.analysis_df.release_number
        self.key = self.analysis_df.scenario_key
        self.sim_start_year = self.analysis_df.sim_start_year
        self.lon_bin_vals = self.analysis_vardict['lon_bin_vals'][:]
        self.lat_bin_vals = self.analysis_vardict['lat_bin_vals'][:]
        self.save_prefix = (self.key + '_' + str(self.sim_start_year) + '_')
        self.save_path = self.f_path.compile_path + self.save_prefix

        if counter == 0:
            self.init_variables()
        self.lon_init = self.analysis_vardict['lon_init'][:]
        self.lat_init = self.analysis_vardict['lat_init'][:]
        self.time = num2date(self.analysis_vardict['time'], self.analysis_vardict['time'].units)
        return

    def init_variables(self):
        self.shp_r = np.shape(self.f_path.file_list)[0]
        # initialise time series variables:
        self.CG_lon = np.zeros([self.shp_t, self.shp_r])
        self.CG_lat = np.zeros([self.shp_t, self.shp_r])
        self.o2_exp = np.zeros([self.shp_t, self.shp_r])
        self.chl_exp = np.zeros([self.shp_t, self.shp_r])
        self.temp_exp = np.zeros([self.shp_t, self.shp_r])

        # initialise for 1D data
        self.site_recruits = np.zeros([self.shp_p])
        self.dom_paths = np.zeros(np.shape(self.analysis_vardict['dom_paths']))
        return

    def update_environment(self, counter):
        self.CG_lon[:, counter] = self.analysis_vardict['CG_lon'][:]
        self.CG_lat[:, counter] = self.analysis_vardict['CG_lat'][:]
        self.o2_exp[:, counter] = self.analysis_vardict['o2_exp'][:]
        self.chl_exp[:, counter] = self.analysis_vardict['chl_exp'][:]
        self.temp_exp[:, counter] = self.analysis_vardict['temp_exp'][:]
        if counter + 1 == self.shp_r:
            self.save_environment()
        return

    def save_environment(self):
        np.save(self.save_path + 'CG_lon.npy', np.array(self.CG_lon))
        np.save(self.save_path + 'CG_lat.npy', np.array(self.CG_lat))
        np.save(self.save_path + 'o2_exp.npy', np.array(self.o2_exp))
        np.save(self.save_path + 'chl_exp.npy', np.array(self.chl_exp))
        np.save(self.save_path + 'temp_exp.npy', np.array(self.temp_exp))
        return

    def update_dom_paths(self, counter):
        self.dom_paths = self.dom_paths + self.analysis_vardict['dom_paths'][:]
        if counter + 1 == self.shp_r:
            np.save(self.save_path + 'dom_paths.npy', np.array(self.dom_paths))
        return

    def update_recruits(self, counter):
        if counter == 0:
            headers = ['release_id', 'date', 'recruit_number', 'recruit_time']
            self.recruit_filename = 'recruit_SG.csv'
            self.init_file(headers, filename=self.recruit_filename)
        release_id = counter + 0
        recruit_t = self.analysis_vardict['recruit_SG_north'][:, 0]
        recruit_index = self.analysis_vardict['recruit_SG_north'][:, 1]
        site_vals = recruit_t > 0
        self.site_recruits = self.site_recruits + site_vals * 1
        id1 = np.where(recruit_t > 0)
        recruit_number = np.shape(id1)[1]
        recruit_time = np.mean(recruit_t[recruit_t > 0])
        day_r = self.sim_start_day
        month_r = self.sim_start_month
        year_r = self.sim_start_year
        dates = ("{:02d}".format(day_r.astype(int)) + '/' + "{:02d}".format(month_r.astype(int)) + '/' +
                 str(year_r.astype(int)))
        df_row = [release_id, dates, recruit_number, recruit_time]
        self.write_file(filename=self.recruit_filename, row=df_row)
        if counter + 1 == self.shp_r:
            np.save(self.save_path + 'site_recruits.npy',
                    np.array([self.lon_init, self.lat_init, self.site_recruits])) #todo: add lat_init, lon_init etc.
        return

    #todo: save table with lat_init, lon_init, recruit numbers, recruit_time etc.
    def init_file(self, headers, filename):
        csv_filename = self.f_path.compile_path + self.save_prefix + filename
        with open(csv_filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(headers)
        return

    def write_file(self, filename, row):
        csv_filename = self.f_path.compile_path + self.save_prefix + filename
        with open(csv_filename, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(row)
        return

    # def read_csv_file(self, filename):
    #     csv_filename = self.f_path.compile_path + self.save_prefix + filename
    #     csv_data = pd.read_csv(csv_filename)
    #     breakpoint()









