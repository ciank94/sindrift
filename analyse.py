# Firstly, specify location of netcdf input files
import numpy as np
from configure import FileExplorer
from post_process import PostProcess
from compile_releases import StoreReleases

# Ready the file paths for analysis, specifying whether running code and accessing files locally or remotely
node_name = 'local'
file_loc = 'local'
fpath = FileExplorer(node=node_name, model_name='cmems', key="SOIN")
if node_name == 'local':
    if file_loc == 'remote':
        fpath.mounted_paths()
    else:
        fpath.local_phys_states()
fpath.search_path(year=2016, release_start=1, release_end=1)  # select files which should be analysed

# Post-process simulation file, saving intermediate data (unique particle visits, transit times ...)
#todo: make an instance of the process object that accepts kwargs- keyword list carrying a key value;
for file_i in fpath.file_list:
    print('Analysing file: ' + file_i)
    pp = PostProcess(fpath, file_i, test=False)
    pp.init_ncfile()
    pp.trajectory_analysis()
    pp.trajectory_df.close()
    pp.analysis_df.close()
    print('Closing: ' + pp.analysis_file)

# first initialise files for storage
n_releases = np.shape(fpath.file_list)[0]
st = StoreReleases(file_explorer=fpath, release_number=n_releases)
counter = -1
for file_i in fpath.file_list:
    counter = counter + 1
    st.read_analysis_df(fpath, analysis_file=file_i, counter=counter)
    st.update_recruits(counter)
    st.update_environment(counter)
    st.update_dom_paths(counter)
    st.analysis_df.close()






