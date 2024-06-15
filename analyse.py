# Firstly, specify location of netcdf input files
from configure import FileExplorer
from post_process import PostProcess

# Ready the file paths for analysis
fpath = FileExplorer(node='local', model_name='cmems', key="APSO")
fpath.search_path(year=2017)  # select files which should be analysed

# Post-process simulation file, saving intermediate data (unique particle visits, transit times ...)
#todo: make an instance of the process object that accepts kwargs- keyword list carrying a key value;
for file_i in fpath.file_list:
    print('Analysing file: ' + file_i)
    pp = PostProcess(fpath, file_i, test=False)
    pp.trajectory_analysis()






