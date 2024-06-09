# Firstly, specify location of netcdf input files
from configure import FileExplorer
from post_process import PostProcess

# Ready the file paths for analysis
fpath = FileExplorer(node='local', model_name='sinmod')
fpath.search_path()  # select files which should be analysed

# Post-process simulation file, saving intermediate data (unique particle visits, transit times ...)
#todo: make an instance of the process object that accepts kwargs- keyword list carrying a key value;
for file_v in fpath.file_list:
    analysis_file = fpath.analysis_path + file_v
    print('Analysing file: ' + analysis_file)
    pp = PostProcess(analysis_file)
    pp.trajectory_analysis(test=False)






