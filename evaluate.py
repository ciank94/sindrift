#todo: add functionality for plotting diagnostics from simulation- add options so that the worms plot can be done alone as a check
from configure import FileExplorer
from plot_data import PlotData

# Ready the file paths for analysis
fpath = FileExplorer(node='local', model_name='sinmod')
fpath.search_path()  # select files which should be plotted

for file_v in fpath.file_list:
    analysis_file = fpath.analysis_path + file_v
    print('Analysing file: ' + analysis_file)
    pld = PlotData(fpath, analysis_file)
    pld.plot_trajectory()
    dom_paths = pld.analysis_df.variables['dom_paths'][:]
    pld.plot_dom_paths(dom_paths)
    breakpoint()