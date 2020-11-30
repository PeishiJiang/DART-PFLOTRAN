"""Read in the file and directory paths from file_paths.nml"""

# Author: Peishi Jiang

from os.path import dirname, abspath, join
import f90nml


def read_filepaths_nml(app_dir, dart_pf_dir):
    """
    Inputs:
        app_dir     -- the path of the application folder
        dart_pf_dir -- the path of the dart-pflotran folder
    Outputs:
        dirs_cfg   -- a dictionary of locations for directories
        files_cfg  -- a dictionary of locations for files
    """

    # Get the location of file_paths.nml
    file_paths_nml = join(dart_pf_dir, "file_paths.nml")

    # Read in the file_paths.nml
    nml = f90nml.read(file_paths_nml)

    # Get the directory paths
    dirs_cfg = {}
    for key, value in nml["app_dir_nml"].items():
        dirs_cfg[key] = value.replace("[app_dir]", app_dir)
    for key, value in nml["dart_dir_nml"].items():
        dirs_cfg[key] = value.replace("[dart_pf_dir]", dart_pf_dir)

    # Get the file paths
    files_cfg = {}
    for nml_item in nml.keys():
        if nml_item in ["app_dir_nml", "dart_dir_nml"]:
            continue
        for key, value in nml[nml_item].items():
            if key in ['pflotran_out_prefix','use_obs_tecfile_for_prior']:
                continue
            elif "[obs_type_dir]" in value:
                files_cfg[key] = value.replace("[obs_type_dir]",
                                               dirs_cfg["obs_type_dir"])
            elif "[app_work_dir]" in value:
                files_cfg[key] = value.replace("[app_work_dir]",
                                               dirs_cfg["app_work_dir"])
            elif "[pflotran_in_dir]" in value:
                files_cfg[key] = value.replace("[pflotran_in_dir]",
                                               dirs_cfg["pflotran_in_dir"])
            elif "[pflotran_out_dir]" in value:
                files_cfg[key] = value.replace("[pflotran_out_dir]",
                                               dirs_cfg["pflotran_out_dir"])
            elif "[dart_data_dir]" in value:
                files_cfg[key] = value.replace("[dart_data_dir]",
                                               dirs_cfg["dart_data_dir"])
            elif "[obs_kind_dir]" in value:
                files_cfg[key] = value.replace("[obs_kind_dir]",
                                               dirs_cfg["obs_kind_dir"])
            elif "[utils_dir]" in value:
                files_cfg[key] = value.replace("[utils_dir]",
                                               dirs_cfg["utils_dir"])
            elif "[dart_work_dir]" in value:
                files_cfg[key] = value.replace("[dart_work_dir]",
                                               dirs_cfg["dart_work_dir"])
            else:
                raise Exception("Unknown value: {}".format(value))

    return dirs_cfg, files_cfg
