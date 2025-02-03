import os, socket

# Un-comment to render in Markdown
# Universal experiment name configuration
# /`root`/`model_name`/`root_work_designator`/`exp_name`/`output_designator`

# | Component field name   | Component long name | Component description |
# |:---------------------- |:------------------- |:--------------------- |
# | `root_dir`             | Root directory      | Name of root directory, usually `project` or `scratch` followed by a pointer to user-dependent subdirectory based on Princeton netIDs |
# | `model_name`           | Model name          | name of the model being run (e.g., AM2.5, HIRAM, FLOR) |
# | `root_work_designator` | Experiment parent directory name | name of the parent directory for all experiments (e.g., `work`, `exp`) |
# | `exp_name`             | Experiment run name | name of the experiment run |
# | `output_designator`    | Output directory name | name of the directory containing output netCDF files |

def get_dir(root, model_name, experiment, experiment_base='CTL1990s', group='GEOCLIM', username='gr7610'):
    """
    Method to obtain absolute path to directory of interest for a given model and experiment.

    Args:
        root (str): name of base filesystem (typically 'projects' or 'scratch')
        model_name (str): name of model of interest (typically 'AM2.5', 'HIRAM', or 'FLOR')
        experiment_base (str, optional): basename of experiments of interest. Defaults to 'CTL1990s'.
        group (str, optional): name of group under which files are stored. Defaults to 'GEOCLIM'.
        username (str, optional): username that file access is desired for. Defaults to 'gr7610'.

    Returns:
        path (str): absolute path to directory of interest for a given model and experiment.
    """

    # Get hostname. Useful for accessing files on Jupyter when worknig on 'tigressdata'.
    hostname = socket.gethostname()
    # Pre-define platform name
    platform = 'tigercpu_intelmpi_18'
    # Pre-define number of processors
    if model_name in ['FLOR']:
        npes = 576
    elif model_name in ['AM2.5C360']:
        npes = 1080
    else:
        npes = 540
    # Get root name. Append /tiger to directory path if working on Tigress.
    if root == 'projects':
        root_dir = '/projects/GEOCLIM/{0}'.format(username) if group == 'GEOCLIM' else '/projects/{0}'.format(username)
    elif root == 'scratch':
        if 'tigress' in hostname:
            root_dir = '/tiger/scratch/gpfs/GEOCLIM/{0}'.format(username) if group == 'GEOCLIM' else '/scratch/gpfs/{0}'.format(username)
        else:
            root_dir = '/scratch/gpfs/GEOCLIM/{0}'.format(username) if group == 'GEOCLIM' else '/scratch/gpfs/{0}'.format(username)    
    # Define auxiliary designators.
    root_work_designator = 'exp' if root == 'projects' else 'work'
    output_designator = 'work/POSTP' if root == 'projects' else 'POSTP'
    # Re-define path to aid search for experiments
    path = os.path.join(root_dir, model_name, root_work_designator)
    # Initialize experiment name.
    exp_name = None
    for subdir in os.listdir(path):
        # Assume control experiments have no suffixes to the base name (experiment_base)
        # Else, look for everything after the first underscore
        if experiment_base in subdir:
            # Define suffix for platform metadata
            platform_suffix = '_{0}_{1}PE'.format(platform, npes) if root == 'scratch' else ''
            # Get experiment name from full string
            subdir_ = subdir.replace(experiment_base + '_', '').replace(platform_suffix, '')
            # Assign experiment names
            if experiment == 'control' and subdir_ == experiment_base:
                # Re-define experiment name if on /scratch
                exp_name = subdir_
            elif experiment is not 'control' and subdir_ == experiment:
                exp_name = experiment_base + '_' + subdir_ + platform_suffix
    # Use conditional to direct output flow. If exp_name is defined, use it. Else, return None.
    if exp_name is None:
        print('Experiment {0} not found in {1}...'.format(experiment, path))
        return None
    else:
        # Reconstruct path with experiment name and output deisgnator.
        # Use 'realpath' if a symlink is picked up on /projects.
        path = os.path.realpath(os.path.join(path, exp_name, output_designator))
        # Return path if hosting it on 
        if root == 'projects' and 'scratch' in path and 'tigress' in hostname:\
            path = '/tiger' + path

        return path