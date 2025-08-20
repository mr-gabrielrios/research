'''
Utility script to create directories for a model-experiment configuration in a long-term filesystem.

Program flow:
1. Define an experiment configuration to create a directory for from `scratch` to a storage directory `projects`
2. Locate the configuration root directory on `projects`
3. Scrape the link to the corresponding work directory on `scratch` and save the link to local memory
4. Create a storage directory for the configuration, if it doesn't already exist (allow an overwrite option)
5. Create symbolic links to the experiment parent directory on `projects`, the TC tracking directory on `scratch`, and the postprocessed climate model output on `scratch`

'''

import argparse
import os
import sys

def generate_symlinks(experiment_root_dirname: str,
                      experiment_work_dirname: str,
                      storage_root_dirname: str,
                      overwrite: bool,
                      diagnostic: bool):
    ''' 
    Create symbolic links to the experiment parent directory on projects, 
    the TC tracking directory on scratch, 
    and the postprocessed climate model output on scratch. 
    '''
    
    experiment_subdir_names = ['analysis_lmh', 'exp', 'POSTP']
    
    for experiment_subdir_name in experiment_subdir_names:
            
        experiment_subdir_root = experiment_root_dirname if experiment_subdir_name == 'exp' else experiment_work_dirname
        print(f'Source: base directory {experiment_subdir_root} and subdirectory specified {experiment_subdir_name}...')
        print(f'Destination: base directory {storage_root_dirname} and subdirectory specified {experiment_subdir_name}...')
    
        experiment_subdir_srcdir = experiment_subdir_root if experiment_subdir_name == 'exp' else os.path.join(experiment_subdir_root, experiment_subdir_name)
        experiment_subdir_dstdir = os.path.join(storage_root_dirname, experiment_subdir_name)
        
        if diagnostic:
            print(f'Generating subdirectory for {experiment_subdir_name} at {experiment_subdir_srcdir}...')

        if not os.path.islink(experiment_subdir_dstdir) or not os.path.isdir(experiment_subdir_dstdir):
            if overwrite:
                os.remove(experiment_subdir_dstdir)
                
            os.symlink(experiment_subdir_srcdir, experiment_subdir_dstdir)

def main(model_name: str,
         experiment_name: str,
         overwrite: bool=False,
         diagnostic: bool=True):

    # 1. Define an experiment configuration to create a directory for from scratch to a storage directory projects
    temporary_root_dirname = '/scratch/gpfs/GEOCLIM/gr7610/tiger3'
    permanent_root_dirname = '/projects/GEOCLIM/gr7610/MODEL_OUT'
    compiler_handle = 'tiger3_intelmpi_24_540PE'

    # 2. Locate the configuration root directory on projects
    
    supported_model_names = ['HIRAM', 'AM2.5', 'FLOR']
    assert model_name in supported_model_names, f'Model name must be in {supported_model_names}.'
    root_dirname_model_name = 'CM2.5' if model_name in ['AM2.5', 'FLOR'] else 'HIRAM'
    
    experiment_root_dirname = os.path.join(temporary_root_dirname, root_dirname_model_name, 'exp', experiment_name)
    if diagnostic:
        print(f'Experiment root directory defined as {experiment_root_dirname}...')
    
    # Ensure the experiment directory (experiment_root_dirname) exists
    assert os.path.isdir(experiment_root_dirname) is True, f'Experiment root directory {experiment_root_dirname} not found, please try again.'

    # 3. Scrape the link to the corresponding work directory on scratch
    experiment_work_dirname = os.readlink(os.path.join(experiment_root_dirname, 'work'))
    if diagnostic:
        print(f'Experiment working directory defined as {experiment_work_dirname}...')
    
    # Ensure the experiment working directory (experiment_work_dirname) exists
    assert os.path.isdir(experiment_work_dirname) is True, f'Experiment working directory {experiment_work_dirname} not found, please try again.'

    # 4. Create a storage directory for the configuration, if it doesn't already exist (allow an overwrite option)
    storage_root_dirname = os.path.join(permanent_root_dirname, model_name, experiment_name)
    
    if (not os.path.isdir(storage_root_dirname) or len(os.listdir(storage_root_dirname)) == 0) or overwrite:
        if diagnostic:
            print(f'Creating storage root directory {storage_root_dirname}...')

        if not os.path.isdir(storage_root_dirname):
            os.mkdir(storage_root_dirname)
        
        generate_symlinks(experiment_root_dirname=experiment_root_dirname,
                          experiment_work_dirname=experiment_work_dirname,
                          storage_root_dirname=storage_root_dirname,
                          overwrite=overwrite,
                          diagnostic=diagnostic)
    else:
        print(f'Storage root directory {storage_root_dirname} already exists, not overwriting...')
        sys.exit()

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--model_name", type=str, help="Name of climate model to create directories for.")
    parser.add_argument("--experiment_name", type=str, help="Name of experiment to create directories for.")
    parser.add_argument("--overwrite", type=bool, default=False, help="Boolean to control whether any existing paths are overwritten.")
    args = parser.parse_args()
    
    main(model_name=args.model_name, experiment_name=args.experiment_name, overwrite=args.overwrite)