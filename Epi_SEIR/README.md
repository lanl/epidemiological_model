### *Standard Workflow*
- Configure conda environment (below)
- Pull changes
	- `git pull`
- Create new git branch
	- `git branch <branch_name>`
	- `git checkout <branch_name>`
- Develop changes
- Commit changes `git commit -m <commit_message>`
- Verify continuous integration job passed
	- left sidebar -> CI/CD -> Pipelines
- Merge branch with master
	- `git checkout master`
	- `git merge <branch_name>`

### *Conda Environment Creation*
To create environment from text file specification, run
- `conda create --name human-epi-env --file conda_environment/human-epi-env.txt`
 
To create environment manually, run the following  
- `conda create --name human-epi-env python=3.8.3`
- `conda activate human-epi-env`
- `conda install --channel conda-forge numpy pyyaml pandas scipy pyarrow matplotlib sphinx pytest`

### *Running Model*
- `./run_human_epi_model.sh <config_file_path> <[1 to run clean script]|[0 to not run clean script]>`

### *Building Documentation*
- `pip install sphinx_rtd_theme`
- `./build_docs.sh`

### *Unit Testing*
UNDER CONSTRUCTION

### *Files*

#### Python Scripts
- **models_main.py**: main function for model - runs each general stage of model.
- **vbdm.py**: main class for model - defines general vector borne disease model.
- **dengue.py**: defines behavior for dengue specific model.
- **wnv.py**: defines behavior for WNV specific model.
- **utils.py**: contains functions for creating loggers, parsing arguments, and function timing.
- **generate_inputs.py**: for testing purposes - generates dummy input for mosquito and human population inputs.
- **plotting.py**: CURRENTLY UNDEVELOPED - for plotting the output of model run.
- **test_human_epi_model.py**: contains code for unit testing the model.

#### Directories
- **config**: contains main configuration files and unit testing configuration files.
- **docs**: sphinx-autodoc documentation.
- **human_model_output**: contains output of the model.
- **initial_states_input**: contains yaml file with initial population states.
- **logs**: contains logfiles for each model run.
- **mosquitoes**: contains mosquito population vector, 1 point for each day.
- **unit_testing**: contains blueprint function for unit testing and documentation on unit testing error checks.

#### Shell Scripts
- **build_docs.sh**: builds documentation using sphinx-autodoc.
- **clean.sh**: for testing purposes - moves old logfiles to logs/logfile_archive and deletes old model output.
	- usage: `./clean.sh <config_file_path>`
- **run_human_epi_model.sh**: generates dummy mosquito population input and runs the model for dengue and WNV.
	- usage: `./run_human_epi_model.sh <config_file_path> <[1 to run clean script]|[0 to not run clean script]>`
