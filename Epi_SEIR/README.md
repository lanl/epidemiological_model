Human Disease Model README

### Conda Environment
- `conda create --name human-epi python=3.8.3`
- `conda install --channel conda-forge numpy pyyaml pandas scipy pyarrow matplotlib sphinx`
- `conda activate human-epi`

### Files

#### Python Scripts
- **models_main.py**: main function for model - runs each general stage of model.
- **vbdm.py**: main class for model - defines general vector borne disease model.
- **dengue.py**: defines behavior for dengue specific model.
- **wnv.py**: defines behavior for WNV specific model. IMPORTANT: currently implements same function as dengue.py.
- **utils.py**: contains functions for creating loggers, parsing arguments, and function timing.
- **generate_inputs.py**: for testing purposes - generates dummy input for mosquito and human population inputs.
- **plotting.py**: CURRENTLY UNDEVELOPED - for plotting the output of model run.

#### Directories
- **config**: contains main configuration file config.yaml
- **docs**: sphinx-autodoc documentation.
- **human_model_output**: contains output of the model.
- **initial_states_input**: contains yaml file with initial population states.
- **logs**: contains logfiles for each model run.
- **mosquitoes**: contains mosquito population vector, 1 point for each day.

#### Shell Scripts
- **build_docs.sh**: builds documentation using sphinx-autodoc.
- **clean.sh**: for testing purposes - moves old logfiles to logs/logfile_archive and deletes old model output.
- **run_human_epi_model.sh**: generates dummy mosquito population input and runs the model for dengue and WNV.
