### *Standard Workflow*
- Configure conda environment (below)
- Pull changes
	- `git pull`
- Create new git branch
	- `git branch <branch_name>`
	- `git checkout <branch_name>`
- Develop changes
- Add changes `git add <specific file> or -A for all files`
- Commit changes `git commit -m <commit_message>`
- Verify continuous integration job passed
	- left sidebar -> CI/CD -> Pipelines
- Merge branch with master by creating Merge Request at gitlab.lanl.gov
- Can merge in the command line as well, but recommend creating a Merge Request
	- `git checkout master`
	- `git merge <branch_name>`

### *Conda Environment Creation*
To create environment from text file specification, run
- `conda env create -f conda_environment/environment.yml`
 
To create environment manually, run the following  
- `conda create --name human-epi-env python=3.9`
- `conda activate human-epi-env`
- `conda install --channel conda-forge numpy pyyaml pandas scipy pyarrow matplotlib sphinx pytest lmfit`

### *Running Model from Shell Scripts*
- `./run_human_epi_model.sh <config_file_path> <[1 to run clean script]|[0 to not run clean script]>`

### *Running Model from Python Scripts*
- `python models_main.py -c <config_file_path> -d <dengue or wnv> -f` (to print model output figures) `-sf` (to save model output figures)
- `python models_params.py -c <config_file_path> -d <denuge or wnv> -l` (to label model outut files by parameter values) `-g` (to generate parameter data file of parameters to sequence through from config/local_param_config.yaml) `-p` <parameter_values_data_file_path> (for if `-g` is not used) `-sf` (to save model output figures)
- `python models_fit.py -c <config_file_path> -d <dengue or wnv> -rm` (to run model with fit parameters) `-f` (to print model output and fitting data figure) `-sf` (to save model output and fitting data figure)

### *Building Documentation*
- `pip install sphinx_rtd_theme`
- `./build_docs.sh`

### *Unit Testing*
- `pytest` is used for unit testing the code.
- run the command `pytest` from the *Epi_SEIR* directory.

### *Gitlab Continuous Integration*
#### SETUP: Building and Pushing the Image
- Install [Docker Desktop](https://www.docker.com/products/docker-desktop).
- Create a GitLab personal access token. [This documentation](https://gitlab.lanl.gov/help/user/profile/personal_access_tokens) walks you through the process. In short, visit https://gitlab.lanl.gov/profile/personal_access_tokens, and create a new personal access token that has `read_registry` and `write_registry` permissions. Copy this token.
- Open a terminal, and change directories to the root of this project.
- Login using your newly-created personal access token: `docker login gitlab.lanl.gov:5050`
    * Your username is your LANL moniker (identical to your GitLab username), and your password is the token.
    * To have Docker remember your credentials locally on a Mac, edit `~/.docker/config.json` to modify the `credsStore` setting to be `osxkeychain` prior to logging in (alternatively, log out using `docker logout gitlab.lanl.gov:5050` and log back in after modifying this setting).
- There is a known [Docker issue](https://github.com/docker/for-mac/issues/2723) where NO_PROXY is not honored. You will need to enable the manual proxy in Docker Desktop when building the image, and disable it before pushing the image. Details are below.
- Go to Docker Desktop GUI -> settings -> resources -> proxies.
- Toggle the "Manual proxy configuration" switch
- Enter "http://proxyout.lanl.gov:8080" for the first two boxes and "\*.lanl.gov" for the third box.
- Build the image: `docker build -t gitlab.lanl.gov:5050/cimmid/disease_and_human_modeling/human_epi_models .`
- Toggle the "Manual proxy configuration" switch again
- Push the image: `docker push gitlab.lanl.gov:5050/cimmid/disease_and_human_modeling/human_epi_models`
- Verify that the image is now present at gitlab.lanl.gov/cimmid/disease_and_human_modeling/human_epi_models/container_registry.
#### Adding a package to the environment
NOTE: There is probably a better way to do this, but this works for now.
- Edit *Dockerfile* in the *human_epi_models* directory.
- Add the name of the package to be added in the `conda install` line of the file.
- Build the image and push it to the container registry as outlined in the **SETUP** steps above.

### *Files*

#### Python Scripts
- **models_main.py**: main function for model - runs each general stage of model.
- **models_params.py**: main function to run model through a series of parameter values.
- **modesl_fit.py**: main function to fit model parameters to data.
- **vbdm.py**: main class for model - defines general vector borne disease model.
- **fit.py**: class for model fitting - defines fitting functions.
- **dengue.py**: defines behavior for dengue specific model.
- **wnv.py**: defines behavior for WNV specific model.
- **utils.py**: contains functions for creating loggers, parsing arguments, and function timing.
- **generate_inputs.py**: for testing purposes - generates dummy input for mosquito and human population inputs.
- **plotting.py**: For plotting the output of model run from an output csv.
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
- **run_human_epi_model.sh**: singe run of the model for dengue and WNV.
	- usage: `./run_human_epi_model.sh <config_file_path> <[1 to run clean script]|[0 to not run clean script]>`
