Human Disease Model README

Directories
- config: contains main configuration file config.yaml
- docs: sphinx-autodoc documentation
- human_model_output: contains output of the model
- initial_states_input: contains yaml file with initial population states
- logs: contains logfiles for each model run
- mosquitoes: contains mosquito population vector, 1 point for each day

Shell Scripts
- build_docs.sh: builds documentation using sphinx-autodoc.
- clean.sh: for model testing purposes - moves old logfiles to logs/logfile_archive and deletes old model output
- run_human_epi_model.sh: generates dummy mosquito population input and runs the model for dengue and WNV
