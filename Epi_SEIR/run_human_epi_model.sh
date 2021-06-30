#!/bin/sh

# Takes two arguments:
# 1. Config file (e.g., ./config/human_epi_config.yaml)
# 2. A flag that is set to 1 if one needs to run clean.sh and 0 otherwise

if [ "$1" == "" ]; then
    echo "Error: Must specify configuration file path"
    exit 1
fi

CONFIG_FILE_PATH=$1
NEED_TO_CLEANUP=$2

if [ "$NEED_TO_CLEANUP" -eq 1 ]; then
    ./clean.sh $CONFIG_FILE_PATH
fi

# generate mosquito input for development of model
# python generate_inputs.py -c $CONFIG_FILE_PATH 

python models_main.py -c $CONFIG_FILE_PATH -d dengue
python models_main.py -c $CONFIG_FILE_PATH -d wnv
