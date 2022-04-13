#!/bin/sh

# Takes two arguments:
# 1. Config file (e.g., ./config/human_epi_config.yaml)
# 2. Directory for the mosquito population
# 3. A flag that is set to 1 if one needs to run clean.sh and 0 otherwise

if [ "$1" == "" ]; then
    echo "Error: Must specify configuration file path"
    exit 1
fi

if [ "$2" == "" ]; then
    echo "Error: Must specify mosquito data directory"
    exit 1
fi

CONFIG_FILE_PATH=$1
MOSQ_DIR=$2
NEED_TO_CLEANUP=$3

if [ "$NEED_TO_CLEANUP" -eq 1 ]; then
    ./clean.sh $CONFIG_FILE_PATH
fi


#python models_main.py -c $CONFIG_FILE_PATH -m $MOSQ_DIR -d dengue
python models_LLM_human.py -c $CONFIG_FILE_PATH -m $MOSQ_DIR -d wnv
