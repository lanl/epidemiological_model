#!/bin/sh

if [ "$1" == "" ]; then
	echo "Error: Must specify configuration file path"
	exit 1
fi

CONFIG_FILE_PATH=$1

./clean.sh $CONFIG_FILE_PATH

# generate mosquito input for development of model
python generate_inputs.py -c $CONFIG_FILE_PATH 

#python models_main.py -c $CONFIG_FILE_PATH -d dengue
python models_main.py -c $CONFIG_FILE_PATH -d wnv

echo " "
cat logs/*.log
