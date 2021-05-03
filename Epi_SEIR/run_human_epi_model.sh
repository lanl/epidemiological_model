#!/bin/sh

./clean.sh

# generate mosquito input for development of model
python generate_inputs.py

python models_main.py -c /Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml -d dengue
#python models_main.py -c /Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml -d wnv

echo " "
cat logs/*.log
