#!/bin/sh

./clean.sh

# generate mosquito input for development of model
python mosquitoes/generate_mosquito_input.py

python models_main.py -c /Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml -d dengue
python models_main.py -c /Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml -d wnv

echo " "
cat logs/*.log
