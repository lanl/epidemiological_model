#!/bin/sh

./clean.sh

python models_main.py -c /Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml -d dengue
python models_main.py -c /Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml -d wnv

cat logs/*.log
