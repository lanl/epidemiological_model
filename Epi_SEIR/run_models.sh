#!/bin/sh

#if [ -f "logs/*.log" ]
#then
./clean_logs.sh
#fi

#python models_main.py -f /Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml -d wnv
python models_main.py -f /Users/jkeithley/Documents/CIMMID/human/dengue_model/Epi_SEIR/config/config.yaml -d dengue
