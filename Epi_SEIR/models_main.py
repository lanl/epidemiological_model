"""Main script to execute vector borne disease models simulation

Run from run_main.sh script. Instantiates VBDM (Vector Borne Disease Model)
given paramters from config.yaml configuration file.
    
    Typical usage example:

    python models_main.py <absolute path to config.yaml>
""" 

import logging
import os
import sys
import VBDM
import argparse

#path = os.environ['PATH']

def main():
    config_file = VBDM.args.config_file

    logger = VBDM.create_logger('main', config_file)

    den = VBDM.DengueSEIRModel(config_file, days=14)
    den.run_model()
    #den.graph_model()

if __name__ == "__main__":
    main()
