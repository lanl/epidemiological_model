import logging
import os
import sys
import VBDM

def main():
    if len(sys.argv) < 2:
        print("Please provide absolute path to root configuration file")
        exit()
    
    config_dir = sys.argv[1]
    config_file = config_dir + 'config.yaml'

    logger = VBDM.create_logger('main', config_file)

    den = VBDM.DengueSEIRModel(config_file, days=14)
    den.run_model()
    #den.graph_model()

if __name__ == "__main__":
    main()
