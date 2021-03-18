import logging
import os
import VBDM

def main():
    config_dir = 'config/'
    config_file = config_dir + 'config.yaml'

    logger = VBDM.create_logger(__name__, config_file)

    den = VBDM.DengueSEIRModel(config_file, days=14)
    den.run_model()
    #den.graph_model()

if __name__ == "__main__":
    main()
