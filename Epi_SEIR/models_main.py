"""Main script to execute vector borne disease models simulation
Run from run_main.sh script. Instantiates VBDM (Vector Borne Disease Model)
given paramters from config.yaml configuration file.

    Typical usage example:

    python models_main.py -f <absolute path to config.yaml> -d <disease name
    [dengue][wnv]>
"""

from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
import vbdm


def main():
    """Runs vector borne disease simulations.

    Instaniates each disease from configuration file
    """
    config_file = vbdm.args.config_file
    disease_name = vbdm.args.disease_name.lower()

    if disease_name == 'dengue':
        den = DengueSEIRModel(config_file)
        den.logger.info(den)
        den.run_model()
        den.save_output(disease_name)
        den.logger.info('SUCCESS')
    elif disease_name == 'wnv':
        wnv = WNVSEIRModel(config_file)
        wnv.run_model()
        wnv.save_output(disease_name)
        wnv.logger.info('SUCCESS')


if __name__ == "__main__":
    main()
