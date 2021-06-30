"""Main script to execute vector borne disease models simulation
Run from run_main.sh script. Instantiates VBDM (Vector Borne Disease Model)
given paramters from config.yaml configuration file.

    Typical usage example:

    python models_main.py -f <absolute path to config.yaml> -d <disease name
    [dengue][wnv]>
"""

from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
#import vbdm
from utils import create_arg_parser
import sys
from logging import shutdown


def main():
    """Runs vector borne disease simulations.

    Instaniates each disease from configuration file
    """
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser()
        args, unknown = parser.parse_known_args()

    # config_file = args.config_file
    disease_name = args.disease_name.lower()

    if disease_name == 'dengue':
        den = DengueSEIRModel(args.config_file)
        den.logger.info(den)
        den.run_model()
        den.save_output(disease_name)
        den.logger.info('SUCCESS')
    elif disease_name == 'wnv':
        wnv = WNVSEIRModel(args.config_file)
        wnv.logger.info(wnv)
        wnv.run_model()
        wnv.save_output(disease_name)
        wnv.logger.info('SUCCESS')


if __name__ == "__main__":
    main()
