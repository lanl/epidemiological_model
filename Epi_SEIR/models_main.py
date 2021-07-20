"""Main script to execute vector borne disease models simulation
Run from run_main.sh script. Instantiates VBDM (Vector Borne Disease Model)
given paramters from config.yaml configuration file.

    Typical usage example:

    python models_main.py -f <absolute path to config.yaml> -d <disease name
    [dengue][wnv]>
"""

from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser
import sys

# TODO remove this comment

def main():
    """Runs vector borne disease simulations.

    Instaniates each disease from configuration file
    """
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser()
        args, unknown = parser.parse_known_args()

    disease_name = args.disease_name.lower()

    if disease_name == 'dengue':
        disease = DengueSEIRModel(args.config_file)
    elif disease_name == 'wnv':
        disease = WNVSEIRModel(args.config_file)

    disease.logger.info(disease)
    disease.run_model()
    disease.save_output(disease_name)
    disease.logger.info('SUCCESS')


if __name__ == "__main__":
    main()
