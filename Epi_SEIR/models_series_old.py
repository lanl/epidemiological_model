"""Script to execute vector borne disease models simulation through a series of parameter values
Instantiates VBDM (Vector Borne Disease Model)
given paramters from config.yaml configuration file.

    Typical usage example:

    python models_series.py -c <absolute path to config.yaml> -d <disease name
    [dengue][wnv]> -l <rename labels [True]>
"""
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
from utils import create_arg_parser
import sys
import numpy as np

# TODO remove this comment

def main_params():
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
    
    param_values = np.arange(disease.start_value, disease.stop_value, disease.step_value)
    
    for k in param_values:
        disease.params[disease.param_of_interest] = k
        disease.logger.info(disease)
        disease.run_model(disease_name)
        disease.save_output(disease_name, args.sim_labels)
        disease.logger.info('SUCCESS')


if __name__ == "__main__":
    main_params()
