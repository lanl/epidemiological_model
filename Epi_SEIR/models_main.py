"""Main script to execute vector borne disease models simulation
Run from run_main.sh script. Instantiates VBDM (Vector Borne Disease Model)
given paramters from config.yaml configuration file.

    Typical usage example:

    python models_main.py -f <absolute path to config.yaml> -d <disease name
    [dengue][wnv]>
"""

from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
import VBDM


def main():
    config_file = VBDM.args.config_file
    disease_name = VBDM.args.disease_name.lower()

    if disease_name == 'dengue':
        den = DengueSEIRModel(config_file)
        den.run_model()
        den.save_output(disease_name)
        print("SUCCESS")
    elif disease_name == 'wnv':
        wnv = WNVSEIRModel(config_file)
        wnv.run_model()
        wnv.save_output(disease_name)
        print("SUCCESS")


if __name__ == "__main__":
    main()
