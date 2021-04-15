"""Main script to execute vector borne disease models simulation
Run from run_main.sh script. Instantiates VBDM (Vector Borne Disease Model)
given paramters from config.yaml configuration file.

    Typical usage example:

    python models_main.py <absolute path to config.yaml>
"""

from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
import VBDM


def main():
    config_file = VBDM.args.config_file
    disease = VBDM.args.disease.lower()

    # logger = VBDM.create_logger('main', config_file)

    # print("DENGUE SUCCESS IS:", den.success)

    if disease == 'dengue':
        den = DengueSEIRModel(config_file, days=3)
        den.run_model()
    elif disease == 'wnv':
        wnv = WNVSEIRModel(config_file, days=3)
        wnv.run_model()

    # wnv.save_model()
    # print("WNV SUCCESS IS:", wnv.success)
    # den.graph_model()

    # if den.success and wnv.success:
    #    print("SUCCESS")
    # else:
    #    print("FAILURE")


if __name__ == "__main__":
    main()
