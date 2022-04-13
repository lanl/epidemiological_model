"""Main script to execute Logistic Link Model
"""

from LLM import LogisticLinkModel
from utils import create_arg_parser
import sys


def main():
    """Runs vector borne disease simulations.

    Instaniates each disease from configuration file
    """
    parser = create_arg_parser()
    args, unknown = parser.parse_known_args()

    mosq = LogisticLinkModel(args.config_file)

    mosq.run_model()
    mosq.save_output()
    mosq.logger.info('SUCCESS')


if __name__ == "__main__":
    main()
