import pytest
import models_main
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
# import vbdm
from utils import create_arg_parser
import sys

arg_list = ['config/pytest_config_1.yaml', 'config/pytest_config_2.yaml']
# arg_list = [('config/local_test_config.yaml', 'wnv'),
#             ('config/local_test_config.yaml', 'dengue')]


# DENGUE
@pytest.fixture
@pytest.mark.parametrize("config_file", arg_list)
def setup():
    #parser = create_arg_parser()
    #args, unknown = parser.parse_known_args()

    #config_file = args.config_file

    disease = DengueSEIRModel(config_file, args)

    disease.logger.info(disease)
    disease.run_model()
    disease.save_output('dengue')
    disease.logger.info('SUCCESS')

    return disease


# DENGUE
# @pytest.mark.skip
@pytest.mark.parametrize("config_file", arg_list)
def test_monkey(setup, monkeypatch, config_file):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'dengue'])
        # models_main.main()
        print("FLAG ----------", setup)

        assert True
