import pytest
#import models_main
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
def setup_dengue(config_file):
    disease = DengueSEIRModel(config_file)

    disease.logger.info(disease)
    disease.run_model()
    disease.save_output('dengue')
    disease.logger.info('SUCCESS')

    return disease


# WNV
@pytest.fixture
def setup_wnv(config_file):
    disease = WNVSEIRModel(config_file)

    disease.logger.info(disease)
    disease.run_model()
    disease.save_output('wnv')
    disease.logger.info('SUCCESS')

    return disease


# DENGUE
# @pytest.mark.skip
@pytest.mark.parametrize("config_file", arg_list)
def test_dengue(setup_dengue, monkeypatch, config_file):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'dengue'])
        # models_main.main()
        # print('\033[7m' + "FLAG ----------" + '\033[0m', setup_dengue.states)

        assert True


# WNV
# @pytest.mark.skip
@pytest.mark.parametrize("config_file", arg_list)
def test_wnv(setup_wnv, monkeypatch, config_file):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'wnv'])
        # models_main.main()
        # print('\033[7m' + "FLAG ----------" + '\033[0m', setup_wnv.states)

        assert True
