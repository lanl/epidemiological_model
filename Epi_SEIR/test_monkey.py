import pytest
import models_main
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
# import vbdm
from utils import create_arg_parser
import sys

arg_list = [('config/local_test_config.yaml', 'wnv'),
            ('config/local_test_config.yaml', 'dengue')]


@pytest.fixture
# @pytest.mark.parametrize("config_file, disease_name", arg_list)
def setup(config_file, disease_name):
    #with monkeypatch.context() as m:
    #    m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', disease_name])

    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser()
        args, unknown = parser.parse_known_args()

    config_file = args.config_file
    disease_name = args.disease_name.lower()

    if disease_name == 'dengue':
        disease = DengueSEIRModel(config_file, args)
    elif disease_name == 'wnv':
        disease = WNVSEIRModel(config_file, args)

    disease.logger.info(disease)
    disease.run_model()
    disease.save_output(disease_name)
    disease.logger.info('SUCCESS')

    return disease


# @pytest.mark.skip
@pytest.mark.parametrize("config_file, disease_name", arg_list)
def test_monkey(setup, monkeypatch, config_file, disease_name):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', disease_name])
        # models_main.main()

        assert True
