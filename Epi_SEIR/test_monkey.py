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
def setup():
    if not sys.argv[0].endswith('sphinx-build'):
        parser = create_arg_parser()
        args = parser.parse_args()

    config_file = args.config_file
    disease_name = args.disease_name.lower()

    if disease_name == 'dengue':
        den = DengueSEIRModel(config_file, args)
        den.logger.info(den)
        den.run_model()
        den.save_output(disease_name)
        den.logger.info('SUCCESS')
    elif disease_name == 'wnv':
        wnv = WNVSEIRModel(config_file, args)
        wnv.logger.info(wnv)
        wnv.run_model()
        wnv.save_output(disease_name)
        wnv.logger.info('SUCCESS')


# @pytest.mark.skip
@pytest.mark.parametrize("config_file, disease_name", arg_list)
def test_monkey(setup, monkeypatch, config_file, disease_name):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', disease_name])
        # models_main.main()
        assert True
