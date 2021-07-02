import pytest
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
import sys

arg_list = ['config/pytest_config_1.yaml',
            'config/pytest_config_2.yaml']


class TestDengue:
    @pytest.fixture
    def setup_dengue(self, config_file):
        disease = DengueSEIRModel(config_file)

        disease.logger.info(disease)
        disease.run_model()
        disease.save_output('dengue')
        disease.logger.info('SUCCESS')

        return disease

    # @pytest.mark.skip
    @pytest.mark.parametrize("config_file", arg_list)
    def test_dengue(self, setup_dengue, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'dengue'])
            # print('\033[7m' + "FLAG ----------" + '\033[0m', setup_dengue.states)

            assert True


class TestWNV:
    @pytest.fixture
    def setup_wnv(self, config_file):
        disease = WNVSEIRModel(config_file)

        disease.logger.info(disease)
        disease.run_model()
        disease.save_output('wnv')
        disease.logger.info('SUCCESS')

        return disease

    # @pytest.mark.skip
    @pytest.mark.parametrize("config_file", arg_list)
    def test_wnv(self, setup_wnv, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'wnv'])
            # print('\033[7m' + "FLAG ----------" + '\033[0m', setup_wnv.states)

            assert True
