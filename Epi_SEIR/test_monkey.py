import pytest
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
import sys

arg_list = ['config/unit_testing/pytest_config_1.yaml']
            #'config/unit_testing/pytest_config_2.yaml']

value_error_arglist = ['config/unit_testing/positive_duration.yaml',
                       'config/unit_testing/positive_resolution.yaml',
                       'config/unit_testing/duration_le_mosq.yaml',
                       'config/unit_testing/unique_position.yaml',
                       'config/unit_testing/positive_position.yaml',
                       'config/unit_testing/positive_states.yaml']
type_error_arglist = []

"""
TODO: have set of config files with things that should raise a type of error
such as IO errors, ValueError, etc. each of these error types gets it own test function
with the list of these config files that triggers them as an arglist. Have each one applicable
to both diseases.
"""


class TestDengue:
    @pytest.fixture
    def setup_dengue(self, config_file):
        disease = DengueSEIRModel(config_file)

        disease.logger.info(disease)
        disease.run_model()
        disease.save_output('dengue')
        disease.logger.info('SUCCESS')

        return disease

    # TESTING BLUEPRINT FUNCTION
    @pytest.mark.skip
    @pytest.mark.parametrize("config_file", arg_list)
    # def test_dengue(self, setup_dengue, monkeypatch, config_file):
    def test_dengue(self, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'dengue'])
            # print('\033[7m' + "FLAG ----------" + '\033[0m', setup_dengue.states)
            with pytest.raises(ValueError):
                disease = DengueSEIRModel(config_file)
            # models_main.main()
            # assert True

    @pytest.mark.parametrize("config_file", value_error_arglist)
    def test_value_error(self, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'dengue'])
            with pytest.raises(ValueError):
                disease = DengueSEIRModel(config_file)


class TestWNV:
    @pytest.fixture
    def setup_wnv(self, config_file):
        return 0
        disease = WNVSEIRModel(config_file)

        disease.logger.info(disease)
        disease.run_model()
        disease.save_output('wnv')
        disease.logger.info('SUCCESS')

        return disease

    # TESTING BLUEPRINT FUNCTION
    @pytest.mark.skip
    @pytest.mark.parametrize("config_file", arg_list)
    def test_wnv(self, setup_wnv, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'wnv'])
            # print('\033[7m' + "FLAG ----------" + '\033[0m', setup_wnv.states)

            assert True

    @pytest.mark.parametrize("config_file", value_error_arglist)
    def test_value_error(self, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'wnv'])
            with pytest.raises(ValueError):
                disease = WNVSEIRModel(config_file)
