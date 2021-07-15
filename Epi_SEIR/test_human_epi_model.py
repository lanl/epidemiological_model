import pytest
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
import sys

value_error_arglist = ['config/unit_testing/positive_duration.yaml',
                       'config/unit_testing/positive_resolution.yaml',
                       'config/unit_testing/duration_le_mosq.yaml',
                       'config/unit_testing/unique_position.yaml',
                       'config/unit_testing/positive_position.yaml',
                       'config/unit_testing/positive_states.yaml',
                       'config/unit_testing/mosq_positive.yaml',
                       'config/unit_testing/output_types.yaml']
type_error_arglist = ['config/unit_testing/strings.yaml',
                      'config/unit_testing/position_integers.yaml',
                      'config/unit_testing/numerical_states.yaml',
                      'config/unit_testing/mosq_numerical.yaml',
                      'config/unit_testing/output_is_string.yaml']

# TODO delete after gitlab ci is fixed
value_error_arglist = ['config/unit_testing/positive_duration.yaml',
                       'config/unit_testing/positive_states.yaml',
                       'config/unit_testing/mosq_positive.yaml']

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

        return disease

    # @pytest.mark.parametrize("config_file", ['config/local_test_config.yaml'])
    # def test_stuff():
    #     assert True

    # test ValueError response
    @pytest.mark.parametrize("config_file", value_error_arglist)
    def test_value_error(self, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'dengue'])
            with pytest.raises(ValueError):
                disease = DengueSEIRModel(config_file)

    # test TypeError response
    @pytest.mark.skip
    @pytest.mark.parametrize("config_file", type_error_arglist)
    def test_type_error(self, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'dengue'])
            with pytest.raises(TypeError):
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

    # test ValueError response
    @pytest.mark.skip
    @pytest.mark.parametrize("config_file", value_error_arglist)
    def test_value_error(self, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'wnv'])
            with pytest.raises(ValueError):
                disease = WNVSEIRModel(config_file)

    # test TypeError response
    @pytest.mark.skip
    @pytest.mark.parametrize("config_file", type_error_arglist)
    def test_type_error(self, monkeypatch, config_file):
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'wnv'])
            with pytest.raises(TypeError):
                disease = WNVSEIRModel(config_file)
