"""Unit Testing for Human Epi model code

Contains unit tests for CIMMID Human Epi code.

usage: pytest

"""

import pytest
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
import sys

value_error_arglist = ['config/unit_testing/positive_duration.yaml',
                       'config/unit_testing/positive_resolution.yaml',
                       #'config/unit_testing/duration_le_mosq.yaml',
                       'config/unit_testing/unique_position.yaml',
                       'config/unit_testing/positive_position.yaml',
                       'config/unit_testing/positive_states.yaml',
                       #'config/unit_testing/mosq_positive.yaml',
                       'config/unit_testing/output_types.yaml']
type_error_arglist = ['config/unit_testing/strings.yaml',
                      'config/unit_testing/position_integers.yaml',
                      'config/unit_testing/numerical_states.yaml',
                      #'config/unit_testing/mosq_numerical.yaml',
                      'config/unit_testing/output_is_string.yaml']


class TestDengue:
    """Defines a class to test dengue code.

    """

    @pytest.mark.parametrize("config_file", value_error_arglist)
    def test_value_error(self, monkeypatch, config_file):
        """tests ValueError exceptions.

        attributes:
            value_error_arglist: list of configuration file with values to trip ValueError.

        """
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'dengue'])
            with pytest.raises(ValueError):
                disease = DengueSEIRModel(config_file)

    @pytest.mark.parametrize("config_file", type_error_arglist)
    def test_type_error(self, monkeypatch, config_file):
        """tests TypeError exceptions.

        attributes:
            type_error_arglist: list of configuration file with values to trip TypeError.

        """
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'dengue'])
            with pytest.raises(TypeError):
                disease = DengueSEIRModel(config_file)


class TestWNV:
    """Defines a class to test WNV code.

    """

    @pytest.mark.parametrize("config_file", value_error_arglist)
    def test_value_error(self, monkeypatch, config_file):
        """Tests ValueError exceptions.

        Attributes:
            value_error_arglist: list of configuration file with values to trip ValueError.

        """
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'wnv'])
            with pytest.raises(ValueError):
                disease = WNVSEIRModel(config_file)

    @pytest.mark.parametrize("config_file", type_error_arglist)
    def test_type_error(self, monkeypatch, config_file):
        """tests TypeError exceptions.

        attributes:
            type_error_arglist: list of configuration file with values to trip TypeError.

        """
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', 'wnv'])
            with pytest.raises(TypeError):
                disease = WNVSEIRModel(config_file)
