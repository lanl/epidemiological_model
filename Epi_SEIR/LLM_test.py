import pytest
from LLM import LogisticLinkModel
import sys



value_error_arglist = ['config/unit_testing/date_length.yaml',
                       'config/unit_testing/date_nan.yaml',
                       'config/unit_testing/pop_active_nan.yaml',
                       'config/unit_testing/pop_total_nan.yaml',
                       'config/unit_testing/bad_data_file.yaml',
                       'config/unit_testing/data_file_not_csv.yaml',
                       'config/unit_testing/date_not_string.yaml',
                       'config/unit_testing/wrong_date_format.yaml',
                       'config/unit_testing/wrong_data_col_names.yaml',
                       'config/unit_testing/wrong_population_location_name.yaml',
                       'config/unit_testing/wrong_population_type_name.yaml']

class TestLLM:
    """Defines a class to test LogisticLinkModel code.

    """
    @pytest.mark.parametrize("config_file", value_error_arglist)
    def test_value_error(self, monkeypatch, config_file):
        """tests ValueError exceptions.

        attributes:
            value_error_arglist: list of configuration file with values to trip ValueError.

        """
        with monkeypatch.context() as m:
            m.setattr(sys, 'argv', ['models_main', '-c', config_file])
            with pytest.raises(ValueError):
                mosq = LogisticLinkModel(config_file)
