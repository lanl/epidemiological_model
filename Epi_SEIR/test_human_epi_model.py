"""Unit Testing for Human Epi model code

Contains unit tests for CIMMID Human Epi code.

usage: pytest

"""

import pytest
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
import sys
import pandas as pd

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

param_dict_list_dengue = [{'nu_h': 0.1}, {'nu_v': 0.1}, {'gamma_h': 0.1}, {'mu_v': 0.05}, 
                          {'beta_h': 0.5}, {'beta_v': 0.5}, {'a_v': 0.5}, {'r_v': 0.1},
                          {'K_v': 200005}, {'nu_h': 0.1, 'nu_v': 0.1, 'gamma_h': 0.1, 'mu_v': 0.05,
                          'beta_h': 0.5, 'beta_v': 0.5, 'a_v': 0.5, 'r_v': 0.1, 'K_v': 200005}]
param_dict_list_wnv = [{'nu_b': 0.1}, {'nu_v': 0.15}, {'mu_b': 0.1}, {'mu_v': 0.05}, {'beta_b': 0.5}, 
                       {'alpha_b': 0.5}, {'alpha_v': 0.5}, {'eta': .005}, {'r_b': 0.005}, {'r_s':0.010},
                       {'K_b': 9500}, {'K_s': 8950}, {'nu_b': 0.1, 'nu_v': 0.15, 'mu_b': 0.1, 'mu_v': 0.05,
                       'beta_b': 0.5, 'alpha_b': 0.5, 'alpha_v': 0.5, 'eta': .005, 'r_b': 0.005, 'r_s':0.010,
                       'K_b': 9500, 'K_s': 8950}]


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
    
    #ask about the hard coding of the config file, and the sequencing of param_dict
    class TestModelOutput:
        """
            Tests for `run_model()` and corresponding output for DengueSEIRModel and DengueSEIRModel.param_dict.
        """

        def test_same_model_out(self):
            """
                For identical model runs, check that each output has correct dimensions and that the outputs are identical 
            """
            disease1 = DengueSEIRModel('config/local_test_config.yaml')
            disease1.run_model('dengue')
            disease2 = DengueSEIRModel('config/local_test_config.yaml')
            disease2.run_model('dengue')
            run1 = pd.DataFrame(dict(zip(list(disease1.state_names_order.values()), disease1.model_output.T)))
            run2 = pd.DataFrame(dict(zip(list(disease2.state_names_order.values()), disease2.model_output.T))) 
           
            #check that we have the correct number of columns and rows
            assert len(run1.columns) == len(disease1.initial_states)
            assert len(run2.columns) == len(disease2.initial_states)
            assert len(run1.index) == disease1.config_dict['DURATION'] * disease1.config_dict['RESOLUTION']
            assert len(run2.index) == disease2.config_dict['DURATION'] * disease2.config_dict['RESOLUTION']
            
            #check that both model run outputs are the same
            assert len(run1.columns) == len(run2.columns)
            assert sum(run1.columns == run2.columns) == len(run1.columns)
            assert len(run1.index) == len(run2.index)
            for k in range(0, len(run1.columns)):
                           assert sum(run1.iloc[:,k] == run2.iloc[:,k]) == len(run1.index)
        
        @pytest.mark.parametrize("param_dict", param_dict_list_dengue)
        def test_same_model_out_param_dict(self, param_dict):
            """
                For identical model runs of param_dict method, check that each output has correct dimensions and that the outputs are identical
            """
            disease1 = DengueSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease1.run_model('dengue')
            disease2 = DengueSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease2.run_model('dengue')
            run1 = pd.DataFrame(dict(zip(list(disease1.state_names_order.values()), disease1.model_output.T)))
            run2 = pd.DataFrame(dict(zip(list(disease2.state_names_order.values()), disease2.model_output.T)))
            
            #check that we have the correct number of columns and rows
            assert len(run1.columns) == len(disease1.initial_states)
            assert len(run2.columns) == len(disease2.initial_states)
            assert len(run1.index) == disease1.config_dict['DURATION'] * disease1.config_dict['RESOLUTION']
            assert len(run2.index) == disease2.config_dict['DURATION'] * disease2.config_dict['RESOLUTION']
            
            #check that both model run outputs are the same
            assert len(run1.columns) == len(run2.columns)
            assert sum(run1.columns == run2.columns) == len(run1.columns)
            assert len(run1.index) == len(run2.index)
            for k in range(0, len(run1.columns)):
                           assert sum(run1.iloc[:,k] == run2.iloc[:,k]) == len(run1.index)
        

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
          
    class TestModelOutput:
        """
            Tests for `run_model()` and corresponding output for WNVSEIRModel and WNVSEIRModel.param_dict.
        """                       
        def test_same_model_out(self):
            """
                For identical model runs, check that each output has correct dimensions and that the outputs are identical 
            """
            disease1 = WNVSEIRModel('config/local_test_config.yaml')
            disease1.run_model('wnv')
            disease2 = WNVSEIRModel('config/local_test_config.yaml')
            disease2.run_model('wnv')
            run1 = pd.DataFrame(dict(zip(list(disease1.state_names_order.values()), disease1.model_output.T)))
            run2 = pd.DataFrame(dict(zip(list(disease2.state_names_order.values()), disease2.model_output.T)))
            
            #check that we have the correct number of columns and rows
            assert len(run1.columns) == len(disease1.initial_states)
            assert len(run2.columns) == len(disease2.initial_states)
            assert len(run1.index) == disease1.config_dict['DURATION'] * disease1.config_dict['RESOLUTION']
            assert len(run2.index) == disease2.config_dict['DURATION'] * disease2.config_dict['RESOLUTION']
            
            #check that both model run outputs are the same
            assert len(run1.columns) == len(run2.columns)
            assert sum(run1.columns == run2.columns) == len(run1.columns)
            assert len(run1.index) == len(run2.index)
            #changed to -1 for now so Ih is not included
            for k in range(0, (len(run1.columns)-1)):
                           assert sum(run1.iloc[:,k] == run2.iloc[:,k]) == len(run1.index)
        
        @pytest.mark.parametrize("param_dict", param_dict_list_wnv)
        def test_same_model_out_param_dict(self, param_dict):
            """
                For identical model runs of param_dict method, check that each output has correct dimensions and that the outputs are identical
            """
            disease1 = WNVSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease1.run_model('wnv')
            disease2 = WNVSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease2.run_model('wnv')
            run1 = pd.DataFrame(dict(zip(list(disease1.state_names_order.values()), disease1.model_output.T)))
            run2 = pd.DataFrame(dict(zip(list(disease2.state_names_order.values()), disease2.model_output.T)))
            
            #check that we have the correct number of columns and rows
            assert len(run1.columns) == len(disease1.initial_states)
            assert len(run2.columns) == len(disease2.initial_states)
            assert len(run1.index) == disease1.config_dict['DURATION'] * disease1.config_dict['RESOLUTION']
            assert len(run2.index) == disease2.config_dict['DURATION'] * disease2.config_dict['RESOLUTION']
            
            #check that both model run outputs are the same
            assert len(run1.columns) == len(run2.columns)
            assert sum(run1.columns == run2.columns) == len(run1.columns)
            assert len(run1.index) == len(run2.index)
            #changed to -1 for now so Ih is not included
            for k in range(0, (len(run1.columns)-1)):
                           assert sum(run1.iloc[:,k] == run2.iloc[:,k]) == len(run1.index)
