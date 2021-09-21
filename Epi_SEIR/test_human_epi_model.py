"""Unit Testing for Human Epi model code

Contains unit tests for CIMMID Human Epi code.

usage: pytest

"""

import pytest
from dengue import DengueSEIRModel
from wnv import WNVSEIRModel
import sys
import pandas as pd
import numpy as np

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

# Create parameter value lise to sequence through
def gen_new_params(disease_name):
    if disease_name == 'dengue':
        params = DengueSEIRModel('config/local_test_config.yaml').params
    elif disease_name == 'wnv':
        params = WNVSEIRModel('config/local_test_config.yaml').params
    rng = np.random.default_rng()
    scalars = rng.uniform(low = .75, high = 1.25, size = len(params))
    
    #make sure none of the scalars are 1 so the output is different
    for i in range(0, len(scalars)):
        if (.99 < scalars[i] < 1.01) == True:
            scalars[i] = 1.011
            
    param_dict_list = []
    for i in range(0, len(params)):
        #make the change in eta big enough that it will create a different model output
        #will need to change this if we get rid of the poisson distribution
        if disease_name == 'wnv' and list(params.keys())[i] == 'eta':
            scalar = rng.uniform(low = 100, high = 1000, size = 1)
            val = {list(params.keys())[i]: (list(params.values())[i] * scalar)[0]}
        #r_v is not very sensitive, so sometimes getting the same output even when changed slightly: this is to fix that
        elif disease_name == 'dengue' and list(params.keys())[i] == 'r_v':
            scalar = rng.uniform(low = 2, high = 10, size = 1)
            val = {list(params.keys())[i]: (list(params.values())[i] * scalar)[0]}
        else:
            val = {list(params.keys())[i]: (list(params.values()) * scalars)[i]}
        param_dict_list.append(val)

    param_dict_list.append(dict(zip(list(params.keys()), list(params.values()) * scalars)))
    return(param_dict_list)

param_dict_list_dengue = gen_new_params('dengue')
param_dict_list_wnv = gen_new_params('wnv')

dengue = DengueSEIRModel('config/local_test_config.yaml')
wnv = WNVSEIRModel('config/local_test_config.yaml')

eq_points_dengue = [{'Sh': dengue.initial_states['Sh'], 'Eh': 0, 'Ih':0, 'Rh': 0, 'Sv': 0, 'Ev': 0, 'Iv': 0}, {'Sh': dengue.initial_states['Sh'], 'Eh': 0, 'Ih':0, 'Rh': 0, 'Sv': dengue.params['K_v'], 'Ev': 0, 'Iv': 0}]
eq_points_wnv = [{'Sv': 0, 'Ev': 0, 'Iv': 0, 'Sb': wnv.initial_states['Sb'], 'Eb': 0, 'Ib': 0, 'Rb': 0, 'Ih': 0}, {'Sv': wnv.params['K_v'], 'Ev': 0, 'Iv': 0, 'Sb': wnv.initial_states['Sb'], 'Eb': 0, 'Ib': 0, 'Rb': 0, 'Ih': 0}]


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

        def test_model_out(self):
            """
                For identical model runs, check that each output has correct dimensions, that the outputs are identical, and that the output changes over time.
            """
            disease1 = DengueSEIRModel('config/local_test_config.yaml')
            disease1.run_model('dengue')
            disease2 = DengueSEIRModel('config/local_test_config.yaml')
            disease2.run_model('dengue')
            run1 = pd.DataFrame(dict(zip(list(disease1.state_names_order.values()), disease1.model_output.T)))
            run2 = pd.DataFrame(dict(zip(list(disease2.state_names_order.values()), disease2.model_output.T))) 
           
            #check that we have the correct number of columns and rows
            assert len(run1.columns) == len(disease1.initial_states)
            assert len(run1.index) == disease1.config_dict['DURATION'] * disease1.config_dict['RESOLUTION']
            
            #check that both model run outputs are the same
            assert len(run1.columns) == len(run2.columns)
            assert sum(run1.columns == run2.columns) == len(run1.columns)
            assert len(run1.index) == len(run2.index)
            
            col_names = list(run1.columns)
            for k in col_names:
                assert sum(run1[k] == run2[k]) == len(run1.index)
                #make sure output changes not at equilibrium
                assert round(sum(run1[k].diff().iloc[1:,]),3) != 0.000
        
        @pytest.mark.parametrize("param_dict", param_dict_list_dengue)
        def test_model_out_param_dict(self, param_dict):
            """
                For identical model runs of param_dict method, check that each output has correct dimensions, that the outputs are identical, and that the      output changes over time.\n
                Check that param_dict model run is different from standard model run.
            """
            disease1 = DengueSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease1.run_model('dengue')
            disease2 = DengueSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease2.run_model('dengue')
            run1 = pd.DataFrame(dict(zip(list(disease1.state_names_order.values()), disease1.model_output.T)))
            run2 = pd.DataFrame(dict(zip(list(disease2.state_names_order.values()), disease2.model_output.T)))
            
            #add normal run to compare against
            disease = DengueSEIRModel('config/local_test_config.yaml')
            disease.run_model('dengue')
            norm_run = pd.DataFrame(dict(zip(list(disease.state_names_order.values()), disease.model_output.T)))
            
            #check that we have the correct number of columns and rows
            assert len(run1.columns) == len(disease1.initial_states)
            assert len(run1.index) == disease1.config_dict['DURATION'] * disease1.config_dict['RESOLUTION']
            
            #check that both model run outputs are the same
            assert len(run1.columns) == len(run2.columns)
            assert sum(run1.columns == run2.columns) == len(run1.columns)
            assert len(run1.index) == len(run2.index)
            
            #check identical runs are the same, and param_dict run different from normal: having issues with the latter, parameter 7
            col_names = list(run1.columns)
            out_sums = []
            for k in col_names:
                assert sum(run1[k] == run2[k]) == len(run1.index)
                assert round(sum(run1[k].diff().iloc[1:,]),3) != 0.000
                out_sums.append(sum(norm_run[k] == run1[k]))
            assert sum(out_sums) < len(out_sums)*500
        
        @pytest.mark.parametrize("eq", eq_points_dengue)
        def test_eq_points(self, eq):
            """
                For model runs with initial states as equilibrium points, check that output does not change at all.
            """
            disease = DengueSEIRModel('config/local_test_config.yaml')
            disease.initial_states = eq
            disease.run_model('dengue')
            run = pd.DataFrame(dict(zip(list(disease.state_names_order.values()), disease.model_output.T)))
            
            col_names = list(run.columns)
            for k in col_names:
                #rounding due to returning very small numbers
                assert round(sum(run[k].diff().iloc[1:,]),3) == 0.000
        
        @pytest.mark.parametrize("eq", eq_points_dengue)
        @pytest.mark.parametrize("param_dict", param_dict_list_dengue)
        def test_eq_points_param_dict(self, eq, param_dict):
            """
                For model runs of param_dict method with initial states as equilibrium points, check that output does not change at all.
            """
            
            disease = DengueSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease.initial_states = eq
            #need to change eq where Sv = K_v if K_v changes
            if disease.initial_states['Sv'] != 0 and disease.initial_states['Sv'] != disease.params['K_v']:
                disease.initial_states['Sv'] = disease.params['K_v']
                
            disease.run_model('dengue')
            run = pd.DataFrame(dict(zip(list(disease.state_names_order.values()), disease.model_output.T)))
            
            col_names = list(run.columns)
            for k in col_names:
                assert round(sum(run[k].diff().iloc[1:,]),3) == 0.000
        

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

        def test_model_out(self):
            """
                For identical model runs, check that each output has correct dimensions, that the outputs are identical, and that the output changes over time.
            """
            disease1 = WNVSEIRModel('config/local_test_config.yaml')
            disease1.run_model('wnv')
            disease2 = WNVSEIRModel('config/local_test_config.yaml')
            disease2.run_model('wnv')
            run1 = pd.DataFrame(dict(zip(list(disease1.state_names_order.values()), disease1.model_output.T)))
            run2 = pd.DataFrame(dict(zip(list(disease2.state_names_order.values()), disease2.model_output.T))) 
           
            #check that we have the correct number of columns and rows
            assert len(run1.columns) == len(disease1.initial_states)
            assert len(run1.index) == disease1.config_dict['DURATION'] * disease1.config_dict['RESOLUTION']
            
            #check that both model run outputs are the same
            assert len(run1.columns) == len(run2.columns)
            assert sum(run1.columns == run2.columns) == len(run1.columns)
            assert len(run1.index) == len(run2.index)
            
            col_names = list(run1.columns)
            col_names.remove('Infected Humans')
            for k in col_names:
                assert sum(run1[k] == run2[k]) == len(run1.index)
                #make sure output changes not at equilibrium
                assert round(sum(run1[k].diff().iloc[1:,]),3) != 0.000
        
        @pytest.mark.parametrize("param_dict", param_dict_list_wnv)
        def test_model_out_param_dict(self, param_dict):
            """
                For identical model runs of param_dict method, check that each output has correct dimensions, that the outputs are identical, and that the      output changes over time.\n
                Check that param_dict model run is different from standard model run.
            """
            disease1 = WNVSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease1.run_model('wnv')
            disease2 = WNVSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease2.run_model('wnv')
            run1 = pd.DataFrame(dict(zip(list(disease1.state_names_order.values()), disease1.model_output.T)))
            run2 = pd.DataFrame(dict(zip(list(disease2.state_names_order.values()), disease2.model_output.T)))
            
            #add normal run to compare against
            disease = WNVSEIRModel('config/local_test_config.yaml')
            disease.run_model('wnv')
            norm_run = pd.DataFrame(dict(zip(list(disease.state_names_order.values()), disease.model_output.T)))
            
            #check that we have the correct number of columns and rows
            assert len(run1.columns) == len(disease1.initial_states)
            assert len(run1.index) == disease1.config_dict['DURATION'] * disease1.config_dict['RESOLUTION']
            
            #check that both model run outputs are the same
            assert len(run1.columns) == len(run2.columns)
            assert sum(run1.columns == run2.columns) == len(run1.columns)
            assert len(run1.index) == len(run2.index)
            
            #check identical runs are the same, and param_dict run different from normal: having issues with the latter, parameter 7
            col_names = list(run1.columns)
            col_names.remove('Infected Humans')
            out_sums = []
            for k in col_names:
                assert sum(run1[k] == run2[k]) == len(run1.index)
                assert round(sum(run1[k].diff().iloc[1:,]),3) != 0.000
                out_sums.append(sum(norm_run[k] == run1[k]))
            assert sum(out_sums) < len(out_sums)*500
        
        @pytest.mark.parametrize("eq", eq_points_wnv)
        def test_eq_points(self, eq):
            """
                For model runs with initial states as equilibrium points, check that output does not change at all.
            """
            disease = WNVSEIRModel('config/local_test_config.yaml')
            disease.initial_states = eq
            disease.run_model('wnv')
            run = pd.DataFrame(dict(zip(list(disease.state_names_order.values()), disease.model_output.T)))
            
            col_names = list(run.columns)
            for k in col_names:
                #rounding due to returning very small numbers
                assert round(sum(run[k].diff().iloc[1:,]),3) == 0.000
        
        @pytest.mark.parametrize("eq", eq_points_wnv)
        @pytest.mark.parametrize("param_dict", param_dict_list_wnv)
        def test_eq_points_param_dict(self, eq, param_dict):
            """
                For model runs of param_dict method with initial states as equilibrium points, check that output does not change at all.
            """
            
            disease = WNVSEIRModel.param_dict('config/local_test_config.yaml', param_dict)
            disease.initial_states = eq
            #need to change eq where Sv = K_v if K_v changes
            if disease.initial_states['Sv'] != 0 and disease.initial_states['Sv'] != disease.params['K_v']:
                disease.initial_states['Sv'] = disease.params['K_v']
                
            disease.run_model('wnv')
            run = pd.DataFrame(dict(zip(list(disease.state_names_order.values()), disease.model_output.T)))
            
            col_names = list(run.columns)
            for k in col_names:
                assert round(sum(run[k].diff().iloc[1:,]),3) == 0.000
        