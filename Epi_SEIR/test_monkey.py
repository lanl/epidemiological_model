import sys
import pytest
import models_main


@pytest.mark.parametrize("config_file, disease_name", [('config/local_test_config.yaml', 'wnv'), ('config/local_test_config.yaml', 'dengue')])
def test_monkey(monkeypatch, config_file, disease_name):
    with monkeypatch.context() as m:
        m.setattr(sys, 'argv', ['models_main', '-c', config_file, '-d', disease_name])
        models_main.main()
        assert True
