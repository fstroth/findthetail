import pytest
import numpy as np

from FindTheTail.ftt import Ftt


@pytest.fixture
def ftt_with_parameters():
    ftt = Ftt(np.arange(10), 'test_data', 100)
    return ftt


@pytest.fixture
def ftt_with_data():
    ftt = Ftt(np.arange(10), 'test_data')
    return ftt


@pytest.fixture
def ftt_data_with_dublicates():
    ftt = Ftt(np.ones(10, dtype='float'), 'test_data')
    return ftt


def test_initialisation_with_parameters(ftt_with_parameters):
    assert ftt_with_parameters.data_name == 'test_data'
    # order should be changed to descanding
    assert (ftt_with_parameters.data == np.arange(10)[::-1]).all()
    assert ftt_with_parameters.mc_steps == 100


def test_ftt_with_dublicates_in_data(ftt_data_with_dublicates):
    assert (np.abs(ftt_data_with_dublicates.data[:-1]-ftt_data_with_dublicates.data[1:]) > 0).all()


@pytest.mark.parametrize('number, expected', [
    (5., 0),
    (5.1, 1),
    (5.12, 2),
    (5.123, 3),
    (5.123456789, 9),
    (5.0001000, 4),
    (0000.0001, 4),
    (0.000, 0)
])
def test_get_significant_digit(number, expected):
    assert Ftt.get_significant_digit(number) == expected

