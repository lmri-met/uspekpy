import pytest
import numpy as np

import uspeckpy.uspek as usp


class TestRandomUniform:
    def test_random_uniform(self):
        # Define the mean and the standard deviation
        loc, scale = 10, 0.01

        # Perform the random number generation
        result = usp.random_uniform(loc, scale)

        # Calculate the expected range
        low = loc * (1 - scale * np.sqrt(3))
        high = loc * (1 + scale * np.sqrt(3))

        # Assert that the result is within the expected range
        assert isinstance(result, float)
        assert result >= low
        assert result <= high


class TestRandomNormal:
    def test_random_normal(self):
        # Define the mean and the standard deviation
        loc, scale = 10, 0.01

        # Perform the random number generation
        result = usp.random_normal(loc, scale)

        # Calculate the expected range
        low = loc * (1 - 6 * scale)
        high = loc * (1 + 6 * scale)

        # Assert that the result is within the expected range
        assert isinstance(result, float)
        assert result >= low
        assert result <= high

# TODO: USpek
