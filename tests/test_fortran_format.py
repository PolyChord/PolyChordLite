"""Tests for the custom Fortran E24.15E3 and I12 formatters."""

import math
import numpy as np
import pytest
from pypolychord.polychord import _format_fortran_double, _format_fortran_int


class TestFormatFortranDouble:
    """Test _format_fortran_double against Fortran E24.15E3 output."""

    def test_field_width(self):
        """Every formatted value should be exactly 24 characters."""
        values = [1.0, -1.0, 0.0, 1e-100, 1e+300, float('nan'), float('inf')]
        for v in values:
            result = _format_fortran_double(v)
            assert len(result) == 24, f"Width {len(result)} for {v}: '{result}'"

    def test_positive_one(self):
        result = _format_fortran_double(1.0)
        assert result.strip() == '0.100000000000000E+001'

    def test_negative_pi(self):
        result = _format_fortran_double(-3.14159265358979)
        assert result.strip().startswith('-0.314159265358979')
        assert 'E+' in result

    def test_zero(self):
        result = _format_fortran_double(0.0)
        assert result.strip() == '0.000000000000000E+000'

    def test_negative_zero(self):
        result = _format_fortran_double(-0.0)
        assert result.strip() == '-0.000000000000000E+000'

    def test_small_value(self):
        result = _format_fortran_double(1e-100)
        # Fortran 0-based: 1e-100 = 0.1e-99
        assert 'E-099' in result
        assert result.strip().startswith('0.1000000000000')

    def test_large_value(self):
        result = _format_fortran_double(1e+300)
        assert 'E+301' in result
        assert result.strip().startswith('0.1000000000000')

    def test_nan(self):
        for nan_val in [float('nan'), np.float64('nan')]:
            result = _format_fortran_double(nan_val)
            assert len(result) == 24
            assert result.strip() == 'NaN'

    def test_inf(self):
        result = _format_fortran_double(float('inf'))
        assert len(result) == 24
        assert result.strip() == 'Infinity'

    def test_neg_inf(self):
        result = _format_fortran_double(float('-inf'))
        assert len(result) == 24
        assert result.strip() == '-Infinity'

    def test_logzero(self):
        """Test the typical logzero value used in PolyChord."""
        result = _format_fortran_double(-1e30)
        assert len(result) == 24
        assert result.strip().startswith('-0.1000000000000')
        assert 'E+031' in result

    def test_numpy_float64(self):
        """Ensure numpy float64 values work correctly."""
        result = _format_fortran_double(np.float64(2.5))
        assert len(result) == 24
        assert result.strip() == '0.250000000000000E+001'

    def test_normalization(self):
        """Mantissa should always start with 0. (Fortran convention)."""
        values = [1.0, 0.5, 123.456, 1e-10, 9.99]
        for v in values:
            result = _format_fortran_double(v).strip()
            # After optional sign, should start with '0.'
            if result.startswith('-'):
                assert result[1:3] == '0.', f"Bad normalization for {v}: {result}"
            else:
                assert result[0:2] == '0.', f"Bad normalization for {v}: {result}"


class TestFormatFortranInt:
    """Test _format_fortran_int against Fortran I12 output."""

    def test_field_width(self):
        result = _format_fortran_int(42)
        assert len(result) == 12

    def test_right_justified(self):
        result = _format_fortran_int(5)
        assert result == '           5'

    def test_negative(self):
        result = _format_fortran_int(-10)
        assert result == '         -10'

    def test_zero(self):
        result = _format_fortran_int(0)
        assert result == '           0'

    def test_numpy_int(self):
        result = _format_fortran_int(np.int64(100))
        assert result == '         100'
