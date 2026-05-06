"""Tests for xrd_profile.size_distributions — lognormal/normal fits
to W-A column-length A_size data."""
import numpy as np
import pytest

from xrd_profile.size_distributions import (
    fit_size_distribution, lognormal_a_size, normal_a_size)


class TestLognormalASize:
    def test_a_size_at_L0_is_one(self):
        """A_size(0) is normalised to 1 (Fourier coefficient at origin)."""
        D_median, sigma = 100.0, 0.3
        A0 = lognormal_a_size(np.array([0.0]), D_median, sigma)
        assert np.isclose(A0[0], 1.0, atol=1e-10)

    def test_a_size_decreases_monotonically(self):
        D_median, sigma = 100.0, 0.3
        L = np.linspace(0, 500, 100)
        A = lognormal_a_size(L, D_median, sigma)
        diffs = np.diff(A)
        # Monotonically decreasing (allow tiny numerical noise).
        assert np.all(diffs <= 1e-10)


class TestFitSizeDistribution:
    def test_recovers_lognormal_synthesis_within_5pct(self):
        D_median_synth, sigma_synth = 150.0, 0.4
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, D_median_synth, sigma_synth)
        # Add a touch of noise to simulate W-A output.
        rng = np.random.default_rng(seed=42)
        A_size_noisy = A_size + rng.normal(0, 0.005, len(L))

        result = fit_size_distribution(L, A_size_noisy)
        ln = result['lognormal']
        assert abs(ln['D_median'] - D_median_synth) / D_median_synth < 0.05
        assert abs(ln['sigma'] - sigma_synth) / sigma_synth < 0.10
        assert ln['fit_r2'] > 0.95

    def test_returns_none_when_too_few_valid_L(self):
        L = np.array([0, 50, 100])      # only 3 points -> too few
        A_size = np.array([1.0, 0.7, 0.4])
        result = fit_size_distribution(L, A_size)
        assert result is None

    def test_metadata_keys_present(self):
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, 100.0, 0.3)
        result = fit_size_distribution(L, A_size)
        assert set(result) >= {'lognormal', 'normal', 'method',
                                'initial_guess', 'n_valid_L'}
        assert result['method'] == 'curve_fit'
        assert result['initial_guess'] == 'moments'
        assert result['n_valid_L'] >= 4

    def test_lognormal_volume_mean_greater_than_median(self):
        """For a lognormal with sigma>0, D_mean_volume > D_median."""
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, 100.0, 0.4)
        result = fit_size_distribution(L, A_size)
        ln = result['lognormal']
        assert ln['D_mean_volume'] > ln['D_median']

    def test_normal_fit_returns_finite_for_lognormal_data(self):
        """Even when the data is lognormal, the normal fit should
        complete without error and return finite parameters with a
        lower R^2 than the lognormal fit."""
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, 100.0, 0.4)
        result = fit_size_distribution(L, A_size)
        nrm = result['normal']
        assert np.isfinite(nrm['D_mean'])
        assert np.isfinite(nrm['sigma'])
        assert result['lognormal']['fit_r2'] >= nrm['fit_r2']

    def test_negative_a_size_clipped(self):
        """A_size from W-A can dip slightly negative; fit should still
        complete (clipped to >= 0 in the basis function)."""
        L = np.linspace(0, 600, 60)
        A_size = lognormal_a_size(L, 100.0, 0.3)
        A_size[-5:] = -0.001  # noise dip at high L
        result = fit_size_distribution(L, A_size)
        assert result is not None
        assert np.isfinite(result['lognormal']['D_median'])

    def test_invalid_inputs_raise(self):
        with pytest.raises(ValueError):
            fit_size_distribution(np.array([1.0, 2.0]),
                                   np.array([1.0]))  # mismatched len
