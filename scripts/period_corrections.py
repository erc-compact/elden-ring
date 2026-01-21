#!/usr/bin/env python3
"""
Central module for pulsar period correction calculations.

This module provides the canonical implementations of period correction
formulas used throughout the ELDEN-RING pipeline. All other code should
import from this module rather than implementing these formulas locally.

The period from acceleration searches (peasoup/PRESTO) is measured at the
center of the observation. These functions correct the period to the
appropriate epoch for different folding backends.

Reference: Originally from foldx.py

Usage:
    from period_corrections import a_to_pdot, period_correction_for_prepfold

    pdot = a_to_pdot(period, accel)
    fold_period = period_correction_for_prepfold(period, pdot, tsamp, fft_size)
"""

import numpy as np
from typing import Union

# Physical constants
LIGHT_SPEED = 2.99792458e8  # Speed of Light in SI (m/s)


def a_to_pdot(
    period_s: Union[float, np.ndarray],
    accel_ms2: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """
    Convert acceleration to period derivative.

    The relationship between acceleration and period derivative is:
        pdot = P * a / c

    where:
        P = period in seconds
        a = acceleration in m/s^2
        c = speed of light in m/s

    Args:
        period_s: Period in seconds (scalar or array)
        accel_ms2: Acceleration in m/s^2 (scalar or array)

    Returns:
        Period derivative (dimensionless, same shape as input)

    Example:
        >>> pdot = a_to_pdot(0.001, 10.0)  # 1ms pulsar, 10 m/s^2 acceleration
        >>> print(f"pdot = {pdot:.2e}")
        pdot = 3.34e-11
    """
    return period_s * accel_ms2 / LIGHT_SPEED


def period_correction_for_prepfold(
    p0: Union[float, np.ndarray],
    pdot: Union[float, np.ndarray],
    tsamp: float,
    fft_size: int
) -> Union[float, np.ndarray]:
    """
    Correct period to beginning of observation for prepfold.

    The period from peasoup/PRESTO acceleration searches is measured at the
    center of the FFT window (center of observation for full-length FFT).
    This function corrects it to the beginning of the observation, which is
    what prepfold expects when using the -topo flag.

    The correction formula is:
        P_corrected = P0 - pdot * (fft_size * tsamp) / 2

    This shifts the period from mid-observation to the start.

    Args:
        p0: Period in seconds at center of observation (scalar or array)
        pdot: Period derivative from a_to_pdot() (scalar or array)
        tsamp: Sampling time in seconds
        fft_size: FFT size used in the search (number of samples)

    Returns:
        Corrected period at beginning of observation (same shape as input)

    Example:
        >>> period = 0.001  # 1ms
        >>> pdot = a_to_pdot(period, 10.0)
        >>> tsamp = 64e-6   # 64 microseconds
        >>> fft_size = 2**27  # 134M samples
        >>> fold_period = period_correction_for_prepfold(period, pdot, tsamp, fft_size)

    Note:
        When using prepfold with the -topo flag, pass the corrected period
        and the pdot value directly:
            prepfold -topo -p {fold_period} -pd {pdot} -dm {dm} ...
    """
    return p0 - pdot * float(fft_size) * tsamp / 2


def period_correction_for_pulsarx(
    p0: Union[float, np.ndarray],
    pdot: Union[float, np.ndarray],
    no_of_samples: int,
    tsamp: float,
    fft_size: int
) -> Union[float, np.ndarray]:
    """
    Correct period for PulsarX folding.

    This correction accounts for the difference between the FFT size used
    in the search and the actual number of samples in the observation.

    The correction formula is:
        P_corrected = P0 - pdot * (fft_size - nsamples) * tsamp / 2

    Args:
        p0: Period in seconds at center of observation (scalar or array)
        pdot: Period derivative from a_to_pdot() (scalar or array)
        no_of_samples: Number of samples in the observation
        tsamp: Sampling time in seconds
        fft_size: FFT size used in the search

    Returns:
        Corrected period for PulsarX folding (same shape as input)

    Note:
        For PulsarX/psrfold_fil2, typically we use period_correction_for_prepfold
        because we want the epoch to point to tstart. This function is provided
        for cases where a different correction is needed.
    """
    return p0 - pdot * float(fft_size - no_of_samples) * tsamp / 2


def calculate_fold_parameters(
    period: Union[float, np.ndarray],
    accel: Union[float, np.ndarray],
    tsamp: float,
    fft_size: int,
    backend: str = 'prepfold'
) -> dict:
    """
    Calculate all folding parameters from search results.

    This is a convenience function that computes all parameters needed
    for folding in a single call.

    Args:
        period: Period in seconds from search (at center of observation)
        accel: Acceleration in m/s^2
        tsamp: Sampling time in seconds
        fft_size: FFT size used in search
        backend: Folding backend - 'prepfold' or 'pulsarx'

    Returns:
        Dictionary with keys:
            - fold_period: Corrected period for folding
            - pdot: Period derivative
            - frequency: Spin frequency (1/period)
            - fdot: Frequency derivative

    Example:
        >>> params = calculate_fold_parameters(0.001, 10.0, 64e-6, 2**27)
        >>> print(f"Fold period: {params['fold_period']:.10f}")
        >>> print(f"Pdot: {params['pdot']:.2e}")
    """
    pdot = a_to_pdot(period, accel)

    if backend == 'prepfold':
        fold_period = period_correction_for_prepfold(period, pdot, tsamp, fft_size)
    else:
        # For PulsarX, we typically use the same correction as prepfold
        # because we want pepoch to point to tstart
        fold_period = period_correction_for_prepfold(period, pdot, tsamp, fft_size)

    # Calculate frequency and frequency derivative
    frequency = 1.0 / fold_period if np.any(fold_period > 0) else 0.0
    fdot = -pdot / (period * period) if np.any(period > 0) else 0.0

    return {
        'fold_period': fold_period,
        'pdot': pdot,
        'frequency': frequency,
        'fdot': fdot
    }


# Alias for backwards compatibility
def pdot_from_accel(period_s, accel_ms2):
    """Alias for a_to_pdot() for backwards compatibility."""
    return a_to_pdot(period_s, accel_ms2)


if __name__ == '__main__':
    # Simple test/demo
    print("Period Corrections Module - Test")
    print("=" * 50)

    # Test parameters
    period = 0.001  # 1 ms pulsar
    accel = 10.0    # 10 m/s^2
    tsamp = 64e-6   # 64 microseconds
    fft_size = 2**27  # ~134M samples

    # Calculate
    pdot = a_to_pdot(period, accel)
    fold_period = period_correction_for_prepfold(period, pdot, tsamp, fft_size)

    print(f"Input period:     {period*1000:.6f} ms")
    print(f"Acceleration:     {accel:.2f} m/s^2")
    print(f"Sampling time:    {tsamp*1e6:.2f} us")
    print(f"FFT size:         {fft_size} samples")
    print()
    print(f"Period derivative: {pdot:.6e}")
    print(f"Corrected period:  {fold_period*1000:.6f} ms")
    print(f"Period change:     {(fold_period - period)*1e9:.3f} ns")

    # Test with array input
    print()
    print("Array test:")
    periods = np.array([0.001, 0.01, 0.1])
    accels = np.array([10.0, 5.0, 1.0])
    pdots = a_to_pdot(periods, accels)
    print(f"Periods: {periods}")
    print(f"Accels:  {accels}")
    print(f"Pdots:   {pdots}")
