"""Tests for ligandparam.io.coordinates helpers."""

import warnings

import numpy as np
import pytest

import MDAnalysis as mda

from ligandparam.io.coordinates import repair_zero_masses


def _universe(masses):
    u = mda.Universe.empty(len(masses), trajectory=True)
    u.add_TopologyAttr("masses", list(masses))
    u.atoms.positions = np.arange(len(masses) * 3, dtype=float).reshape(-1, 3)
    return u


def test_zero_masses_are_replaced():
    """AtomGroup.masses returns a copy, so an in-place masked assignment silently
    does nothing. This is what left center_of_mass() returning NaN."""
    u = _universe([0.0, 0.0, 12.0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        repair_zero_masses(u)
    assert not np.any(np.isclose(u.atoms.masses, 0, atol=0.1))


def test_center_of_mass_is_finite_after_repair():
    u = _universe([0.0, 0.0, 0.0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        repair_zero_masses(u)
    assert np.all(np.isfinite(u.atoms.center_of_mass()))


def test_real_masses_are_left_alone():
    u = _universe([12.011, 1.008, 15.999])
    original = u.atoms.masses.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        repair_zero_masses(u)
    np.testing.assert_allclose(u.atoms.masses, original)
