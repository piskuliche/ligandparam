"""Tests for the charge normalization arithmetic in StageNormalizeCharge."""

import numpy as np
import pytest

from ligandparam.stages.charge import StageNormalizeCharge


@pytest.fixture
def stage(tmp_path):
    mol2 = tmp_path / "in.mol2"
    mol2.touch()
    return StageNormalizeCharge(
        "Normalize", main_input=mol2, cwd=tmp_path,
        out_mol2=tmp_path / "out.mol2", net_charge=0.0)


def test_defaults(stage):
    assert stage.precision == 0.0001
    assert stage.decimals == 4


def test_already_neutral_is_unchanged(stage):
    charges = np.array([0.5, -0.5])
    rounded, total, diff = stage.check_charge(charges)
    assert total == pytest.approx(0.0)
    assert diff == pytest.approx(0.0)


def test_normalize_reaches_the_target_net_charge(stage):
    charges = np.array([0.3, 0.3, -0.5005])
    rounded, total, diff = stage.check_charge(charges)
    adjusted = stage.normalize(rounded, diff)
    _, new_total, _ = stage.check_charge(adjusted)
    assert new_total == pytest.approx(0.0, abs=1e-9)


def test_normalize_conserves_the_atom_count(stage):
    charges = np.array([0.3, 0.3, -0.5005])
    rounded, _, diff = stage.check_charge(charges)
    assert len(stage.normalize(rounded, diff)) == len(charges)


def test_residual_below_one_precision_step_is_a_no_op(stage):
    """count rounds to 0 here. Dividing by it produced adjust=inf and then a
    misleading "Charge normalization failed"."""
    charges = np.array([0.1, -0.1])
    with np.errstate(divide="raise", invalid="raise"):
        result = stage.normalize(charges, 5e-8)
    assert np.all(np.isfinite(result))
    np.testing.assert_allclose(result, charges)


def test_nonzero_net_charge_target(tmp_path):
    mol2 = tmp_path / "in.mol2"
    mol2.touch()
    stage = StageNormalizeCharge(
        "Normalize", main_input=mol2, cwd=tmp_path,
        out_mol2=tmp_path / "out.mol2", net_charge=-1.0)
    charges = np.array([-0.4, -0.3, -0.2995])
    rounded, _, diff = stage.check_charge(charges)
    adjusted = stage.normalize(rounded, diff)
    _, new_total, _ = stage.check_charge(adjusted)
    assert new_total == pytest.approx(-1.0, abs=1e-9)
