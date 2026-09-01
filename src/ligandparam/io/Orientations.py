"""Deterministic orientation sets for multi-orientation RESP sampling.

``so3_n28`` packs 28 unit quaternions to cover SO(3). Because ``q`` and ``-q``
are the same rotation, pairwise distances use the geodesic
``2 * arccos(|q_i * q_j|)``.

``legacy_euler`` keeps the historical FreeLigand alpha/beta grid (gamma fixed
at 0) for reproducibility. Both protocols use 28 Gaussian ESP jobs.
"""

from typing import Union

import numpy as np

DEFAULT_ORIENTATION_PROTOCOL = "so3_n28"
N_ORIENTATIONS_SO3_N28 = 28

# Historical FreeLigand Euler grid (Rx/Ry sampling only; gamma stays 0).
LEGACY_EULER_ALPHA = [0, 30, 60, 90, 120, 150, 180]
LEGACY_EULER_BETA = [0, 30, 60, 90]
LEGACY_EULER_GAMMA = [0]

# Deterministic 28-point maximin quaternion pack. Identity is the reference
# orientation. Remaining points used farthest-point sampling from 300,000
# Haar-distributed unit-quaternion candidates (seed 20260805). Minimum
# pairwise SO(3) angle is approximately 62.71 degrees. Checked in so every
# platform submits identical Gaussian orientations.
SO3_N28_QUATERNIONS = np.array(
    [
        [1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [0.0000031356414280, -0.9095615636449691, 0.2509231135745325, 0.3312632684194684],
        [0.0017879157977430, 0.1247193741688564, 0.9262611244169686, -0.3556433753886636],
        [0.0078239493834159, 0.3871362166582298, 0.3029334583036826, 0.8708017313991565],
        [0.4961428959522769, 0.5863148435607335, -0.4462152018955720, 0.4593137540001764],
        [0.5041515145269712, 0.7145851118048574, 0.4786771515293781, 0.0778944991043741],
        [0.5005575786310763, 0.2041289980303654, -0.7313779544881572, -0.4157640560784830],
        [0.4952707781158514, -0.5704789551117755, 0.4384520406893776, -0.4868474361999007],
        [0.5092915194176429, -0.7094099404369925, -0.4815873222748422, -0.0737111638939350],
        [0.4964004176413882, -0.2102709024272367, 0.7370556519613690, 0.4075803465204031],
        [0.5003162543711515, -0.3128797119605562, -0.2060168875961550, 0.7806067982542133],
        [0.4975699786441238, 0.3325882924158541, 0.1981652112087048, -0.7762343030081844],
        [0.0097359825817809, 0.9221666839398029, 0.0297240522286650, 0.3855260022277908],
        [0.7112727963216501, 0.2811026109365733, 0.2005467327609979, 0.6122526760385053],
        [0.7099174745023861, 0.0879170535678644, 0.6512524057645994, -0.2532944434327783],
        [0.7052874294278562, -0.6501750443027605, 0.1920690451732453, 0.2072475224054270],
        [0.7093887416672131, -0.2521247939354007, -0.2076878794256816, -0.6245529971255467],
        [0.7062885310760540, 0.6433031022982426, -0.1743648927021977, -0.2385676290632025],
        [0.0107552290884569, 0.3537763202248029, 0.8600470594403679, 0.3674856403227599],
        [0.0144666775263194, -0.5711124821405287, 0.8201403619338488, -0.0314807037985041],
        [0.7107044596330199, -0.0740007821177396, -0.6563589557442694, 0.2421073656830626],
        [0.0138410467288641, -0.7372261886371829, -0.4721514926450980, 0.4830930968311887],
        [0.0233069522220613, 0.3746700516510961, -0.3880300856993315, -0.8417314244853006],
        [0.0005746331841135, -0.1831193286485443, 0.4311633198668584, -0.8834959948258958],
        [0.4899198715287050, -0.1869509736482237, 0.3321010434809690, 0.7840514969389000],
        [0.3700004339350381, -0.6722781449171541, -0.2946198108966455, -0.5695093869128571],
        [0.8432658219406016, -0.3551443875512570, -0.2221332322981222, 0.3367967408489471],
        [0.3056089818521138, -0.8899786809963406, 0.3205340293575676, -0.1086233566371907],
    ],
    dtype=float,
)


def quaternion_to_matrix(quaternion: Union[np.ndarray, list, tuple]) -> np.ndarray:
    """Convert a scalar-first unit quaternion ``(w, x, y, z)`` to a matrix.

    Parameters
    ----------
    quaternion : array-like, shape (4,)
        Quaternion to convert. It is normalized before conversion.

    Returns
    -------
    np.ndarray, shape (3, 3)
        Proper orthogonal rotation matrix.
    """
    q = np.asarray(quaternion, dtype=float)
    if q.shape != (4,):
        raise ValueError(f"Expected quaternion shape (4,), got {q.shape}")
    norm = np.linalg.norm(q)
    if np.isclose(norm, 0.0):
        raise ValueError("A zero quaternion does not define a rotation")
    w, x, y, z = q / norm
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


def get_quaternion_pack(name: str = "so3_n28") -> np.ndarray:
    """Return a copy of a named deterministic quaternion orientation pack.

    Parameters
    ----------
    name : {"so3_n28"}, optional
        Pack identifier.

    Returns
    -------
    np.ndarray, shape (n, 4)
        Scalar-first unit quaternions ``(w, x, y, z)``.
    """
    if name != "so3_n28":
        raise ValueError(f"Unknown quaternion orientation protocol: {name!r}")
    return SO3_N28_QUATERNIONS.copy()


def minimum_pairwise_rotation_angle(quaternions: np.ndarray, degrees: bool = True) -> float:
    """Return the minimum SO(3) geodesic angle among a quaternion set.

    Parameters
    ----------
    quaternions : np.ndarray, shape (n, 4)
        Unit (or near-unit) quaternions.
   degrees : bool, optional
        If True, return degrees; otherwise radians.
    """
    q = np.asarray(quaternions, dtype=float)
    if q.ndim != 2 or q.shape[1] != 4 or len(q) < 2:
        raise ValueError("Expected at least two quaternions with shape (n, 4)")
    q = q / np.linalg.norm(q, axis=1, keepdims=True)
    dots = np.abs(q @ q.T)
    np.fill_diagonal(dots, 0.0)
    angle = 2.0 * np.arccos(np.clip(np.max(dots), -1.0, 1.0))
    return float(np.degrees(angle) if degrees else angle)


def legacy_euler_kwargs() -> dict:
    """Keyword arguments for :class:`~ligandparam.stages.Gaussian.StageGaussianRotation` in legacy mode."""
    return {
        "orientation_protocol": "legacy_euler",
        "alpha": list(LEGACY_EULER_ALPHA),
        "beta": list(LEGACY_EULER_BETA),
        "gamma": list(LEGACY_EULER_GAMMA),
    }
