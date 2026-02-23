from typing import TYPE_CHECKING, Final, Tuple
import enum

if TYPE_CHECKING:
    import numpy.typing as npt
    import numpy as np

class RungeKuttaCoefficients(enum.Enum):
    FE1: Final[int]
    RK3: Final[int]
    RK4a: Final[int]
    RK4b: Final[int]
    RK5a: Final[int]
    RK5b: Final[int]
    RKF12a: Final[int]
    RKF12b: Final[int]
    RKF23a: Final[int]
    RKF23b: Final[int]
    RKF34a: Final[int]
    RKF34b: Final[int]
    RKF45a: Final[int]
    RKF45b: Final[int]
    RKF56: Final[int]
    RKF67: Final[int]
    RKF78: Final[int]
    RKF8: Final[int]
    RKDP54_7M: Final[int]
    RKV65e: Final[int]
    RKV65r: Final[int]
    RKV76e: Final[int]
    RKV76r: Final[int]
    RKV87e: Final[int]
    RKV87r: Final[int]
    RKV98e: Final[int]
    RKV98r: Final[int]
    BE1: Final[int]
    LG2: Final[int]
    LG4: Final[int]
    LG4s: Final[int]
    LG6: Final[int]
    LG6s: Final[int]
    Lobatto3A2: Final[int]
    Lobatto3A4: Final[int]
    Lobatto3A6: Final[int]

def get_butcher_tableau(method: RungeKuttaCoefficients) -> Tuple[npt.NDArray[np.floating], npt.NDArray[np.floating], npt.NDArray[np.floating], npt.NDArray[np.floating], npt.NDArray[np.floating], int, int, bool, bool]: ...

class CrouchGrossmanCoefficients(enum.Enum):
    CG23: Final[int]
    CG4a: Final[int]
    CG5a: Final[int]

def get_crouch_grossman_coefficients(method: CrouchGrossmanCoefficients) -> Tuple[npt.NDArray[np.floating], npt.NDArray[np.floating], npt.NDArray[np.floating], npt.NDArray[np.floating], npt.NDArray[np.floating], int, int, bool, bool]: ...
