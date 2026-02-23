from typing import TYPE_CHECKING, Final, Callable

if TYPE_CHECKING:
    import numpy.typing as npt
    import numpy as np
    import enum


class RootSolvingMethod(enum.Enum):
    Undefined: Final[int]
    Newton: Final[int]
    Hybrid: Final[int]

class RootSolvingOptions:
    method: RootSolvingMethod

    tol: float
    dx: float
    max_iterations: int

class EuclideanRootSystem:
    lower_bound: npt.NDArray[np.floating]
    upper_bound: npt.NDArray[np.floating]

    def __init__(self, objective: Callable[[npt.NDArray[np.floating]], npt.NDArray[np.floating]]) -> None: ...

    @staticmethod
    def objective(x: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]: ...
    @staticmethod
    def jacobian(x: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]: ...
