from typing import TYPE_CHECKING, Final, Callable

if TYPE_CHECKING:
    import numpy.typing as npt
    import numpy as np
    import enum
    from lielab.domain import CompositeManifold, CompositeAlgebra


class ExtremizationMethod(enum.Enum):
    Undefined: Final[int]
    LineSearch: Final[int]
    Newton: Final[int]
    Golden: Final[int]

class ExtremizationOptions:
    method: ExtremizationMethod

    abstol: float
    reltol: float
    max_iterations: int

    # Golden options
    golden_ratio: float

    # LineSearch and GradientDescent options
    contraction_factor: float
    initial_alpha: float
    initial_step_size: float
    sufficient_decrease: float


class EuclideanExtremizationSystem:
    lower_bound: npt.NDArray[np.floating]
    upper_bound: npt.NDArray[np.floating]
    def __init__(self, objective: Callable[[npt.NDArray[np.floating]], float],
                       gradient: Callable[[npt.NDArray[np.floating]], npt.NDArray[np.floating]] | None = None,
                       hessian: Callable[[npt.NDArray[np.floating]], npt.NDArray[np.floating]] | None = None) -> None: ...
    def objective(self, x: npt.NDArray[np.floating]) -> float: ...
    def gradient(self, x: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]: ...
    def hessian(self, x: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]: ...
