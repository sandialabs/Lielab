from typing import TYPE_CHECKING, overload, Callable

if TYPE_CHECKING:
    import numpy.typing as npt
    import numpy as np
    from lielab.integrate import RungeKuttaCoefficients, EuclideanIVPSystem, IVPOptions, IVPSolution, IVPStatus

class RungeKutta:
    status: IVPStatus
    success: bool
    message: str
    abstol: float
    reltol: float
    error_estimate: float
    can_variable_step: bool
    implicit: bool
    order: int
    A: npt.NDArray[np.floating]
    B: npt.NDArray[np.floating]
    Bhat: npt.NDArray[np.floating]
    C: npt.NDArray[np.floating]
    e: npt.NDArray[np.floating]
    n: int
    K: npt.NDArray[np.floating]

    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, coefficients: RungeKuttaCoefficients) -> None: ...

    def __call__(self, vectorfield: Callable[[float, npt.NDArray[np.floating]], npt.NDArray[np.floating]],
                       y0: npt.NDArray[np.floating],
                       t0: float,
                       dt: float) -> npt.NDArray[np.floating]: ...

class RungeKuttaFlow:
    status: IVPStatus
    success: bool
    message: str
    tolerance_not_met: bool
    iterations: int

    def __init__(self) -> None: ...

    def __call__(self, dynamics: EuclideanIVPSystem,
                       tspan: npt.NDArray[np.floating],
                       y0: npt.NDArray[np.floating],
                       options: IVPOptions) -> IVPSolution: ...
