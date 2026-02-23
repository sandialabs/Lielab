from typing import TYPE_CHECKING, Callable, Final, overload

if TYPE_CHECKING:
    import lielab.domain
    import numpy.typing as npt
    import numpy as np
    import enum
    from .Coefficients import *

class IVPMethod(enum.Enum):
    Undefined: Final[int]
    RungeKutta: Final[int]
    CrouchGrossman: Final[int]
    MuntheKaas: Final[int]


class IVPOptions:
    # IVP Meta-options
    method: IVPMethod

    # Multi-method options
    dt: float
    dt_min: float
    dt_max: float
    variable_time_step: bool

    reltol: float
    abstol: float

    small: float
    large: float
    pessimist: float

    max_iterations: int

    # Runge-Kutta specific options
    coefficients: RungeKuttaCoefficients

    # Crouch-Grossman specific options
    crouch_grossman_coefficients: CrouchGrossmanCoefficients

    # Munthe-Kaas specific options
    rebase_every_step: bool


class IVPStatus(enum.Enum):
    ERROR: Final[int]
    ERROR_MAX_ITERATIONS: Final[int]
    ERROR_INFS_IN_EVENT: Final[int]
    ERROR_NANS_IN_EVENT: Final[int]
    ERROR_INFS_IN_VF: Final[int]
    ERROR_NANS_IN_VF: Final[int]
    RUNNING: Final[int]
    SUCCESS: Final[int]
    SUCCESS_EVENT: Final[int]
    SUCCESS_BUT_TOL: Final[int]
    SUCCESS_EVENT_BUT_TOL: Final[int]


class IVPSolution:
    success: bool
    status: int
    message: str
    time_to_solution: float

    chuck_size: int
    current_index: int

    t: npt.NDArray[np.floating]
    y: list[lielab.domain.CompositeManifold]
    ybar: npt.NDArray[np.floating]
    theta: list[lielab.domain.CompositeAlgebra]
    thetabar: npt.NDArray[np.floating]

    debug: npt.NDArray[np.floating]

    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, other: IVPSolution) -> None: ...
    @overload
    def __init__(self, num_eoms: int) -> None: ...

    def copy(self) -> IVPSolution: ...
    @overload
    def trim_chunk(self) -> None: ...
    @overload
    def trim_chunk(self, last_index: int) -> None: ...
    def add_chunk(self) -> None: ...
    
    @overload
    def add_data(self, t_add: float, ybar_add: npt.NDArray[np.floating]) -> None: ...
    @overload
    def add_data(self, t_add: float, ybar_add: npt.NDArray[np.floating], thetabar_add: npt.NDArray[np.floating]) -> None: ...


class EuclideanIVPSystem:
    def __init__(self, vectorfield: Callable[[float, npt.NDArray[np.floating]], npt.NDArray[np.floating]],
                       event: Callable[[float, npt.NDArray[np.floating]], float] | None = None): ...
    def event(self, t: float, y: npt.NDArray[np.floating]) -> float: ...
    def vectorfield(self, t: float, y: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]: ...

class HomogeneousIVPSystem:
    def __init__(self, generator: Callable[[float, lielab.domain.CompositeManifold], lielab.domain.CompositeAlgebra],
                       action: Callable[[lielab.domain.CompositeGroup, lielab.domain.CompositeManifold], lielab.domain.CompositeManifold] | None = None,
                       connection: Callable[[lielab.domain.CompositeAlgebra, lielab.domain.CompositeAlgebra], lielab.domain.CompositeAlgebra] | None = None,
                       coordinates: Callable[[lielab.domain.CompositeAlgebra], lielab.domain.CompositeGroup] | None = None,
                       event: Callable[[float, lielab.domain.CompositeManifold], float] | None = None): ...
    def action(self, g: lielab.domain.CompositeGroup, y: lielab.domain.CompositeManifold) -> lielab.domain.CompositeManifold: ...
    def connection(self, theta: lielab.domain.CompositeAlgebra, xi: lielab.domain.CompositeAlgebra) -> lielab.domain.CompositeAlgebra: ...
    def coordinates(self, x: lielab.domain.CompositeAlgebra) -> lielab.domain.CompositeGroup: ...
    def event(self, t: float, y: lielab.domain.CompositeManifold) -> float: ...
    def generator(self, t: float, y: lielab.domain.CompositeManifold) -> lielab.domain.CompositeAlgebra: ...
