from typing import overload, TYPE_CHECKING, Any

if TYPE_CHECKING:
    from lielab.domain import CompositeManifold
    from lielab.integrate import EuclideanIVPSystem, HomogeneousIVPSystem, IVPOptions, IVPSolution
    import numpy.typing as npt
    import numpy as np


@overload
def solve_ivp(dynamics: EuclideanIVPSystem,
              tspan: npt.NDArray[np.floating[Any]] | list[float],
              y0: npt.NDArray[np.floating[Any]] | list[float],
              options: IVPOptions = IVPOptions()) -> IVPSolution: ...
@overload
def solve_ivp(dynamics: HomogeneousIVPSystem,
              tspan: npt.NDArray[np.floating[Any]] | list[float],
              y0: CompositeManifold,
              options: IVPOptions = IVPOptions()) -> IVPSolution: ...
