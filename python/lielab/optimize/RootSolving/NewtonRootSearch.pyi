from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from lielab.optimize import EuclideanRootSystem, RootSolvingOptions
    import numpy.typing as npt
    import numpy as np

class NewtonRootSearch:
    iteration: int
    success: bool
    message: str

    def __init__(self) -> None: ...
    def __call__(self, problem: EuclideanRootSystem, x_guess: npt.NDArray[np.floating], options: RootSolvingOptions) -> npt.NDArray[np.floating]: ...
