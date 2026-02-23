from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from lielab.optimize import EuclideanExtremizationSystem, ExtremizationOptions
    import numpy.typing as npt
    import numpy as np

class NewtonMinimize:
    iteration: int
    success: bool
    message: str

    def __init__(self) -> None: ...
    def __call__(self, problem: EuclideanExtremizationSystem, x_guess: npt.NDArray[np.floating], options: ExtremizationOptions) -> npt.NDArray[np.floating]: ...
