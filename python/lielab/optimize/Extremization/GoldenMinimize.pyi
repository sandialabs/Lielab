from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    import numpy.typing as npt
    import numpy as np
    from lielab.optimize import ExtremizationOptions

class GoldenMinimize:
    # status: ExtremizationStatus
    success: bool
    message: str

    def __init__(self) -> None: ...
    def __call__(self, objective: Callable[[npt.NDArray[np.fdloating]], float],
                       options: ExtremizationOptions) -> npt.NDArray[np.floating]: ...
