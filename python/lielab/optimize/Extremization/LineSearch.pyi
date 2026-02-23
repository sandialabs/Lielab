from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    import numpy.typing as npt
    import numpy as np

class LineSearch:
    # status: ExtremizationStatus
    success: bool
    message: str

    alpha: float

    def __init__(self) -> None: ...

    # def __call__(self, objective: Callable[[npt.NDArray[np.floating]], float], x, y, z) -> npt.NDArray[np.floating]: ...
