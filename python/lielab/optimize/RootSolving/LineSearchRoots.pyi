from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from lielab.optimize import EuclideanRootSystem, RootSolvingOptions

class LineSearchRoots:
    num_objective_evals: int
    iteration: int
    success: bool
    message: str

    def __init__(self) -> None: ...
    def __call__(self, problem: EuclideanRootSystem, x_guess: float, options: RootSolvingOptions) -> float: ...
