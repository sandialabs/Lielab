from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import lielab

def left_Lie_group_action(g: lielab.domain.CompositeGroup, y: lielab.domain.CompositeManifold) -> lielab.domain.CompositeManifold: ...
def right_Lie_group_action(g: lielab.domain.CompositeGroup, y: lielab.domain.CompositeManifold) -> lielab.domain.CompositeManifold: ...
