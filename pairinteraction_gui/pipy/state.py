"""Adding some methods to the StateOne and StateTwo classes of pairinteraction.
"""
import logging

from pipy import picomplex
from pipy import pireal

logger = logging.getLogger(__name__)


L_string = ["S", "P", "D", "F", "G", "H", "I", "J", "K", "L"]


def StateOne_getKey(self):
    """Returning a informative string describing the state."""
    if self.getJ() % 1 == 0:
        text = f"{self.getN()}{L_string[self.getL()]}{int(self.getJ())},{int(self.getM())}"
    elif (2 * self.getJ()) % 1 == 0:
        text = f"{self.getN()}{L_string[self.getL()]}{int(2*self.getJ())}/2,{int(2*self.getM())}/2"
    else:
        raise ValueError("J is not an integer or half-integer")
    return text


def StateTwo_getKey(self):
    """Returning a informative string describing the state."""
    return f"{self.getFirstState().getKey()};{self.getSecondState().getKey()}"


for pi in [picomplex, pireal]:
    pi.StateOne.getKey = StateOne_getKey
    pi.StateTwo.getKey = StateTwo_getKey
