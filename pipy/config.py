"""Defining class Config.
This class is used to store all the configurations that are needed for the calculations.
"""
from hashlib import blake2b

import numpy as np

from pipy.misc import CustomDict

InvalidKey = object()


class Config:
    """Class for saving all the configurations that are needed for one pairinteraction calculation.
    This class is used for single atoms as well as pairstates.

    The idea is to provide one user dictionary (dic), which is sufficient to determine all the parameters,
    but might need some calculation/tweeking to get parameters used to initialize a pairinteraction calculation.
    Therefore all parameters, which are needed for the calculation are properties of this class and
    should be accessed via the properties and not via the user dictionary!
    You can also get a complete dict of all the parameters for the calculation via the output_dict.
    """

    def __init__(self, dic):
        """Initialize the Config class with the user dictionary dic.

        Args:
            dic: Dictionary with all the user parameters.
        """
        if not isinstance(dic, CustomDict):
            dic = CustomDict(dic, flattened=True)
        self._dic = dic  # private attribute, you should not access this directly from outside!

    def _get(self, key, *args):
        return self._dic.get(key, *args)

    @property
    def nAtoms(self):
        if hasattr(self, "force_nAtoms"):
            return self.force_nAtoms
        nAtoms = 1
        if self.useSameAtom or any(k in self._dic for k in ["atom2", "species2", "atom2.species"]):
            nAtoms = 2
        return self._get("nAtoms", nAtoms)

    @property
    def iAtom(self):
        if hasattr(self, "force_iAtom"):
            return self.force_iAtom
        return self._get("iAtom", 1)

    @property
    def isReal(self):
        isReal = all(atom.isReal for atom in self.atoms)
        return self._get("isReal", isReal)

    @property
    def method(self):
        return self._get("method", "NUMEROV").upper()

    @property
    def diagonalizeThreshold(self):
        return self._get("diagonalizeThreshold", 1e-10)

    @property
    def diamagnetism(self):
        return self._get("diamagnetism", False)

    @property
    def pathCache(self):
        return self._get("pathCache", "./cache/")

    @property
    def useSameAtom(self):
        return self._get("useSameAtom", self._get("samebasis", False))

    @property
    def atoms(self):
        if self.nAtoms == 1 or self.useSameAtom:
            return [self.atom1]
        return [self.atom1, self.atom2]

    @property
    def atom1(self):
        if not hasattr(self, "_atom1"):
            self._atom1 = Atom(self, 1)
        return self._atom1

    @property
    def atom2(self):
        if not hasattr(self, "_atom2"):
            if self.useSameAtom:
                self._atom2 = self.atom1
            else:
                self._atom2 = Atom(self, 2)
        return self._atom2

    @property
    def atom(self):
        return getattr(self, "atom" + str(self.iAtom))

    @property
    def pair(self):
        if not hasattr(self, "_pair"):
            self._pair = Pair(self)
        return self._pair

    def toDict(self):
        """Get a dictionary with all the output parameters.

        Returns:
            dict: Flat dictionary with all the output parameters.
        """
        output = {
            k: getattr(self, k)
            for k in ["nAtoms", "isReal", "method", "diagonalizeThreshold", "diamagnetism", "useSameAtom"]
        }
        for atom in self.atoms:
            output.update({f"atom{atom.iAtom}.{k}": v for k, v in atom.toDict().items()})
        if self.nAtoms == 2:
            output.update({f"pair.{k}": v for k, v in self.pair.toDict().items()})
        return output

    def toHash(self):
        """Return a hash of the toDict."""
        pairs = sorted((k, v) for k, v in self.toDict().items())
        h = blake2b(digest_size=10)
        h.update(str(pairs).encode("utf-8"))
        return h.hexdigest()

    def shallowCopy(self):
        return Config(self._dic)

    def deepCopy(self):
        return Config(self._dic.copy())


class Atom:
    """'Sub'Config class for one atom."""

    def __init__(self, config, i):
        self.iAtom = i
        self.config = config

    def _get(self, key, *args):
        return self.config._get(
            f"atom{self.iAtom}." + key,
            self.config._get("atom." + key, self.config._get(key + str(self.iAtom), self.config._get(key, *args))),
        )

    @property
    def isReal(self):
        isReal = self.Ey == 0 and self.By == 0
        return self._get("isReal", isReal)

    @property
    def diamagnetism(self):
        return self.config.diamagnetism

    @property
    def species(self):
        return self._get("species", None)

    @property
    def Efield(self):
        return self._get("Efield", [self._get("E" + x, 0) for x in "xyz"])

    @property
    def Ex(self):
        return self.Efield[0]

    @property
    def Ey(self):
        return self.Efield[1]

    @property
    def Ez(self):
        return self.Efield[2]

    @property
    def Bfield(self):
        return self._get("Bfield", [self._get("B" + x, 0) for x in "xyz"])

    @property
    def Bx(self):
        return self.Bfield[0]

    @property
    def By(self):
        return self.Bfield[1]

    @property
    def Bz(self):
        return self.Bfield[2]

    @property
    def qunumbers(self):
        if not hasattr(self, "_qunumbers"):
            self._qunumbers = self.getQunumbers()
        return self._qunumbers

    def getQunumbers(self):
        """Get the quantum numbers for the atom.

        Returns:
            list: List of quantum numbers for the atom.
            List of tuples (n, l, j, m).
        """
        qunumbers = []
        iAtom_list = [self.iAtom]
        if self.iAtom == 1 and self.config.useSameAtom:
            iAtom_list = [1, 2]
        for iAtom in iAtom_list:
            self.iAtom = iAtom
            qns = self._get("qunumbers", InvalidKey)
            if qns is not InvalidKey:
                shape = np.shape(qns)
                if shape == (4,):
                    qunumbers.append(qns)
                elif len(shape) == 2 and shape[1] == 4:
                    qunumbers += qns
                else:
                    raise ValueError("Atom qunumbers must be a list of 4-tuples or a Nx4 array.")
            else:
                n, l, j, m = (self._get(q, InvalidKey) for q in "nljm")
                if all(x is not InvalidKey for x in [n, l, j, m]):
                    qunumbers.append((n, l, j, m))
        self.iAtom = iAtom_list[0]
        qunumbers = [(int(qs[0]), int(qs[1]), float(qs[2]), float(qs[3])) for qs in qunumbers]
        return qunumbers

    @property
    def conserveMomenta(self):
        conserveMomenta = self._get("momenta", None) is not None
        return self._get("conserveMomenta", self._get("conserveM", conserveMomenta))

    @property
    def momenta(self):
        momenta = self._get("momenta", InvalidKey)
        if momenta is not InvalidKey:
            return momenta
        elif self.conserveMomenta:
            Ms = np.array(self.qunumbers)[:, 3]
            if not np.all(Ms == Ms[0]):
                raise ValueError("Not all atomstates have the same M, cannot conserve momenta.")
            return Ms[0]
        else:
            return None

    def restrictQN(self, Q):
        """Allowed formatting in config:
        - dQ, deltaQ, dQSingle, deltaQSingle
        - minQ, minQSingle, maxQ, maxQSingle
        """
        Q = Q.upper()
        dQ = self._get(
            "d" + Q,
            self._get("delta" + Q, self._get("delta" + Q + "Single", self._get("delta" + Q + "Single", InvalidKey))),
        )
        minQ = self._get("min" + Q, self._get("min" + Q + "Single", InvalidKey))
        maxQ = self._get("max" + Q, self._get("max" + Q + "Single", InvalidKey))

        if dQ is not InvalidKey and (minQ is not InvalidKey or maxQ is not InvalidKey):
            raise ValueError(f"d{Q} and min{Q}/max{Q} are given, only use one of them!")
        if minQ is InvalidKey or maxQ is InvalidKey:
            if dQ is InvalidKey or dQ is None:
                return None, None
            else:
                values = np.array(self.qunumbers)[:, ["N", "L", "J", "M"].index(Q)]
                minQ, maxQ = min(values) - dQ, max(values) + dQ

        if Q in "NL":
            return int(minQ), int(maxQ)
        else:  # Q in "JM"
            return float(minQ), float(maxQ)

    @property
    def minN(self):
        return self.restrictQN("N")[0]

    @property
    def maxN(self):
        return self.restrictQN("N")[1]

    @property
    def minL(self):
        return self.restrictQN("L")[0]

    @property
    def maxL(self):
        return self.restrictQN("L")[1]

    @property
    def minJ(self):
        return self.restrictQN("J")[0]

    @property
    def maxJ(self):
        return self.restrictQN("J")[1]

    @property
    def minM(self):
        return self.restrictQN("M")[0]

    @property
    def maxM(self):
        return self.restrictQN("M")[1]

    # no property for energies, because it is potentially complex to calculate
    def minEnergy(self, **kwargs):
        return self.restrictEnergy(**kwargs)[0]

    def maxEnergy(self, **kwargs):
        return self.restrictEnergy(**kwargs)[1]

    def restrictEnergy(self, **kwargs):
        if not hasattr(self, "_restrictEnergy"):
            dQ = self._get("dEnergy", self._get("deltaEnergy", self._get("deltaESingle", InvalidKey)))
            minQ = self._get("minEnergy", self._get("minEnergySingle", InvalidKey))
            maxQ = self._get("maxEnergy", self._get("maxEnergySingle", InvalidKey))

            if dQ is InvalidKey and (minQ is InvalidKey or maxQ is InvalidKey):
                raise ValueError("dEnergy and minEnergy/maxEnergy are given, only use one of them!")
            elif minQ is InvalidKey or maxQ is InvalidKey:
                if dQ is InvalidKey or dQ is None:
                    self._restrictEnergy = None, None
                else:
                    values = self.getEnergies(**kwargs)
                    self._restrictEnergy = min(values) - dQ, max(values) + dQ
            else:  # minQ and maxQ are given
                self._restrictEnergy = minQ, maxQ

        return self._restrictEnergy

    def getEnergies(self, StateOne=None, **kwargs):
        if not hasattr(self, "_energies"):
            if StateOne is None:
                from pipy import pireal

                StateOne = pireal.StateOne
            states = [StateOne(self.species, *qn) for qn in self.qunumbers]
            self._energies = [s.getEnergy() for s in states]
        return self._energies

    def setEfield(self, value):
        self.config._dic[f"atom{self.iAtom}.Efield"] = value

    def setBfield(self, value):
        self.config._dic[f"atom{self.iAtom}.Bfield"] = value

    def toDict(self):
        """Get a dictionary with all the output parameters.

        Returns:
            dict: Flat dictionary with all the output parameters.
        """
        output = {
            **{
                k: getattr(self, k)
                for k in [
                    "species",
                    "momenta",
                    "minN",
                    "maxN",
                    "minL",
                    "maxL",
                    "minJ",
                    "maxJ",
                    "minM",
                    "maxM",
                    "Bx",
                    "By",
                    "Bz",
                    "Ex",
                    "Ey",
                    "Ez",
                ]
            },
            "minEnergy": self.minEnergy(),
            "maxEnergy": self.maxEnergy(),
        }
        return output


class Pair:
    """'Sub'Config class for pair state."""

    def __init__(self, config):
        self.config = config

    def _get(self, key, *args):
        return self.config._get("pair." + key, self.config._get(key, *args))

    @property
    def order(self):
        return self._get("order", self._get("exponent", 3))

    @property
    def minimalNorm(self):
        return self._get("minimalNorm", 5e-2)

    @property
    def angle(self):
        return self._get("angle", self._get("theta", 0))

    @property
    def distance(self):
        return self._get("distance", self._get("R", np.inf))

    def symmetry(self, sym):
        val = self._get(sym, self._get(sym[:3], None))
        val = None if val is None else val.upper()
        dic = {
            None: [None, "NONE"],
            "ODD": ["ODD", "O"],
            "EVEN": ["EVEN", "E"],
        }
        for k, v in dic.items():
            if val in v:
                return k
        raise ValueError(f"Invalid symmetry value for {sym}: {val}")

    @property
    def inversion(self):
        return self.symmetry("inversion")

    @property
    def permutation(self):
        return self.symmetry("permutation")

    @property
    def reflection(self):
        return self.symmetry("reflection")

    @property
    def useQunumbersSinglePair(self):
        return self._get("useQunumbersSinglePair", False)

    @property
    def qunumbers(self):
        if not hasattr(self, "_qunumbers"):
            self._qunumbers = self.getQunumbers()
        return self._qunumbers

    def getQunumbers(self):
        """Get the quantum numbers for the atom.

        Returns:
            list: List of quantum numbers for the atom.
            List of tuples (n1, l1, j1, m1, n2, l2, j2, m2).
        """
        qunumbers = self._get("qunumbers", InvalidKey)
        if qunumbers is not InvalidKey:
            shape = np.shape(qunumbers)
            if shape == (8,):
                qunumbers = [qunumbers]
            elif len(shape) == 2 and shape[1] == 8:
                pass
            elif len(shape) == 3 and shape[1:] == (4, 2):
                # reshape [(N1, N2), (L1, L2), (J1, J2), (M1, M2)] to (N1, L1, J1, M1, N2, L2, J2, M2)
                qunumbers = np.reshape(qunumbers, (shape[0], 8), order="F")
            elif len(shape) == 3 and shape[1:] == (2, 4):
                # reshape [(N1, L1, J1, M1), (N2, L2, J2, M2)] to (N1, L1, J1, M1, N2, L2, J2, M2)
                qunumbers = np.reshape(qunumbers, (shape[0], 8), order="C")
            else:
                raise ValueError("Pair qunumbers must be a list of 8-tuples, Nx8, Nx4x2 or Nx2x4 array.")
        elif self.useQunumbersSinglePair:
            qunumbers = []
            for atom in self.config.atoms:
                qunumbers += atom.qunumbers
            if len(qunumbers) != 2:
                raise ValueError("When useQunumbersSinglePair there must be exactly two single atom qunumbers.")
            qunumbers = [(*qunumbers[0], *qunumbers[1])]
        else:
            qunumbers = []
            for i, qn1 in enumerate(self.config.atom1.qunumbers):
                for qn2 in self.config.atom2.qunumbers[i * self.config.useSameAtom :]:
                    qunumbers.append((*qn1, *qn2))

        qunumbers = [
            (int(qs[0]), int(qs[1]), float(qs[2]), float(qs[3]), int(qs[4]), int(qs[5]), float(qs[6]), float(qs[7]))
            for qs in qunumbers
        ]
        return qunumbers

    @property
    def conserveMomenta(self):
        conserveMomenta = self._get("momenta", None) is not None
        return self._get("conserveMomenta", self._get("conserveM", conserveMomenta))

    @property
    def momenta(self):
        momenta = self._get("momenta", InvalidKey)
        if momenta is not InvalidKey:
            return momenta
        elif self.conserveMomenta:
            Ms = np.array(self.qunumbers)[:, 3] + np.array(self.qunumbers)[:, 7]
            if not np.all(Ms == Ms[0]):
                raise ValueError("Not all atomstates have the same M, cannot conserve momenta.")
            return int(Ms[0])
        else:
            return None

    def restrictQN(self, Q):
        """Allowed formatting in config:
        - dQ, deltaQ, dQPair, deltaQPair
        - minQ, minQPair, maxQ, maxQPair
        """
        Q = Q.upper()
        dQ = self._get(
            "d" + Q,
            self._get("delta" + Q, self._get("delta" + Q + "Pair", self._get("delta" + Q + "Pair", InvalidKey))),
        )
        minQ = self._get("min" + Q, self._get("min" + Q + "Pair", InvalidKey))
        maxQ = self._get("max" + Q, self._get("max" + Q + "Pair", InvalidKey))

        if dQ is not InvalidKey and (minQ is not InvalidKey or maxQ is not InvalidKey):
            raise ValueError(f"d{Q} and min{Q}/max{Q} are given, only use one of them!")
        if minQ is InvalidKey or maxQ is InvalidKey:
            if dQ is InvalidKey or dQ is None or dQ < 0:
                return None, None
            else:
                raise NotImplementedError("Not sure if you should use this!")
                values = (
                    np.array(self.qunumbers)[:, ["N", "L", "J", "M"].index(Q)]
                    + np.array(self.qunumbers)[:, ["N", "L", "J", "M"].index(Q) + 4]
                )
                minQ, maxQ = min(values) - dQ, max(values) + dQ

        if Q in "NL":
            return int(minQ), int(maxQ)
        else:  # Q in "JM"
            return float(minQ), float(maxQ)

    # TODO add qun properties, min, max and limits for each qunumber?

    # no property for energies, because it is potentially complex to calculate
    def minEnergy(self, **kwargs):
        return self.restrictEnergy(**kwargs)[0]

    def maxEnergy(self, **kwargs):
        return self.restrictEnergy(**kwargs)[1]

    def restrictEnergy(self, **kwargs):
        if not hasattr(self, "_restrictEnergy"):
            dQ = self._get("dEnergy", self._get("deltaEnergy", self._get("deltaEPair", InvalidKey)))
            minQ = self._get("minEnergy", self._get("minEnergyPair", InvalidKey))
            maxQ = self._get("maxEnergy", self._get("maxEnergyPair", InvalidKey))

            if dQ is InvalidKey and (minQ is InvalidKey or maxQ is InvalidKey):
                raise ValueError("dEnergy and minEnergy/maxEnergy are given, only use one of them!")
            elif minQ is InvalidKey or maxQ is InvalidKey:
                if dQ is InvalidKey or dQ is None:
                    self._restrictEnergy = None, None
                else:
                    values = self.getEnergies(**kwargs)
                    self._restrictEnergy = min(values) - dQ, max(values) + dQ
            else:  # minQ and maxQ are given
                self._restrictEnergy = minQ, maxQ

        return self._restrictEnergy

    def getEnergies(self, atom1=None, atom2=None, **kwargs):
        if not hasattr(self, "_energies"):
            if atom1 is None or atom2 is None:
                raise NotImplementedError("For now atom1 and atom2 must be given to calculate Pair energies.")
                # TODO: implement this via import AtomOne ...
            qns1 = np.array(self.qunumbers)[:, :4]
            qns2 = np.array(self.qunumbers)[:, 4:]
            E1 = atom1.getEnergiesStates(qns1)
            E2 = atom2.getEnergiesStates(qns2)
            self._energies = list(E1 + E2)
        return self._energies

    def setSymmetry(self, sym, val):
        self.config._dic["pair." + sym] = val

    def setAngle(self, val):
        self.config._dic["pair.angle"] = val

    def setDistance(self, val):
        self.config._dic["pair.distance"] = val

    def toDict(self):
        """Get a dictionary with all the output parameters.

        Returns:
            dict: Flat dictionary with all the output parameters.
        """
        output = {
            **{
                k: getattr(self, k)
                for k in [
                    "order",
                    "minimalNorm",
                    "momenta",
                    #  "minN", "maxN", "minL", "maxL", "minJ", "maxJ", "minM", "maxM",
                    "inversion",
                    "permutation",
                    "reflection",
                    "angle",
                    "distance",
                ]
            },
            "minEnergy": self.minEnergy(),
            "maxEnergy": self.maxEnergy(),
        }
        return output
