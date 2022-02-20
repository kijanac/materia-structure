from __future__ import annotations
from typing import Generator, Iterable, Optional, TextIO, Union, TYPE_CHECKING

import contextlib
import functools
import itertools
from materia.utils import expand, rotation_matrix
import numpy as np
from openbabel import openbabel as ob
from openbabel import pybel
import pathlib
import pubchempy as pcp
import scipy.linalg
import tempfile
import unyt

if TYPE_CHECKING:
    import cclib
    from materia.schema import Schema

__all__ = ["Structure"]


class Structure:
    """Atomistic structure.

    Parameters
    ----------
    pybelmol : pybel.Molecule
        Pybel Molecule to use as internal structure representation.
    """

    def __init__(self, pybelmol: pybel.Molecule) -> None:
        self.pybelmol = pybelmol

    @staticmethod
    def from_obmol(obmol: ob.OBMol) -> Structure:
        """Build structure from Openbabel OBMol.

        Parameters
        ----------
        obmol : ob.OBMol
            Openbabel OBMol.

        Returns
        -------
        Structure
            Structure from OBMol.
        """
        return Structure(pybel.Molecule(obmol))

    @staticmethod
    def from_coordinates(
        elements: Iterable[Union[str, int]], coordinates: Iterable[Union[int, float]]
    ) -> Structure:
        """Build structure from coordinates.

        Parameters
        ----------
        elements : Iterable[Union[str, int]]
            Elements specified by atomic number or atomic symbols.
        coordinates : Iterable[Union[int, float]]
            Atomic coordinates given as a flat iterable.

        Returns
        -------
        Structure
            Structure from coordinates.
        """
        obmol = ob.OBMol()
        print(coordinates)

        for Z, (x, y, z) in zip(
            elements, np.split(np.array(coordinates), len(elements))
        ):
            obatom = ob.OBAtom()
            if isinstance(Z, str):
                obatom.SetAtomicNum(ob.GetAtomicNum(Z))
            else:
                obatom.SetAtomicNum(int(Z))
            obatom.SetVector(x, y, z)
            obmol.AddAtom(obatom)

        obmol.ConnectTheDots()
        obmol.PerceiveBondOrders()

        return Structure.from_obmol(obmol)

    @staticmethod
    def from_cclib(parsed: cclib.parser.ccData) -> Structure:
        """Build structure from cclib parsed data.

        Parameters
        ----------
        parsed : cclib.parser.ccData
            Parsed data from cclib.

        Returns
        -------
        Structure
            Structure from parsed data.
        """
        return Structure.from_coordinates(
            parsed.atomnos, parsed.atomcoords[-1].flatten()
        )

    @staticmethod
    def from_schema(s: Schema) -> Structure:
        """Build structure from Materia schema.

        Parameters
        ----------
        s : materia.schema.Schema
            Schema.

        Returns
        -------
        Structure
            Structure from schema.
        """
        return Structure.from_coordinates(s.molecule.symbols, s.molecule.geometry)

    @staticmethod
    def generate(
        identifier_type: str,
        identifier: str,
        forcefield: Optional[str] = "mmff94",
        steps: Optional[int] = 500,
    ) -> Structure:
        """Generate structure from string identifier.

        Parameters
        ----------
        identifier_type : str
            Identifier type.
        identifier : str
            Molecular identifier.
        forcefield : str, optional
            Forcefield for structure optimization.
            By default `mmff94`.
        steps : int, optional
            Number of structure optimization steps.
            By default 500.

        Returns
        -------
        Structure
            Generated structure.
        """
        pybelmol = pybel.readstring(identifier_type, identifier)
        pybelmol.addh()
        pybelmol.localopt(forcefield, steps)

        return Structure(pybelmol)

    @staticmethod
    def read(
        filepath: Union[str, pathlib.Path], filetype: Optional[str] = None
    ) -> Structure:
        """Read structure from a file.

        Parameters
        ----------
        filepath : Union[str, pathlib.Path]
            Path to file from which the structure will be read.
        filetype : str, optional
            File extension.
            If `None`, inferred from `filepath`.
            By default `None`.

        Returns
        -------
        Structure
            Structure from file.
        """
        filetype = filetype or pathlib.Path(filepath).suffix[1:]
        pybelmol = next(pybel.readfile(filetype, expand(filepath)))
        return Structure(pybelmol)

    @staticmethod
    def retrieve(
        identifier_type: str,
        identifier: str,
        use_first_match: Optional[bool] = False,
    ) -> Structure:
        """Retrieve structure from PubChemPy database using identifier.

        Parameters
        ----------
        identifier_type : str
            Identifier type.
        identifier : str
            Molecular identifier.
        use_first_match : bool, optional
            Use first match if multiple matches are found.
            By default False.

        Returns
        -------
        Structure
            Retrieved structure.

        Raises
        ------
        ValueError
        """
        if identifier_type not in ("name", "smiles", "inchi", "inchikey"):
            raise ValueError(
                "Provide name, SMILES, InChi, or InChiKey to retrieve structure."
            )

        try:
            cids = pcp.get_cids(identifier, identifier_type)
        except OSError:
            raise ValueError(f"Structure retrieval for {identifier} failed.")

        if not (use_first_match or len(cids) == 1):
            msg = f"Structure retrieval for {identifier} returned multiple compounds."
            raise ValueError(msg)
        else:
            cid = cids[0]

        if cid == 0:
            raise ValueError(f"Structure retrieval for {identifier} failed.")

        try:
            compound = pcp.Compound.from_cid(cid, record_type="3d")

            elements = []
            coordinates = []

            for a in compound.atoms:
                elements.append(a.element)
                coordinates.extend([a.x, a.y, a.z])

            return Structure.from_coordinates(elements, coordinates)
        except pcp.NotFoundError:
            # no 3d structure from pubchem
            # there must be a 2d structure since a cid was found
            [property_dict] = pcp.get_properties(
                properties="IsomericSMILES", identifier=cid, namespace="cid"
            )
            return Structure.generate("smi", property_dict["IsomericSMILES"])

    def rotate(self, axis: np.ndarray, angle: float) -> None:
        """Rotate structure about given axis by given angle.

        Parameters
        ----------
        axis: numpy.ndarray
            Axis of rotation.
            Shape: :math:`(3,1)`
        angle: float
            Angle of rotation (in radians) about the axis of rotation.
        """
        R = rotation_matrix(axis=axis, theta=angle)
        elements = self.atomic_numbers
        coordinates = (self.coords @ R.T).flatten().tolist()
        s = Structure.from_coordinates(elements, coordinates)
        self.pybelmol = s.pybelmol

    def translate(self, v: unyt.unit_array) -> None:
        """Translate structure using given vector.

        Parameters
        ----------
        v: unyt.unit_array
            Translation vector.
            Shape: :math:`(1,3)`
        """
        elements = self.atomic_numbers
        coordinates = (self.coords + v).flatten().tolist()
        s = Structure.from_coordinates(elements, coordinates)
        self.pybelmol = s.pybelmol

    def write(
        self,
        filepath: str,
        filetype: Optional[str] = None,
        overwrite: Optional[bool] = False,
    ) -> None:
        """Write structure to a file.

        Parameters
        ----------
        filepath : str
            Path to file to which the structure will be written.
        filetype : str, optional
            File extension.
            By default `None`.
        overwrite : bool, optional
            If `True`, overwrite `filepath` if it already exists.
            Ignored if `filepath` is a file-like object.
            By default `False`.

        Raises
        ------
        ValueError
            Raised if file extension is not recognized
        """

        filetype = filetype or pathlib.Path(filepath).suffix[1:]

        self.pybelmol.write(filetype, expand(filepath), overwrite=overwrite)

    @contextlib.contextmanager
    def tempfile(
        self, suffix: str, dir: Optional[str] = None
    ) -> Generator[Optional[TextIO], None, None]:
        """Write structure to temporary file.

        Parameters
        ----------
        suffix : str
            Extension for temporary file.
        dir : Optional[str], optional
            Directory containing temporary file.
            By default None.
        """
        with tempfile.NamedTemporaryFile(
            dir=expand(dir) if dir is not None else None, suffix=suffix
        ) as fp:
            try:
                self.write(fp.name, overwrite=True)
                yield fp
            finally:
                pass

    @property
    @functools.cache
    def atomic_masses(self) -> unyt.unyt_array:
        """unyt.unyt_array: Isotopic atomic masses in amu."""
        return [atom.exactmass for atom in self.pybelmol.atoms] * unyt.amu

    @property
    @functools.cache
    def atomic_numbers(self) -> list[int]:
        """list[int]: Atomic numbers."""
        return [atom.atomicnum for atom in self.pybelmol.atoms]

    @property
    def atomic_symbols(self) -> list[str]:
        "list[str]: Atomic symbols."
        return [ob.GetSymbol(atom.atomicnum) for atom in self.pybelmol.atoms]

    @property
    @functools.cache
    def atomic_weights(self) -> unyt.unyt_array:
        """unyt.unyt_array: Elemental atomic masses in amu."""
        return [atom.atomicmass for atom in self.pybelmol.atoms] * unyt.amu

    @property
    @functools.cache
    def centered_coords(self) -> unyt.unyt_array:
        """unyt.unyt_array: Atomic coordinates in COM frame (in angstroms)."""
        return self.coords - self.center_of_mass

    @property
    @functools.cache
    def center_of_mass(self) -> unyt.unyt_array:
        """unyt.unyt_array: Molecular center of mass in angstroms."""
        return ((self.atomic_masses * self.coords.T).T).sum(0) / self.mass

    @property
    def charge(self) -> int:
        """int: Molecular charge."""
        return self.pybelmol.charge

    @charge.setter
    def charge(self, charge: int) -> None:
        return self.obmol.SetTotalCharge(charge)

    @property
    def connectivity(self) -> list[list[int, int, int]]:
        """list[list[int, int, int]]: Bond data (atom index pairs and bond order)."""
        return [
            [bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1, bond.GetBondOrder()]
            for bond in ob.OBMolBondIter(self.obmol)
        ]

    @property
    @functools.cache
    def coords(self) -> unyt.unyt_array:
        """unyt.unyt_array: Atomic coordinates in angstroms."""
        return np.vstack([atom.coords for atom in self.pybelmol.atoms]) * unyt.angstrom

    @property
    @functools.cache
    def distance_matrix(self) -> unyt.unyt_array:
        """unyt.unyt_array: Matrix of squared interatomic distances (in angstrom^2)."""
        # NOTE: equation taken from https://arxiv.org/pdf/1804.04310.pdf
        pp = self.coords @ self.coords.T

        pp_repeat = np.tile(np.diag(pp.value), (self.num_atoms, 1)) * pp.units

        return pp_repeat + pp_repeat.T - 2 * pp

    def fragment(
        self,
        indices: Iterable[Iterable[int]],
        charges: Optional[Iterable[int]] = None,
        multiplicities: Optional[Iterable[int]] = None,
    ) -> list[Structure]:
        """[summary]

        Parameters
        ----------
        indices : Iterable[Iterable[int]]
            [description]
        charges : Optional[Iterable[int]]
            [description]
        multiplicities : Optional[Iterable[int]]
            [description]

        Returns
        -------
        list[Structure]
            [description]
        """
        fragments = []
        charges = charges or []
        multiplicities = multiplicities or []
        for inds, charge, mult in itertools.zip_longest(
            indices, charges, multiplicities
        ):
            elements = np.take(self.atomic_numbers, inds)
            coordinates = self.coords[inds, :].flatten().tolist()
            s = Structure.from_coordinates(elements, coordinates)
            if charge is not None:
                s.charge = charge

            if mult is not None:
                s.multiplicity = mult
            fragments.append(s)

        return fragments

    @property
    def formula(self) -> str:
        """str: Molecular formula."""
        return self.pybelmol.formula

    @property
    @functools.cache
    def inertia_tensor(self) -> unyt.unyt_array:
        """unyt.unyt_array: Moment of inertia tensor (in amu * angstrom^2)."""
        ms = self.atomic_masses
        rs = self.centered_coords
        A = np.einsum("ij, kl -> ikl", rs ** 2, np.eye(3))
        B = np.einsum("...i, ...j -> ...ij", rs, rs)
        return np.einsum("i, ijk -> jk", ms, A - B) * ms.units * rs.units ** 2

    @property
    def isotopes(self) -> list[int]:
        """list[int]: Atomic isotopes."""
        return [atom.isotope for atom in self.pybelmol.atoms]

    @property
    @functools.cache
    def mass(self) -> unyt.unyt_quantity:
        """unyt.unyt_quantity: Molecular mass in amu. """
        return self.pybelmol.exactmass * unyt.amu

    @property
    @functools.cache
    def molar_mass(self) -> unyt.unyt_quantity:
        """unyt.unyt_quantity: Molar mass in amu."""
        return self.pybelmol.molwt * unyt.amu

    @property
    def multiplicity(self) -> int:
        """int: Molecular spin multiplicity."""
        return self.pybelmol.spin

    @multiplicity.setter
    def multiplicity(self, multiplicity: int) -> None:
        return self.obmol.SetTotalSpinMultiplicity(multiplicity)

    @property
    def name(self) -> str:
        """str: Structure name."""
        return self.pybelmol.title

    @name.setter
    def name(self, s: str) -> None:
        self.pybelmol.title = s

    @property
    def num_atoms(self) -> int:
        """int: Number of atoms in structure."""
        return len(self.pybelmol.atoms)

    @property
    def obmol(self) -> ob.OBMol:
        """openbabel.OBMol: Openbabel OBMol representation."""
        return self.pybelmol.OBMol

    @property
    @functools.cache
    def principal_axes(self) -> np.ndarray:
        """numpy.ndarray: Principal axes of moment of inertia tensor."""
        _, axes = scipy.linalg.eigh(self.inertia_tensor.value)
        return axes

    @property
    @functools.cache
    def principal_moments(self) -> unyt.unyt_array:
        """unyt.unyt_array: Principal moments of inertia (in amu * angstrom^2)."""
        return (
            scipy.linalg.eigvalsh(self.inertia_tensor.value) * self.inertia_tensor.units
        )
