from materia.structure import Structure
import numpy as np
import pubchempy as pcp
import unittest
import unittest.mock as mock
import unyt


class MockCCData:
    def __init__(self, atomnos, atomcoords):
        self.atomnos = atomnos
        self.atomcoords = atomcoords


class TestStructure(unittest.TestCase):
    def setUp(self):
        self.elements = [1, 1, 8]
        self.water_coords = np.array(
            [[[0.757, 0.586, 0.000], [-0.757, 0.586, 0.000], [0.000, 0.000, 0.000]]]
        )
        self.flat_water_coords = self.water_coords.flatten().tolist()

    def test_from_coordinates(self):
        Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertTrue(True)

    def test_from_cclib(self):
        parsed = mock.MagicMock()
        parsed.atomnos = self.elements
        parsed.atomcoords = self.water_coords

        Structure.from_cclib(parsed)
        self.assertTrue(True)

    def test_from_schema(self):
        schema = mock.MagicMock()
        schema.molecule.symbols = ["H", "H", "O"]
        schema.molecule.geometry = self.flat_water_coords

        Structure.from_schema(schema)
        self.assertTrue(True)

    def test_read_xyz(self):
        with mock.patch("openbabel.pybel.readfile") as m1:
            Structure.read("/home/user/some_file.xyz")
            m1.assert_called_with("xyz", "/home/user/some_file.xyz")

    def test_write_xyz(self):
        pybelmol = mock.MagicMock()
        pybelmol.write = mock.MagicMock()
        s = Structure(pybelmol)
        s.write("/home/user/ok.xyz", overwrite=True)
        pybelmol.write.assert_called_with("xyz", "/home/user/ok.xyz", overwrite=True)
        s.write("/home/user/ok.xyz", overwrite=False)
        pybelmol.write.assert_called_with("xyz", "/home/user/ok.xyz", overwrite=False)

    def test_tempfile_xyz(self):
        pybelmol = mock.MagicMock()
        pybelmol.write = mock.MagicMock()
        s = Structure(pybelmol)
        with s.tempfile("xyz") as f:
            pybelmol.write.assert_called_with("", f.name, overwrite=True)

    def test_retrieve(self):
        a1 = mock.MagicMock()
        a1.element = 6
        a1.x = 0
        a1.y = 0
        a1.z = 0

        a2 = mock.MagicMock()
        a2.element = 1
        a2.x = 0.2774
        a2.y = 0.8929
        a2.z = 0.2544

        a3 = mock.MagicMock()
        a3.element = 1
        a3.x = 0.6068
        a3.y = -0.2383
        a3.z = -0.7169

        compound = mock.MagicMock()
        compound.atoms = [a1, a2, a3]

        with self.assertRaises(ValueError):
            Structure.retrieve("fakeidtype", "benzene")

        with mock.patch("pubchempy.get_cids") as m:
            m.side_effect = OSError
            with self.assertRaises(ValueError):
                Structure.retrieve("smiles", "ap3r8jpf9aj3")

        with mock.patch("pubchempy.get_cids") as m1, mock.patch(
            "pubchempy.Compound.from_cid"
        ) as m2:
            m1.return_value = [0]
            with self.assertRaises(ValueError):
                Structure.retrieve("smiles", "O")

            m1.return_value = [962, 962, 962]
            m2.return_value = compound
            with self.assertRaises(ValueError):
                Structure.retrieve("smiles", "O")

            Structure.retrieve("smiles", "O", use_first_match=True)

            m1.return_value = [962]
            Structure.retrieve("smiles", "O")

            m2.side_effect = pcp.NotFoundError
            with mock.patch("pubchempy.get_properties") as m3:
                m3.return_value = [{"IsomericSMILES": "O"}]
                Structure.retrieve("smiles", "O")

    def test_atomic_symbols(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        self.assertEqual(he.atomic_symbols, ["He"])

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertEqual(water.atomic_symbols, ["H", "H", "O"])

    def test_atomic_weights(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        self.assertTrue(unyt.allclose_units(he.atomic_weights, [4.002602] * unyt.amu))

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertTrue(
            unyt.allclose_units(
                water.atomic_weights, [1.00794, 1.00794, 15.9994] * unyt.amu
            )
        )

    def test_center_of_mass(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        print(he.center_of_mass)
        self.assertTrue(
            unyt.allclose_units(he.center_of_mass, [0.0, 0.0, 0.0] * unyt.angstrom)
        )

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertTrue(
            unyt.allclose_units(
                water.center_of_mass, [0.0, 0.06558212, 0.0] * unyt.angstrom
            )
        )

    def test_connectivity(self):
        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertEqual(water.connectivity, [[0, 2, 1], [1, 2, 1]])

    def test_generate_charge(self):
        s = Structure.generate("smi", "O")
        self.assertEqual(s.charge, 0)

        s = Structure.generate("smi", "c1ccccc1")
        self.assertEqual(s.charge, 0)

        s = Structure.generate("smi", "[NH4+]")
        self.assertEqual(s.charge, 1)

        s = Structure.generate("smi", "[Ti+4]")
        self.assertEqual(s.charge, 4)

        s = Structure.generate("smi", "[OH-]")
        self.assertEqual(s.charge, -1)

    def test_distance_matrix(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        self.assertTrue(
            unyt.allclose_units(he.distance_matrix, [[0.0]] * unyt.angstrom ** 2)
        )

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertTrue(
            unyt.allclose_units(
                water.distance_matrix,
                [
                    [0.0, 2.292196, 0.916445],
                    [2.292196, 0.0, 0.916445],
                    [0.916445, 0.916445, 0.0],
                ]
                * unyt.angstrom ** 2,
            )
        )

    def test_fragment(self):
        s = Structure.from_coordinates(self.elements, self.flat_water_coords)

        indices = [[0, 2], [1]]
        fragments = s.fragment(indices)

        self.assertEqual(len(fragments), 2)
        self.assertEqual(fragments[0].charge, 0)
        self.assertEqual(fragments[0].num_atoms, 2)
        self.assertEqual(fragments[1].charge, 0)
        self.assertEqual(fragments[1].num_atoms, 1)

        charges = [-1, 1]
        multiplicities = [2, 2]
        fragments = s.fragment(indices, charges, multiplicities)

        self.assertEqual(len(fragments), 2)
        self.assertEqual(fragments[0].charge, -1)
        self.assertEqual(fragments[0].num_atoms, 2)
        self.assertEqual(fragments[0].multiplicity, 2)
        self.assertEqual(fragments[1].charge, 1)
        self.assertEqual(fragments[1].num_atoms, 1)
        self.assertEqual(fragments[1].multiplicity, 2)

    def test_isotopes(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        self.assertEqual(he.isotopes, [0])

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertEqual(water.isotopes, [0, 0, 0])

    def test_formula(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        self.assertEqual(he.formula, "He")

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertEqual(water.formula, "H2O")

    def test_inertia_tensor(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        self.assertTrue(
            unyt.allclose_units(
                he.inertia_tensor,
                np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
                * unyt.amu
                * unyt.angstrom ** 2,
            )
        )

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertTrue(
            unyt.allclose_units(
                water.inertia_tensor,
                np.array(
                    [
                        [0.61470248, 0.0, 0.0],
                        [0.0, 1.15506625, 0.0],
                        [0.0, 0.0, 1.76976873],
                    ]
                )
                * unyt.amu
                * unyt.angstrom ** 2,
            )
        )

        water.rotate([0, 0, 1], np.pi / 2)
        self.assertTrue(
            unyt.allclose_units(
                water.inertia_tensor,
                np.array(
                    [
                        [0.61470248, 0.0, 0.0],
                        [0.0, 1.15506625, 0.0],
                        [0.0, 0.0, 1.76976873],
                    ]
                )
                * unyt.amu
                * unyt.angstrom ** 2,
            )
        )

        water.translate([0, 0, 1] * unyt.angstrom)
        self.assertTrue(
            unyt.allclose_units(
                water.inertia_tensor,
                np.array(
                    [
                        [0.61470248, 0.0, 0.0],
                        [0.0, 1.15506625, 0.0],
                        [0.0, 0.0, 1.76976873],
                    ]
                )
                * unyt.amu
                * unyt.angstrom ** 2,
            )
        )

    def test_molar_mass(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        self.assertTrue(unyt.allclose_units(he.molar_mass, 4.002602 * unyt.amu))

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertTrue(unyt.allclose_units(water.molar_mass, 18.01528 * unyt.amu))

    def test_name(self):
        s = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        s.name = "He"
        self.assertEqual(s.name, "He")

    def test_obmol(self):
        pybelmol = mock.MagicMock(OBMol="dummy")

        s = Structure(pybelmol)
        self.assertEqual(s.obmol, "dummy")

    def test_principal_axes(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        self.assertTrue(np.allclose(he.principal_axes, np.eye(3)))

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertTrue(np.allclose(water.principal_axes, np.eye(3)))

    def test_principal_moments(self):
        he = Structure.from_coordinates(["He"], [0.000, 0.000, 0.000])
        self.assertTrue(
            unyt.allclose_units(
                he.principal_moments, [0.0, 0.0, 0.0] * unyt.amu * unyt.angstrom ** 2
            )
        )

        water = Structure.from_coordinates(self.elements, self.flat_water_coords)
        self.assertTrue(
            unyt.allclose_units(
                water.principal_moments,
                [0.61470248, 1.15506625, 1.76976873] * unyt.amu * unyt.angstrom ** 2,
            )
        )
