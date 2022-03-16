import os
import unittest
import chem
from e3fp.fingerprint.fprint import Fingerprint

glucose = open(os.path.join(os.path.dirname(
    __file__), 'chemicals', 'glucose.mol'), 'r').read()

fluoxetine = open(os.path.join(os.path.dirname(
    __file__), 'chemicals', 'fluoxetine.mol'), 'r').read()

benzene = open(os.path.join(os.path.dirname(
    __file__), 'chemicals', 'benzene.mol'), 'r').read()


class ChemInfo(unittest.TestCase):
    def test_transform_pubchem(self):
        compound_info = {
            "id": {
                "id": {
                    "cid": 241
                }
            },
            "props": [
                {
                    "urn": {
                        "label": "IUPAC Name",
                        "name": "Preferred",
                        "datatype": 1,
                        "version": "2.7.0",
                        "software": "Lexichem TK",
                        "source": "OpenEye Scientific Software",
                        "release": "2021.05.07"
                    },
                    "value": {
                        "sval": "benzene"
                    }
                },
                {
                    "urn": {
                        "label": "Log P",
                        "name": "XLogP3",
                        "datatype": 7,
                        "version": "3.0",
                        "source": "sioc-ccbg.ac.cn",
                        "release": "2021.05.07"
                    },
                    "value": {
                        "fval": 2.1
                    }
                },
            ],
            "count": {
                "heavy_atom": 6,
                "atom_chiral": 0,
                "atom_chiral_def": 0,
                "atom_chiral_undef": 0,
                "bond_chiral": 0,
                "bond_chiral_def": 0,
                "bond_chiral_undef": 0,
                "isotope_atom": 0,
                "covalent_unit": 1,
                "tautomers": -1
            }
        }
        props = chem.transform_pubchem_props(compound_info)
        expected = {
            "Preferred IUPAC Name": 'benzene',
            "XLogP3 Log P": 2.1
        }
        self.assertEqual(props, expected)

    def test_extract_value(self):
        string_value = {
            "ReferenceNumber": 20,
            "Reference": [
                "HANSCH,C ET AL. (1995)"
            ],
            "Value": {
                "StringWithMarkup": [
                    {
                        "String": "log Kow = -2.13 (LogP)"
                    }
                ]
            }
        }

        number_value = {
            "ReferenceNumber": 32,
            "Reference": [
                "Hansch, C., Leo, A., D. Hoekman. Exploring QSAR - Hydrophobic, Electronic, and Steric Constants. Washington, DC: American Chemical Society., 1995., p. 18"
            ],
            "Value": {
                "Number": [-2.13]
            }
        }

        self.assertEqual(chem.extract_value(string_value), -2.13)
        self.assertEqual(chem.extract_value(number_value), -2.13)

    def test_chem_info(self):
        """Should find chemical info"""
        glucose_info = chem.get_chemical_info(
            chem.chem_3d_from_mol_block(glucose))
        self.assertTrue(glucose_info['logP'] < 2.0)

        benzene_info = chem.get_chemical_info(
            chem.chem_3d_from_mol_block(benzene))
        self.assertTrue(2.0 < benzene_info['logP'] < 2.5)

        fluoxetine_info = chem.get_chemical_info(
            chem.chem_3d_from_mol_block(fluoxetine))
        self.assertTrue(4.4 < fluoxetine_info['logP'] < 5.5)


if __name__ == '__main__':
    unittest.main()
