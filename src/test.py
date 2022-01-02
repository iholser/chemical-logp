import os
import unittest
import chem

glucose = open(os.path.join(os.path.dirname(
    __file__), 'glucose.mol'), 'r').read()


class ChemInfo(unittest.TestCase):
    def test_chem_info(self):
        """Should find chemical info"""
        mol_block, fp_on_bits, mol = chem.chem_info_from_smiles(
            'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O')
        self.assertEqual(mol_block, glucose)


if __name__ == '__main__':
    unittest.main()
