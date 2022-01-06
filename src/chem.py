from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors

from e3fp.fingerprint.fprint import Fingerprint

# smiles -> mol -> 3d -> molBlock, 3d fingerprint, daylight fingerprint, macs, inchi, inchi key
# addHs


# def chem_info(mol):
#     """take in smiles, return mol_block (3d), fp_on_bits (daylight fingerprint) """
#     try:
#         mol.UpdatePropertyCache()
#         mol = Chem.AddHs(mol)
#         Chem.EmbedMolecule(mol, randomSeed=0xf00d)
#         mol_block = Chem.MolToMolBlock(mol)
#         # fp_on_bits = list(Chem.RDKFingerprint(mol).GetOnBits())
#         # fp = Fingerprint.from_rdkit(
#         #     Chem.GetMorganFingerprintAsBitVect(mol, 2), name='fp1')  # TODO: use name property for fingerprint database
#         return mol_block, mol
#     except Exception as e:
#         print(e)
#         pass
#     return None, None, None


def chem_3d(mol):
    try:
        mol.UpdatePropertyCache()
        mol = Chem.AddHs(mol)
        # seed is necessary to ensure that 3d coordinates are the same for testing
        Chem.EmbedMolecule(mol, randomSeed=0xf00d)
        return mol
    except Exception as e:
        print(e)
        pass
    return None


def get_fingerprint(mol, logP):
    # regenerate smiles to ensure that they are cononical
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    return Fingerprint.from_rdkit(
        Chem.GetMorganFingerprintAsBitVect(mol, 2), name=smiles, props={'logP': logP})


def get_rdkit_logp(mol):
    return Descriptors.MolLogP(mol)


def chem_3d_from_smiles(smiles):
    return chem_3d(Chem.MolFromSmiles(smiles))
