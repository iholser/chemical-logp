from rdkit.Chem import AllChem as Chem

# smiles -> mol -> 3d -> molBlock, 3d fingerprint, daylight fingerprint, macs, inchi, inchi key
# addHs


def chem_info(mol):
    """take in smiles, return mol_block (3d), fp_on_bits (daylight fingerprint) """
    try:
        mol.UpdatePropertyCache()
        mol = Chem.AddHs(mol)
        Chem.EmbedMolecule(mol, randomSeed=0xf00d)
        mol_block = Chem.MolToMolBlock(mol)
        fp_on_bits = list(Chem.RDKFingerprint(mol).GetOnBits())
        return mol_block, fp_on_bits, mol
    except Exception as e:
        pass


def chem_info_from_smiles(smiles):
    return chem_info(Chem.MolFromSmiles(smiles))
