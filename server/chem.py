from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from e3fp.fingerprint.fprint import Fingerprint
from scipy.sparse import csr_matrix
import pickle

model = pickle.load(open("clf.p", "rb"))


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


def get_wildman_crippen_logp(mol):
    return Descriptors.MolLogP(mol)


def chem_3d_from_smiles(smiles):
    return chem_3d(Chem.MolFromSmiles(smiles))


def chem_3d_from_mol_block(mol_block):
    return chem_3d(Chem.MolFromMolBlock(mol_block))


def get_logP(mol):
    fp = Fingerprint.from_rdkit(Chem.GetMorganFingerprintAsBitVect(mol, 2))
    arr = csr_matrix(fp.to_vector(sparse=True, dtype=Fingerprint.vector_dtype))
    fold_arr = csr_matrix(
        (arr.data, arr.indices % 1024, arr.indptr),
        shape=arr.shape,
    )
    fold_arr.sum_duplicates()
    fold_arr = fold_arr[:, :1024].tocsr()

    return model.predict(fold_arr)[0]


def get_chemical_info(mol):

    molar_refractivity = Chem.Crippen.MolMR(molecule)
    topological_surface_area_mapping = Chem.QED.properties(molecule).PSA
    formal_charge = Chem.rdmolops.GetFormalCharge(molecule)
    heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(molecule)
    num_of_rings = Chem.rdMolDescriptors.CalcNumRings(molecule)
    return {
        "mol": Chem.MolToMolBlock(mol),
        "logP": get_logP(mol),

        "wildman_crippen_logP": get_wildman_crippen_logp(mol),
        'mw': Descriptors.MolWt(mol),
        'h_donor': Lipinski.NumHDonors(mol),
        'h_acceptor': Lipinski.NumHAcceptors(mol),
        'heavy_atom_count': Lipinski.HeavyAtomCount(mol),
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        'molar_refractivity': Chem.Crippen.MolMR(mol),
        'topological_surface_area': Chem.QED.properties(mol).PSA,
        'formal_charge': Chem.rdmolops.GetFormalCharge(molecule),
        'ring_count': Chem.rdMolDescriptors.CalcNumRings(molecule),
    }


# Lipinski:
#     Moleculer Weight <= 500
#     LogP <= 5
#     H-Bond Donor Count <= 5
#     H-Bond Acceptor Count <= 10
# Ghose:
#     Molecular weight between 180 and 480
#     LogP between -0.4 and +5.6
#     Atom count between 20 and 70
#     Molar refractivity between 40 and 130
# Veber:
#     Rotatable bonds <= 10
#     Topological polar surface area <= 140
# REOS:
#     Molecular weight between 200 and 500
#     LogP between -5.0 and +5.0
#     H-bond donor count between 0 and 5
#     H-bond acceptor count between 0 and 10
#     Formal charge between -2 and +2
#     Rotatable bond count between 0 and 8
#     Heavy atom count between 15 and 50
# Rule of 3:
#     Molecular weight <= 300
#     LogP <= 3
#     H-bond donor <= 3
#     H-bond acceptor count <= 3
#     Rotatable bond count <= 3
# Quantitative Estimate of Druglikeness:
#     mass < 400
#     ring count > 0
#     rotatable bond count < 5
#     h-bond donor count <= 5
#     h-bond acceptor count <= 10
#     logP < 5

# Lipinski CA, Lombardo F, Dominy BW, Feeney PJ (March 2001). "Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings". Adv. Drug Deliv. Rev. 46 (1–3): 3–26. doi:10.1016/S0169-409X(00)00129-0. PMID 11259830.
# Ghose AK, Viswanadhan VN, Wendoloski JJ (January 1999). "A knowledge-based approach in designing combinatorial or medicinal chemistry libraries for drug discovery. 1. A qualitative and quantitative characterization of known drug databases". J Comb Chem. 1 (1): 55–68. doi:10.1021/cc9800071. PMID 10746014.
# Veber DF, Johnson SR, Cheng HY, Smith BR, Ward KW, Kopple KD (June 2002). "Molecular properties that influence the oral bioavailability of drug candidates". J. Med. Chem. 45 (12): 2615–23. CiteSeerX 10.1.1.606.5270. doi:10.1021/jm020017n. PMID 12036371.
# Walters, W. Patrick, Matthew T. Stahl, and Mark A. Murcko. "Virtual screening—an overview." Drug discovery today 3.4 (1998): 160-178.
# Congreve, Miles, et al. "A'rule of three'for fragment-based lead discovery?." Drug discovery today 8.19 (2003): 876-877.
# Bickerton, G. Richard, et al. "Quantifying the chemical beauty of drugs." Nature chemistry 4.2 (2012): 90-98.
