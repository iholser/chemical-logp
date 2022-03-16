from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
from rdkit.Chem import QED
from rdkit.Chem import rdmolops
from rdkit.Chem import rdMolDescriptors
from e3fp.fingerprint.fprint import Fingerprint
from scipy.sparse import csr_matrix
import pickle
import statistics
import urllib.request
import urllib.parse
import json
import re

model = pickle.load(open("svr.p", "rb"))


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
    return mol


def get_fingerprint(mol, logP=""):
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


def mol_from_smiles(smiles):
    return Chem.MolToMolBlock(chem_3d_from_smiles(smiles))


def get_logP(mol):
    fp = Fingerprint.from_rdkit(Chem.GetMorganFingerprintAsBitVect(mol, 2))
    return model.predict(fp.fold(1024).to_vector(sparse=True, dtype=Fingerprint.vector_dtype))[0]


def transform_pubchem_props(compound):
    results = {}
    props = compound['props']
    for p in props:
        if ('name' in p['urn']):
            label = '%s %s' % (p['urn']['name'], p['urn']['label'])
            if ('fval' in p['value']):
                results[label] = float(p['value']['fval'])
            elif ('sval' in p['value']):
                results[label] = str(p['value']['sval'])
            elif ('ival' in p['value']):
                results[label] = int(p['value']['ival'])
    return results


def extract_value(info):
    if ('Number' in info['Value']):
        return float(info['Value']['Number'][0])
    elif ('StringWithMarkup' in info['Value']):
        match = re.search(
            '(-?\d+\.?\d{0,4})', info['Value']['StringWithMarkup'][0]['String'])
        if (match):
            return float(match.group())
    return None


def get_experimental_logp(cid):
    try:
        response = urllib.request.urlopen(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON/?heading=LogP" % str(cid)).read()
        data = json.loads(response.decode('utf-8'))
        values = map(extract_value, data['Record']['Section']
                     [0]['Section'][0]['Section'][0]['Information'])
        return round(statistics.mean(filter(lambda v: v is not None, values)), 2)
    except Exception as e:
        print('%s (%s)' % (str(cid), e))
        pass
    return ''


def get_pubchem_values(mol):
    try:
        response = urllib.request.urlopen(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/json" % urllib.parse.quote(Chem.MolToSmiles(mol, isomericSmiles=True), safe="")).read()
        data = json.loads(response.decode('utf-8'))
        cid = data['PC_Compounds'][0]['id']['id']['cid']
        props = transform_pubchem_props(data['PC_Compounds'][0])
        experimental_logp = get_experimental_logp(cid)
        return {"cid": cid, "experimental_logp": experimental_logp, **props}
    except Exception as e:
        print(e)
        pass
    return None


def get_chemical_info(mol):
    pubchem_result = get_pubchem_values(mol)
    XLogP3 = pubchem_result['XLogP3 Log P'] if 'XLogP3 Log P' in pubchem_result else None

    'XLogP3-AA Log P'
    print(pubchem_result)
    return {
        "mol": Chem.MolToMolBlock(mol),
        "logP": get_logP(mol),
        "XLogP3": XLogP3,
        "experimental_logp": pubchem_result['experimental_logp'] if 'experimental_logp' in pubchem_result else None,
        "wildman_crippen_logP": Descriptors.MolLogP(mol),
        'mw': Descriptors.MolWt(mol),
        'h_donor': Lipinski.NumHDonors(mol),
        'h_acceptor': Lipinski.NumHAcceptors(mol),
        'heavy_atom_count': Lipinski.HeavyAtomCount(mol),
        'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        'molar_refractivity': Crippen.MolMR(mol),
        'topological_surface_area': QED.properties(mol).PSA,
        'formal_charge': rdmolops.GetFormalCharge(mol),
        'ring_count': rdMolDescriptors.CalcNumRings(mol),
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
