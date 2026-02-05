from typing import Dict
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.SASA import ShrakeRupley
from amino_acid import AminoAcid
import os
import sys
import pandas as pd

amino_acids = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

MAX_SURFACE_MAP = {
        "ALA": 1.29,
        "ARG": 2.74,
        "ASN": 1.95,
        "ASP": 1.93,
        "CYS": 1.67,
        "GLN": 2.23,
        "GLU": 2.25,
        "GLY": 1.04,
        "HIS": 2.24,
        "ILE": 1.97,
        "LEU": 2.01,
        "LYS": 2.36,
        "MET": 2.24,
        "PHE": 2.40,
        "PRO": 1.59,
        "SER": 1.55,
        "THR": 1.55,
        "TRP": 2.85,
        "TYR": 2.63,
        "VAL": 1.74,
}
MAX_SURFACE_DEFAULT = 2.01 # mean value

def get_max_surf(aa_three):
    return MAX_SURFACE_MAP.get(aa_three, MAX_SURFACE_DEFAULT)


pdb_path = sys.argv[1]
chain_input = sys.argv[2]
#pdb_path = '../1a4y.pdb'

# Parse PDB file
pdb_name = os.path.basename(pdb_path).replace(".pdb", "")
pdb_parser = PDBParser(QUIET=True)
structure = pdb_parser.get_structure(pdb_name, pdb_path)

# Compute ASA
shrake_rupley = ShrakeRupley(
    #probe_radius=1.40, # radius of the probe in A. Default is 1.40, roughly the radius of a water molecule.
    #n_points=200,      # resolution of the surface of each atom. Default is 100. A higher number of points results in more precise measurements, but slows down the calculation.
    #radii_dict=None,   # user-provided dictionary of atomic radii to use in the calculation. Values will replace/complement those in the default ATOMIC_RADII dictionary.
)
shrake_rupley.compute(structure, level="R")


nb_chain = 0


rsa_map = []
for chain_obj in structure[0]:
    nb_chain += 1
    chain = chain_obj.id
    chain_structure = structure[0][chain]
    for residue in chain_structure:
        aa_three_raw = residue.resname
        aa_three = AminoAcid._NON_STANDARD_AAS.get(aa_three_raw, aa_three_raw)

        hetero_flag, resseq, icode = residue.id

        # keep standard residues or mapped residues only
        if hetero_flag != ' ' and aa_three not in amino_acids:
            continue

        aa_one = amino_acids.get(aa_three, 'X')

        asa = residue.sasa
        if isinstance(asa, float):
            rsa_map.append({
                "chain": chain,
                "residue": int(resseq),
                "icode": icode.strip() if isinstance(icode, str) else "",
                "aa": aa_one,
                "aa3": aa_three,
                "rsa": asa / get_max_surf(aa_three),
            })

df = pd.DataFrame(rsa_map)
df["resid"] = df["chain"] + df["aa"] + df["residue"].astype(str) + df["icode"]


if nb_chain > 1:
    df['rsa_monomer'] = None
    ###unbound rsa calculation for the chain of interest
    io = PDBIO()
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            if chain_input == chain_id:
                io.set_structure(chain)
                io.save(f"{pdb_path[:-4]}_{chain_input}.pdb")

    pdb_path = f"{pdb_path[:-4]}_{chain_input}.pdb" 
    # Parse PDB file
    pdb_name = os.path.basename(pdb_path).replace(".pdb", "")
    pdb_parser = PDBParser(QUIET=True)
    structure_monomer = pdb_parser.get_structure(pdb_name, pdb_path)
    shrake_rupley.compute(structure_monomer, level="R")


    # Convert to RSA and format
    rsa_map_monomer = []
    for chain_obj in structure_monomer[0]:
        chain = chain_obj.id
        chain_structure = structure_monomer[0][chain]
        for residue in chain_structure:
            aa_three_raw = residue.resname
            aa_three = AminoAcid._NON_STANDARD_AAS.get(aa_three_raw, aa_three_raw)

            hetero_flag, resseq, icode = residue.id

            # keep standard residues or mapped residues only
            if hetero_flag != ' ' and aa_three not in amino_acids:
                continue

            aa_one = amino_acids.get(aa_three, 'X')

            asa = residue.sasa
            if isinstance(asa, float):
                rsa_val = asa / get_max_surf(aa_three)
                resid = chain + aa_one + str(int(resseq)) + (icode.strip() if isinstance(icode, str) else "")
                df.loc[(df['chain'] == chain_input) & (df['resid'] == resid), 'rsa_monomer'] = rsa_val
        
    try:
       os.remove(pdb_path)
    except OSError as e:
       print(f"Warning: could not delete temporary file {pdb_path}: {e}")

df.to_csv('/opt/job/rsa_values.csv',index = False)

