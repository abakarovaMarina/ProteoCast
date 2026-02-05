import sys
from Bio.Align import PairwiseAligner
from Bio import SeqIO
import gemmi as gm
import pandas as pd
## Pairwise Alignment

GAP_CHAR = "-"
amino_acids = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
##read pdb sequence, only of the chain of interest 
def PDBtoSeq(pdb_path, chain_id='A'):
    model = gm.read_pdb(pdb_path)[0]
    chain = model[chain_id]
    sequence = ''.join(amino_acids[residue.name] for residue in chain if residue.name in amino_acids)
    for residue in chain:
        if residue.name in amino_acids:
            first_chain_residue_nb = residue.seqid.num
            break
    return sequence, first_chain_residue_nb

#structure
pdb_seq, first_res_nb = PDBtoSeq(sys.argv[1], sys.argv[2])  

#MSA
sequences = list(SeqIO.parse(sys.argv[3], "fasta"))
fasta_seq = sequences[0].seq


#align
aligner = PairwiseAligner()
aligner.mode = 'global'

aligner.match_score = 1.0
aligner.mismatch_score = -3.0
aligner.open_gap_score = -2.5
aligner.extend_gap_score = -2.0
aligner.query_end_gap_score = -2.0
aligner.target_end_gap_score = -2.0

alignment = aligner.align(fasta_seq, pdb_seq)[0]

def mapping_identity_from_aligned(alignment, seq1, seq2):

    seq1 = str(seq1)
    seq2 = str(seq2)

    blocks1, blocks2 = alignment.aligned
    mapping = {}
    matches = 0

    for (s1, e1), (s2, e2) in zip(blocks1, blocks2):
        block_len = min(e1 - s1, e2 - s2)
        if block_len <= 0:
            continue

        for k in range(block_len):
            i1 = s1 + k
            i2 = s2 + k
            mapping[i1 + 1] = i2 + 1  # 1-based

            if seq1[i1] == seq2[i2]:
                matches += 1

    if aligner.mode == "global":
        denom = len(seq1)
    elif aligner.mode == "local":
        denom = min(len(seq1), len(seq2))
    else:
        raise ValueError("mode must be 'global' or 'local'")
    
    identity = matches / denom if denom > 0 else 0.0
    return mapping, identity


mapping, identity = mapping_identity_from_aligned(alignment, fasta_seq, pdb_seq)

threshold = 0.20
with open("/opt/job/identity.txt", "w") as f:    
    if identity < threshold:
        f.write(f"True")
    else:
        f.write(f"False")

df = pd.DataFrame(list(mapping.items()), columns=['fasta_pos', 'pdb_pos'])
df['pdb_pos'] = df['pdb_pos'] + first_res_nb - 1
df.to_csv('/opt/job/mapping_fasta_pdb.csv', index=False)
