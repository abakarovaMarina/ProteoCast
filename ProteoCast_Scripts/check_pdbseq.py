import argparse
import sys
from Bio import SeqIO
import gemmi as gm

amino_acids = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

def extract_fasta_sequence(fasta_file):
    """Extracts the sequence from a FASTA file."""
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            return str(record.seq).upper()
    return None


def extract_pdb_sequence(pdb_file):
    """Extracts the sequence from a PDB file using SEQRES records."""
    model = gm.read_pdb(pdb_file)[0][0]
    sequence=''.join(amino_acids[residue.name] for residue in model)

    return sequence

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Compare sequences in a FASTA file and a PDB file.")
    parser.add_argument("--fasta_file",  required=True, help="Path to the FASTA file containing the sequence.")
    parser.add_argument("--pdb_file",  required=True, help="Path to the PDB file containing the sequence.")
    args = parser.parse_args()

    if not args.fasta_file:
        print("Error: The fasta file argument is missing.")
        parser.print_help()
    elif not args.pdb_file:
        print("Error: The pdb file argument is missing.")
        parser.print_help()
    else:
        # Extract sequences
        fasta_sequence = extract_fasta_sequence(args.fasta_file)
        pdb_sequence = extract_pdb_sequence(args.pdb_file)
        if fasta_sequence is None:
            print(f"Error: Unable to read the FASTA file: {args.fasta_file}")
            sys.exit(1)

        # Compare sequences
        if fasta_sequence == pdb_sequence:
            print("The sequences in the FASTA file and PDB file are identical.")
            sys.exit(0) 
        else:
            print("The sequences in the FASTA file and PDB file are different.")
            sys.exit(1) 

if __name__ == "__main__":
    main()

