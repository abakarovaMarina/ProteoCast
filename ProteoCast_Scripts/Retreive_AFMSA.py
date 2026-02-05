import argparse
import os
import requests

def download_alphafold_structure(uniprot_id, output_dir):
    """
    Download the a3m for a given UniProt ID.

    Args:
        uniprot_id (str): The UniProt ID of the protein.
        output_dir (str): Directory to save the downloaded structure.
    """
    base_url = "https://alphafold.ebi.ac.uk/files/msa"
    pdb_filename = f"AF-{uniprot_id}-F1-msa_v6.a3m"
    download_url = f"{base_url}/{pdb_filename}"
    
    try:
        response = requests.get(download_url)
        response.raise_for_status()  # Raise an error for HTTP issues
        
        # Save the file to the specified directory
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, pdb_filename)
        with open(output_path, "wb") as pdb_file:
            pdb_file.write(response.content)
        
        print(f"Structure downloaded successfully: {output_path}")
    except requests.HTTPError as e:
        print(f"Failed to download structure for UniProt ID '{uniprot_id}': {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Download a3m using UniProt ID.")
    parser.add_argument("uniprot_id", type=str, help="UniProt ID of the protein.")
    parser.add_argument("-o", "--output_dir", type=str, default=".", help="Directory to save the downloaded a3m.")
    
    args = parser.parse_args()
    
    # Download the AlphaFold structure
    download_alphafold_structure(args.uniprot_id, args.output_dir)

