"""
Author  - Fernando.S.D.H - 15817
Date - 02/03/2025
Name of the project - Protein Analyzer

"""

import os
import time
from Bio import ExPASy, SwissProt, SeqIO
from Bio.ExPASy import ScanProsite
from Bio.PDB import PDBList, PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class ProteinAnalyzer:

    def __init__(self, output_dir):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    """
    Method 01
    Fetches SwissProt sequences for given UniProt IDs and saves each in a separate FASTA file.
    input_file : Path to the text file containing SwissProt IDs (one per line).
    output_dir : Directory where output FASTA files will be saved.

    1. IDs are given in a txt file taken from swissprot database as a list
    2. for each id the sequences are taken in fasta format
    3. save the sequences to an output directory

    """
    # Method 1: Fetch SwissProt Sequences and save as FASTA
    def fetch_swissprot_records(self,input_file):
        fasta_output_dir = os.path.join(self.output_dir, "fasta_sequences")

        # Ensure output directory exists
        os.makedirs(fasta_output_dir, exist_ok=True)

        # Read UniProt IDs from file
        with open(input_file, "r") as text_file:
            ids = [line.strip() for line in text_file if line.strip()]

        for uniprot_id in ids:
            # Fetch SwissProt record
            handle = ExPASy.get_sprot_raw(uniprot_id)
            record = SwissProt.read(handle)
            handle.close()

            # Prepare FASTA content - fasta format
            fasta_content = f">{record.entry_name}|{record.accessions[0]}\n{record.sequence}\n"

            # Save to a separate FASTA file
            fasta_filename = os.path.join(fasta_output_dir, f"{uniprot_id}.fasta")
            with open(fasta_filename, "w") as fasta_file:
                fasta_file.write(fasta_content)

            print(f"Saved: {fasta_filename}")
            time.sleep(1)  # To avoid overwhelming UniProt servers

        print(f"All sequences saved in {fasta_output_dir}")

    """
    Methpd 02
    input : set of protein records in fasta format
    output : text file that prints domain coordinates, domain ids , domain names found for each sequence

    steps:
    1. open the output file and write the headers to it
    2. read the fasta sequences by open the input file from the directory
    3. Scan the expasy prosite database
    4. from the results list which for one sequence multiple domain hits present, print the domain id, name and coordinates
    5.To output file write the extracted results

    """
    # Method 2: Scan FASTA Sequences for Domains using Prosite
    def scan_prosite_from_fasta(self):
        input_dir = os.path.join(self.output_dir, "fasta_sequences")
        output_file = os.path.join(self.output_dir, "domain_scan_results.txt")

        with open(output_file, "w") as output:
            # Write header
            output.write("Sequence_ID\t\t\tDomain_ID\t\tDomain_Name\t\tStart\t\tEnd\n")

            for fasta_file in os.listdir(input_dir):
                if fasta_file.endswith(".fasta"):
                    file_path = os.path.join(input_dir, fasta_file)

                    for record in SeqIO.parse(file_path, "fasta"):
                        sequence = record.seq

                        # Scan the sequence using ScanProsite
                        handle = ScanProsite.scan(sequence)
                        results = ScanProsite.read(handle)

                        # print(f"ScanProsite results:\n{results[:1000]}")
                    if results:
                        for result in results:
                            # Extracting domain information
                            domain_id = result.get('signature_ac', '')  # Domain ID (signature accession)
                            domain_name = result.get('signature_ac', '')  # Use an empty string if no name is available
                            start = result.get('start', '')  # coordinates start
                            end = result.get('stop', '')  # coordinates end

                            # Write results to the output file
                            output.write(f"{record.id}\t\t{domain_id}\t\t{domain_name}\t\t{start}\t\t\t{end}\n")

    """
    Method03
    Input - pdb ids in  atext file
    Output - downloaded pdb files (structures) and sequences of those pdb files in fasta format

    1.Open and read the input file with pdb ids
    2.Using biopython pdblist() module download the pdb files and iterate through the list to download all
    3.At the same time get the sequence of each structure as a fasta format by using pdbparser()
    4.In each structure need to extract the models, chins and each residues
    5.Then save all the sequences in a single fasta file

    """
    # Method 3: Download PDB Files and Extract Sequences
    def download_pdb_files(self, pdb_id_file):
        pdb_output_dir = os.path.join(self.output_dir, "pdb_files")
        fasta_output = os.path.join(self.output_dir, "pdb_sequences.fasta")

        # Ensure output directory exists
        os.makedirs(pdb_output_dir, exist_ok=True)

        # Read PDB IDs from file
        with open(pdb_id_file, "r") as pdb_input:
            pdb_ids = [line.strip() for line in pdb_input if line.strip()]

        pdbl = PDBList()
        parser = PDBParser(QUIET=True)
        fasta_records = []

        for pdb_id in pdb_ids:
            # pdb_id = pdb_id.lower()
            pdb_path = pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_output_dir, file_format="pdb")

            if pdb_path:
                structure = parser.get_structure(pdb_id, pdb_path)
                for model in structure:
                    for chain in model:
                        seq = "".join([residue.resname for residue in chain if residue.id[0] == " "])
                        seq = Seq(seq.replace("UNK", "X"))  # Replace unknown residues with 'X'
                        fasta_records.append(SeqRecord(seq, id=f"{pdb_id}_{chain.id}", description=""))

        # Save all sequences in a single FASTA file
        with open(fasta_output, "w") as fasta_file:
            SeqIO.write(fasta_records, fasta_file, "fasta")

        print(f"FASTA sequences saved to {fasta_output}")

    """
    Method 04
    input - single pdb file
    output - text file with the extracted details

    1. make the output file header 
    2. read the pdb file
    3. extract the number od chains, residues and number of atoms for each residue in protein 
    4. save the results in one text file

    """
    # Method 4: Extract Details from Single PDB File
    def details_from_pdb_file(self,input_file):
        output_file = os.path.join(self.output_dir, "single_file_detail.txt")

        parser = PDBParser(QUIET=True)

        with open(output_file, "w") as output:
            # make the output file header
            output.write("PDB_ID\t\t\tChain_ID\t\tResidue_ID\t\tResidue_Name\t\tNumber_of_atoms \n")

            pdb_id = os.path.splitext(os.path.basename(input_file))[0]

            # parse the PDB structure
            structure = parser.get_structure(pdb_id, input_file)

            num_chains = 0
            num_residues = 0
            num_atoms_total = 0

            for model in structure:
                for chain in model:
                    num_chains += 1
                    for residue in chain:
                        if residue.id[0] == " ":  # Exclude heteroatoms (HETATM)
                            num_residues += 1
                            residue_id = residue.id[1]  # Residue sequence number
                            residue_name = residue.resname  # Residue three-letter code
                            num_atoms = len(residue)  # Number of atoms in this residue
                            output.write(
                                f"{pdb_id}\t\t\t\t {chain.id}\t\t\t {residue_id}\t\t\t\t{residue_name}\t\t\t\t{num_atoms}\n")

                            num_atoms_total += num_atoms

            # write to the output file
            output.write(f"\nTotal Chains: {num_chains}\n")
            output.write(f"Total Residues: {num_residues}\n")
            output.write(f"Total atoms: {num_atoms_total}\n")

        print(f"Details saved in {output_file}")

    """
    Method 05
    Input- The pdb files 
    Output- number of chains , residues and atoms for each residue in separate text files

    1. Read the pdb files one by one in the directory
    2.initialize the variables to 0
    2. extract the details and calculate from the files 
    3.save each output in a separate text files

    """
    # Method 5: Extract Details from Multiple PDB Files
    def details_from_multiple_pdb(self):
        input_dir = os.path.join(self.output_dir, "pdb_files")

        parser = PDBParser(QUIET=True)

        for pdb_file in os.listdir(input_dir):
            if pdb_file.endswith(".ent"):
                pdb_id = pdb_file.split('.')[0]  # Extracting PDB ID
                file_path = os.path.join(input_dir, pdb_file)
                output_file = os.path.join(input_dir, f"{pdb_id}_details.txt")  # Output file for each PDB

                # parse the PDB structure
                structure = parser.get_structure(pdb_id, file_path)

                with open(output_file, "w") as output:
                    # make the output file header
                    output.write("PDB_ID\t\t\tChain_ID\t\tResidue_ID\t\tResidue_Name\t\tNumber_of_atoms \n")

                    num_chains = 0
                    num_residues = 0
                    num_atoms_total = 0

                    for model in structure:
                        for chain in model:
                            num_chains += 1
                            for residue in chain:
                                if residue.id[0] == " ":  # Exclude heteroatoms (HETATM)
                                    num_residues += 1
                                    residue_id = residue.id[1]  # Residue sequence number
                                    residue_name = residue.resname  # Residue three-letter code
                                    num_atoms = len(residue)  # Number of atoms in this residue
                                    output.write(
                                        f"{pdb_id}\t\t\t\t {chain.id}\t\t\t {residue_id}\t\t\t\t{residue_name}\t\t\t\t{num_atoms}\n")

                                    num_atoms_total += num_atoms

                        # write to the output file
                    output.write(f"\nTotal Chains: {num_chains}\n")
                    output.write(f"Total Residues: {num_residues}\n")
                    output.write(f"Total atoms: {num_atoms_total}\n")

                print(f"Details saved in {output_file}")

    """
    Method 06
    input - a pdb file
    output - all the hetero protein residues in text file 

    1.make the output file and write the header
    2.read the given pdb file , parse the structure
    3.loop through all the models, chains , residues to find the hetero atoms 
    Hetero atoms are the ones that residue.id is not equal to " "
    4. save the output to te output file
    """
    # Method 6: Extract Hetero Residues from a PDB File
    def hetero_residues(self,pdb_file):
        output_file = os.path.join(self.output_dir, "hetero_residues.txt")

        parser = PDBParser(QUIET=True)

        with open(output_file, "w") as output:
            output.write(f"Hetero residues in {pdb_file}\n")

            pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
            # parse the structure
            structure = parser.get_structure(pdb_id, pdb_file)

            # loop through the models, chains, residues for find the hetero atoms
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] != " ":  # hetero atoms are the ones that residue.id is not equal to " "(blank)
                            residue_id = residue.id[1]
                            residue_name = residue.resname
                            chain_id = chain.id

                            output.write(
                                f"{pdb_id}\t\t  Chain :{chain_id}\t\t  Residue: {residue_id} \t \t  Residue name {residue_name}\n")

        print(f"Hetero atoms are saved to {output_file}")

    """
    Method 07

    input:pdb file
    Output:text file containing the coordinates for given atom 

    1. Open the output file and write the headers to it
    2. get the pdb id from the input db file
    3. parse the structure from the file
    4. loop through the models, chains , residues and atoms in each residue 
    5. for each atom given in the method find the coordinates
    6. get the residue id, chain id and residue name also
    7. save the output in a text file

    """
    # Method 7: Get Atom Coordinates
    def get_atom_coordinates(self,pdb_file, atom_name):
        output_file = os.path.join(self.output_dir, "Coordinates.txt")

        parser = PDBParser(QUIET=True)

        # open the output file
        with open(output_file, "w") as output:
            output.write(f"Coordinates for atoms\n")

            # get the pdb id
            pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]

            # parse the structure
            structure = parser.get_structure(pdb_id, pdb_file)

            # looping through the models, chains, residues and atoms to get coordinates
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            if atom.get_name() == atom_name:  # match atom by name
                                x, y, z = atom.get_coord()  # get the coordinates
                                residue_id = residue.id[1]  # get residue sequence number
                                residue_name = residue.resname  # get residue name
                                chain_id = chain.id

                                output.write(
                                    f"{pdb_id} \t\t Chain : {chain_id}\t\t Residue_id: {residue_id}\t\t Residue name: {residue_name}\t\t Coordinates{atom_name}\t:{x:.3f}\t\t ,{y:.3f}\t\t ,{z:.3f}\t\n")


# # Main function
# def main():
#     analyzer = ProteinAnalyzer(output_dir="output")
#
#     analyzer.fetch_swissprot_records("swissprotIDs.txt")
#     analyzer.scan_prosite_from_fasta()
#     analyzer.download_pdb_files("pdb_ids.txt")
#     analyzer.details_from_pdb_file("pdb1gjh.ent")
#     analyzer.details_from_multiple_pdb()
#     analyzer.hetero_residues("pdb1jnx.ent")
#     analyzer.get_atom_coordinates("pdb1jnx.ent", "C")
#
# if __name__ == "__main__":
#     main()


def main():
    # Ask for the directory to store output
    output_dir = input("Enter the directory to store output files: ") #default - output
    if not output_dir:
        output_dir = "output"  # Default directory if the user doesn't provide one

    # Create an instance of the ProteinAnalyzer
    analyzer = ProteinAnalyzer(output_dir=output_dir)

    # Ask for SwissProt IDs file and process it
    swissprot_file = input("Enter the path to the SwissProt IDs file (e.g., swissprotIDs.txt): ")
    if os.path.exists(swissprot_file):
        analyzer.fetch_swissprot_records(swissprot_file)
    else:
        print(f"File {swissprot_file} does not exist.")

    # Scan Prosite domains from FASTA files
    print("Scanning Prosite domains from FASTA sequences...")
    analyzer.scan_prosite_from_fasta()

    # Ask for PDB IDs file and process it
    pdb_ids_file = input("Enter the path to the PDB IDs file (e.g., pdb_ids.txt): ")
    if os.path.exists(pdb_ids_file):
        analyzer.download_pdb_files(pdb_ids_file)
    else:
        print(f"File {pdb_ids_file} does not exist.")

    # Ask for a PDB file to extract details from
    pdb_file = input("Enter the path to a single PDB file for detailed extraction (e.g., pdb1gjh.ent): ")
    if os.path.exists(pdb_file):
        analyzer.details_from_pdb_file(pdb_file)
    else:
        print(f"File {pdb_file} does not exist.")

    # Extract details from multiple PDB files in the directory
    print("Extracting details from multiple PDB files in the directory...")
    analyzer.details_from_multiple_pdb()

    # Ask for a specific PDB file to extract hetero residues
    hetero_pdb_file = input("Enter the path to a PDB file to extract hetero residues (e.g., pdb1jnx.ent): ")
    if os.path.exists(hetero_pdb_file):
        analyzer.hetero_residues(hetero_pdb_file)
    else:
        print(f"File {hetero_pdb_file} does not exist.")

    # Ask for a PDB file and atom name to get atom coordinates
    pdb_file_for_coords = input("Enter the path to a PDB file to extract atom coordinates (e.g., pdb1jnx.ent): ")
    atom_name = input("Enter the atom name (e.g., 'C') to get its coordinates: ")
    if os.path.exists(pdb_file_for_coords):
        analyzer.get_atom_coordinates(pdb_file_for_coords, atom_name)
    else:
        print(f"File {pdb_file_for_coords} does not exist.")

if __name__ == "__main__":
    main()
