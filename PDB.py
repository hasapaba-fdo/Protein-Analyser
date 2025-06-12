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

import os
from Bio.PDB import PDBList, PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def download_pdb_files(pdb_id_file, output_dir="pdb_files", fasta_output="pdb_sequences.fasta"):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read PDB IDs from file
    with open(pdb_id_file, "r") as pdb_input:
        pdb_ids = [line.strip() for line in pdb_input if line.strip()]

    pdbl = PDBList()
    parser = PDBParser(QUIET=True)
    fasta_records = []

    for pdb_id in pdb_ids:
        #pdb_id = pdb_id.lower()
        pdb_path = pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format="pdb")

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


# Example usage
#download_pdb_files("pdb_ids.txt")

"""
Method 04
input - single pdb file
output - text file with the extracted details

1. make the output file header 
2. read the pdb file
3. extract the number od chains, residues and number of atoms for each residue in protein 
4. save the results in one text file

"""
def details_from_pdb_file (input_file,output_file):
    parser = PDBParser(QUIET=True)

    with open (output_file, "w") as output:
    #make the output file header
        output.write("PDB_ID\t\t\tChain_ID\t\tResidue_ID\t\tResidue_Name\t\tNumber_of_atoms \n")

        pdb_id = os.path.splitext(os.path.basename(input_file))[0]

        #parse the PDB structure
        structure = parser.get_structure(pdb_id, input_file)

        num_chains = 0
        num_residues = 0
        num_atoms_total =0

        for model in structure:
            for chain in model:
                num_chains += 1
                for residue in chain :
                    if residue.id[0] == " ":  # Exclude heteroatoms (HETATM)
                        num_residues += 1
                        residue_id = residue.id[1]  # Residue sequence number
                        residue_name = residue.resname  # Residue three-letter code
                        num_atoms = len(residue)  # Number of atoms in this residue
                        output.write(f"{pdb_id}\t\t\t\t {chain.id}\t\t\t {residue_id}\t\t\t\t{residue_name}\t\t\t\t{num_atoms}\n")

                        num_atoms_total += num_atoms

        #write to the output file
        output.write(f"\nTotal Chains: {num_chains}\n")
        output.write(f"Total Residues: {num_residues}\n")
        output.write(f"Total atoms: {num_atoms_total}\n")

    print(f"Details saved in {output_file}")

#example usage
#details_from_pdb_file("pdb1gjh.ent", "single_file_detail.txt")


"""
Method 05
Input- The pdb files 
Output- number of chains , residues and atoms for each residue in separate text files

1. Read the pdb files one by one in the directory
2.initialize the variables to 0
2. extract the details and calculate from the files 
3.save each output in a separate text files

"""
def details_from_multiple_pdb (input_dir):
    parser = PDBParser(QUIET=True)

    for pdb_file in os.listdir(input_dir):
        if pdb_file.endswith(".ent"):
            pdb_id = pdb_file.split('.')[0]  # Extracting PDB ID
            file_path = os.path.join(input_dir, pdb_file)
            output_file = os.path.join(input_dir, f"{pdb_id}_details.txt")  # Output file for each PDB

            #parse the PDB structure
            structure = parser.get_structure(pdb_id, file_path)

            with open (output_file, "w") as output:
            #make the output file header
                output.write("PDB_ID\t\t\tChain_ID\t\tResidue_ID\t\tResidue_Name\t\tNumber_of_atoms \n")

                num_chains = 0
                num_residues = 0
                num_atoms_total =0

                for model in structure:
                    for chain in model:
                        num_chains += 1
                        for residue in chain :
                            if residue.id[0] == " ":  # Exclude heteroatoms (HETATM)
                                num_residues += 1
                                residue_id = residue.id[1]  # Residue sequence number
                                residue_name = residue.resname  # Residue three-letter code
                                num_atoms = len(residue)  # Number of atoms in this residue
                                output.write(f"{pdb_id}\t\t\t\t {chain.id}\t\t\t {residue_id}\t\t\t\t{residue_name}\t\t\t\t{num_atoms}\n")

                                num_atoms_total += num_atoms

                    #write to the output file
                output.write(f"\nTotal Chains: {num_chains}\n")
                output.write(f"Total Residues: {num_residues}\n")
                output.write(f"Total atoms: {num_atoms_total}\n")

            print(f"Details saved in {output_file}")

#example usage
#details_from_multiple_pdb("pdb_files")

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

def hetero_residues (pdb_file, output_file):
    parser = PDBParser(QUIET=True)

    with open (output_file, "w") as output:
        output.write(f"Hetero residues in {pdb_file}\n")

        pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
        #parse the structure
        structure = parser.get_structure(pdb_id, pdb_file)

        #loop through the models, chains, residues for find the hetero atoms
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != " ": #hetero atoms are the ones that residue.id is not equal to " "(blank)
                        residue_id = residue.id[1]
                        residue_name = residue.resname
                        chain_id = chain.id

                        output.write(f"{pdb_id}\t\t  Chain :{chain_id}\t\t  Residue: {residue_id} \t \t  Residue name {residue_name}\n")

    print(f"Hetero atoms are saved to {output_file}")

#example usage
#hetero_residues("pdb1jnx.ent", "hetero_residues.txt")

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
def get_atom_coordinates(pdb_file, atom_name, output_file):
    parser = PDBParser(QUIET=True)

    #open the output file
    with open(output_file, "w") as output:
        output.write (f"Coordinates for atoms\n")

        #get the pdb id
        pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]

        #parse the structure
        structure  = parser.get_structure(pdb_id, pdb_file)

        #looping through the models, chains, residues and atoms to get coordinates
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.get_name() == atom_name: #match atom by name
                            x,y,z = atom.get_coord() #get the coordinates
                            residue_id = residue.id[1] #get residue sequence number
                            residue_name = residue.resname #get residue name
                            chain_id = chain.id

                            output.write(f"{pdb_id} \t\t Chain : {chain_id}\t\t Residue_id: {residue_id}\t\t Residue name: {residue_name}\t\t Coordinates{atom_name}\t:{x:.3f}\t\t ,{y:.3f}\t\t ,{z:.3f}\t\n")

#example usage
#get_atom_coordinates("pdb1jnx.ent", "C", "Coordinates.txt")











