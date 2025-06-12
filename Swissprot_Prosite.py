"""
Method 01
Fetches SwissProt sequences for given UniProt IDs and saves each in a separate FASTA file.
input_file : Path to the text file containing SwissProt IDs (one per line).
output_dir : Directory where output FASTA files will be saved.

1.

"""

from Bio import ExPASy, SwissProt
import os
import time

def fetch_swissprot_records(input_file, output_dir):

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

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
        fasta_filename = os.path.join(output_dir, f"{uniprot_id}.fasta")
        with open(fasta_filename, "w") as fasta_file:
            fasta_file.write(fasta_content)

        print(f"Saved: {fasta_filename}")
        time.sleep(1)  # To avoid overwhelming UniProt servers

    print(f"All sequences saved in {output_dir}")

# Example usage
#fetch_swissprot_records("swissprotIDs.txt", "fasta_sequences")

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
from Bio import SeqIO
from Bio.ExPASy import ScanProsite

def scan_prosite_from_fasta(input_dir, output_file):

    with open(output_file, "w") as output:
        # Write header
        output.write("Sequence_ID\t\t\tDomain_ID\t\tDomain_Name\t\tStart\t\tEnd\n")

        for fasta_file in os.listdir(input_dir):
            if fasta_file.endswith (".fasta"):
                file_path =os.path.join (input_dir,fasta_file)

                for record in SeqIO.parse(file_path, "fasta"):
                    sequence = record.seq

                # Scan the sequence using ScanProsite
                    handle = ScanProsite.scan(sequence)
                    results = ScanProsite.read(handle)

                    #print(f"ScanProsite results:\n{results[:1000]}")
                if results:
                    for result in results:
                        # Extracting domain information
                        domain_id = result.get('signature_ac', '')  # Domain ID (signature accession)
                        domain_name = result.get('signature_ac', '')  # Use an empty string if no name is available
                        start = result.get('start', '') #coordinates start
                        end = result.get('stop', '') #coordinates end

                        # Write results to the output file
                        output.write(f"{record.id}\t\t{domain_id}\t\t{domain_name}\t\t{start}\t\t\t{end}\n")


# Example usage
#scan_prosite_from_fasta("fasta_sequences", "domain_scan_results.txt")



