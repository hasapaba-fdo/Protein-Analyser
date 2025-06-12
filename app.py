from flask import Flask, render_template, request, redirect, send_file
import os
from Protein_Analyzer import ProteinAnalyzer  # your existing class

app = Flask(__name__)
analyzer = ProteinAnalyzer(output_dir="web_output")

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/fetch_swissprot', methods=['POST'])
def fetch_swissprot():
    file = request.files['swissprot_file']
    path = os.path.join("uploads", file.filename)
    file.save(path)
    analyzer.fetch_swissprot_records(path)
    return "SwissProt sequences fetched and saved!"

@app.route('/scan_prosite', methods=['GET'])
def scan_prosite():
    analyzer.scan_prosite_from_fasta()
    return send_file("web_output/domain_scan_results.txt")

@app.route('/download_pdb_files', method = ['POST'])
def download_pdb_files():
    file = request.files['pdb_id_file']
    path = os.path.join("uploads", file.filename)
    file.save(path)
    analyzer.download_pdb_files(path)
    return "PDB files structures and sequences download"

@app.route('/')

if __name__ == '__main__':
    os.makedirs("uploads", exist_ok=True)
    app.run(debug=True)
