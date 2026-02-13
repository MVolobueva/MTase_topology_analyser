import os
import tempfile
import urllib.request
import subprocess

def download_structure(identifier, source='pdb'):
    """Download structure and generate DSSP"""
    identifier = identifier.strip()
    temp_dir = tempfile.mkdtemp()
    
    if source == 'pdb':
        url = f"https://files.rcsb.org/download/{identifier}.pdb"
        pdb_file = os.path.join(temp_dir, f"{identifier}.pdb")
        urllib.request.urlretrieve(url, pdb_file)
        
    elif source == 'alphafold':
        url = f"https://alphafold.ebi.ac.uk/files/AF-{identifier}-F1-model_v4.pdb"
        pdb_file = os.path.join(temp_dir, f"AF-{identifier}.pdb")
        urllib.request.urlretrieve(url, pdb_file)
    
    dssp_file = os.path.join(temp_dir, f"{identifier}.dssp")
    subprocess.run(['dssp', pdb_file, dssp_file], capture_output=True)
    
    return dssp_file

def parse_uploaded_file(uploaded_file):
    """Parse uploaded PDB file and generate DSSP"""
    temp_dir = tempfile.mkdtemp()
    
    pdb_file = os.path.join(temp_dir, uploaded_file.name)
    with open(pdb_file, 'wb') as f:
        f.write(uploaded_file.getvalue())
    
    dssp_file = os.path.join(temp_dir, 'uploaded.dssp')
    subprocess.run(['dssp', pdb_file, dssp_file], capture_output=True)
    
    return dssp_file