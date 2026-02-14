import os
import tempfile
import urllib.request
import subprocess

def download_structure(identifier, source='pdb'):
    """Download structure and generate DSSP"""
    identifier = identifier.strip()
    temp_dir = tempfile.mkdtemp()
    
    try:
        if source == 'pdb':
            url = f"https://files.rcsb.org/download/{identifier}.pdb"
            print(f"Downloading PDB: {url}")
        elif source == 'alphafold':
            url = f"https://alphafold.ebi.ac.uk//files/AF-{identifier}-F1-model_v6.pdb"
            print(f"Downloading AlphaFold: {url}")
        
        pdb_file = os.path.join(temp_dir, f"{identifier}.pdb")
        urllib.request.urlretrieve(url, pdb_file)
        
    except Exception as e:
        print(f"Error downloading: {e}")
        return None
    
    dssp_file = os.path.join(temp_dir, f"{identifier}.dssp")
    result = subprocess.run(['dssp', pdb_file, dssp_file], 
                           capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"DSSP error: {result.stderr}")
        return None
        
    return {
        'dssp': dssp_file,
        'pdb': pdb_file,
        'id': identifier,
        'source': source
    }

# 👇 ДОБАВЬ ЭТУ ФУНКЦИЮ
def parse_uploaded_file(uploaded_file):
    """Parse uploaded PDB file and generate DSSP"""
    import tempfile
    import os
    import subprocess
    
    temp_dir = tempfile.mkdtemp()
    
    # Save uploaded file
    pdb_file = os.path.join(temp_dir, uploaded_file.name)
    with open(pdb_file, 'wb') as f:
        f.write(uploaded_file.getvalue())
    
    # Generate DSSP
    dssp_file = os.path.join(temp_dir, 'uploaded.dssp')
    result = subprocess.run(['dssp', pdb_file, dssp_file], 
                           capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"DSSP error: {result.stderr}")
        return None
    
    # 👇 ВОЗВРАЩАЕМ СЛОВАРЬ!
    return {
        'dssp': dssp_file,
        'pdb': pdb_file
    }