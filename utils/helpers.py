def download_structure(identifier, source='pdb'):
    """Download structure and generate DSSP"""
    import os
    import tempfile
    import urllib.request
    import subprocess
    
    identifier = identifier.strip()
    temp_dir = tempfile.mkdtemp()
    
    if source == 'pdb':
        url = f"https://files.rcsb.org/download/{identifier}.pdb"
        pdb_file = os.path.join(temp_dir, f"{identifier}.pdb")
        print(f"Downloading PDB: {url}")
        urllib.request.urlretrieve(url, pdb_file)
        
    elif source == 'alphafold':
        # AlphaFold URL format
        url = f"https://alphafold.ebi.ac.uk/files/AF-{identifier}-F1-model_v4.pdb"
        pdb_file = os.path.join(temp_dir, f"AF-{identifier}.pdb")
        print(f"Downloading AlphaFold: {url}")
        urllib.request.urlretrieve(url, pdb_file)
    
    # Run DSSP
    dssp_file = os.path.join(temp_dir, f"{identifier}.dssp")
    result = subprocess.run(['dssp', pdb_file, dssp_file], 
                           capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"DSSP error: {result.stderr}")
        return None
        
    return dssp_file