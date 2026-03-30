import os
import tempfile
import urllib.request
import subprocess
import shutil

def get_dssp_command():
    """Определяет, какая команда доступна в системе: mkdssp (новое) или dssp (старое)"""
    if shutil.which("mkdssp"):
        return "mkdssp"
    return "dssp"

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

            print(f"Downloading AlphaFold: {url}")
        
        pdb_file = os.path.join(temp_dir, f"{identifier}.pdb")
        urllib.request.urlretrieve(url, pdb_file)
        
    except Exception as e:
        print(f"Error downloading: {e}")
        return None
    
    dssp_file = os.path.join(temp_dir, f"{identifier}.dssp")
    
    # Используем автоматическое определение команды
    cmd = get_dssp_command()
    
    # Флаг --classic нужен для совместимости с парсерами, если версия DSSP >= 4.0
    # Добавьте флаг '--not-use-dictionary' внутрь списка аргументов:
    result = subprocess.run([cmd, '--classic', '--not-use-dictionary', pdb_file, dssp_file], capture_output=True, text=True)

    
    # Если --classic не поддерживается (совсем старая версия), пробуем без него
    if result.returncode != 0:
        result = subprocess.run([cmd, pdb_file, dssp_file], 
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

def parse_uploaded_file(uploaded_file):
    """Parse uploaded PDB file and generate DSSP"""
    temp_dir = tempfile.mkdtemp()
    
    # Save uploaded file
    pdb_file = os.path.join(temp_dir, uploaded_file.name)
    with open(pdb_file, 'wb') as f:
        f.write(uploaded_file.getvalue())
    
    # Generate DSSP
    dssp_file = os.path.join(temp_dir, 'uploaded.dssp')
    
    cmd = get_dssp_command()
    
    # Пробуем запустить с флагом --classic
    result = subprocess.run([cmd, '--classic', '--not-use-dictionary', pdb_file, dssp_file], 
                           capture_output=True, text=True)
    
    if result.returncode != 0:
        result = subprocess.run([cmd, pdb_file, dssp_file], 
                               capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"DSSP error: {result.stderr}")
        return None
    
    return {
        'dssp': dssp_file,
        'pdb': pdb_file
    }
