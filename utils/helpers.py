import os
import tempfile
import urllib.request
import subprocess
import shutil
import streamlit as st  # Добавили для вывода ошибок на экран

def get_dssp_command():
    """Определяет, какая команда доступна в системе: mkdssp (новое) или dssp (старое)"""
    if shutil.which("mkdssp"):
        return "mkdssp"
    return "dssp"

def download_structure(identifier, source='pdb'):
    """Download structure and generate DSSP"""
    identifier = identifier.strip().upper() # Переводим в верхний регистр
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
        st.error(f"Ошибка загрузки файла: {e}")
        return None
    
    dssp_file = os.path.join(temp_dir, f"{identifier}.dssp")
    cmd = get_dssp_command()
    
    # Запуск с флагами для версии 4.4.10 (Streamlit Cloud)
    result = subprocess.run([cmd, '--classic', '--not-use-dictionary', pdb_file, dssp_file], 
                           capture_output=True, text=True)
    
    # Если команда не сработала, выводим ошибку прямо в Streamlit
    if result.returncode != 0:
        st.warning(f"Попытка запуска {cmd} без флагов...")
        result = subprocess.run([cmd, pdb_file, dssp_file], 
                               capture_output=True, text=True)
    
    if result.returncode != 0:
        # ВЫВОД ОШИБКИ НА ЭКРАН
        st.error(f"DSSP Error (Code {result.returncode}):")
        st.code(result.stderr) 
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
    
    pdb_file = os.path.join(temp_dir, uploaded_file.name)
    with open(pdb_file, 'wb') as f:
        f.write(uploaded_file.getvalue())
    
    dssp_file = os.path.join(temp_dir, 'uploaded.dssp')
    cmd = get_dssp_command()
    
    # Запуск с флагом --not-use-dictionary
    result = subprocess.run([cmd, '--classic', '--not-use-dictionary', pdb_file, dssp_file], 
                           capture_output=True, text=True)
    
    if result.returncode != 0:
        st.error(f"DSSP Error on Upload:")
        st.code(result.stderr)
        return None
    
    return {
        'dssp': dssp_file,
        'pdb': pdb_file
    }
