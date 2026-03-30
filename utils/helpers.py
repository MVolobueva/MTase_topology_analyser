import os
import tempfile
import urllib.request
import subprocess
import stat
import streamlit as st

# 1. Указываем путь к твоему файлу mkdssp, который лежит в корне проекта
DSSP_BIN = os.path.join(os.getcwd(), "mkdssp")

def download_structure(identifier, source='pdb'):
    """Загрузка структуры и запуск локального DSSP"""
    identifier = identifier.strip().upper()
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
        st.error(f"Ошибка загрузки: {e}")
        return None
    
    dssp_file = os.path.join(temp_dir, f"{identifier}.dssp")

    # 2. Даем права на запуск файла mkdssp (обязательно для Linux/Streamlit)
    if os.path.exists(DSSP_BIN):
        try:
            current_mode = os.stat(DSSP_BIN).st_mode
            os.chmod(DSSP_BIN, current_mode | stat.S_IEXEC)
        except Exception as e:
            print(f"Предупреждение chmod: {e}")
    else:
        st.error("Файл mkdssp не найден в корне проекта!")
        return None

    # Создаем копию окружения и добавляем текущую папку в путь поиска библиотек
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = os.getcwd() + ":" + env.get("LD_LIBRARY_PATH", "")



    result = subprocess.run(
        [DSSP_BIN, pdb_file, dssp_file], 
        capture_output=True, 
        text=True, env=env
    )
    
    if result.returncode != 0:
        st.error(f"DSSP Error: {result.stderr}")
        return None
        
    return {
        'dssp': dssp_file,
        'pdb': pdb_file,
        'id': identifier,
        'source': source
    }

def parse_uploaded_file(uploaded_file):
    """Обработка загруженного PDB файла"""
    temp_dir = tempfile.mkdtemp()
    
    pdb_file = os.path.join(temp_dir, uploaded_file.name)
    with open(pdb_file, 'wb') as f:
        f.write(uploaded_file.getvalue())
    
    dssp_file = os.path.join(temp_dir, 'uploaded.dssp')

    # Даем права на запуск
    if os.path.exists(DSSP_BIN):
        os.chmod(DSSP_BIN, os.stat(DSSP_BIN).st_mode | stat.S_IEXEC)
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = os.getcwd() + ":" + env.get("LD_LIBRARY_PATH", "")
    # Запуск
    result = subprocess.run(
        [DSSP_BIN, pdb_file, dssp_file], 
        capture_output=True, 
        text=True, env=env
    )
    
    if result.returncode != 0:
        st.error(f"DSSP Error: {result.stderr}")
        return None
    
    return {
        'dssp': dssp_file,
        'pdb': pdb_file
    }
