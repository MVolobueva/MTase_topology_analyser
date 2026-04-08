#!/usr/bin/env python3
"""
MTase Batch Analyzer
Supports:
- PDB ID (type=pdb)
- AlphaFold ID (type=alphafold)
- Local PDB file (type=file, ID=path to file)

Outputs topology and structural class (A-F) for each chain
"""

import pandas as pd
import sys
import io
import re
import os
import tempfile
import urllib.request
import subprocess
import stat
import shutil
from analyzer import MTaseAnalyzer
from classifier import classify_topology

# Path to DSSP executable
DSSP_BIN = os.path.join(os.getcwd(), "mkdssp")

# Minimum strand length for reliable direction determination
MIN_STRAND_LENGTH = 3


def clean_pdb_file(pdb_path):
    """Removes DBREF and REMARK lines from PDB file"""
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    cleaned = [line for line in lines if not line.startswith(('DBREF', 'REMARK'))]
    
    with open(pdb_path, 'w') as f:
        f.writelines(cleaned)
    
    return pdb_path


def get_structure(id_value, type_value):
    """Downloads or opens structure depending on type"""
    temp_dir = tempfile.mkdtemp()
    
    if type_value == 'pdb':
        # Download from PDB
        url = f"https://files.rcsb.org/download/{id_value}.pdb"
        pdb_file = os.path.join(temp_dir, f"{id_value}.pdb")
        urllib.request.urlretrieve(url, pdb_file)
        
    elif type_value == 'alphafold':
        # Download from AlphaFold
        url = f"https://alphafold.ebi.ac.uk/files/AF-{id_value}-F1-model_v6.pdb"
        pdb_file = os.path.join(temp_dir, f"{id_value}.pdb")
        urllib.request.urlretrieve(url, pdb_file)
        
    elif type_value == 'file':
        # Use local file
        pdb_file = id_value
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"File not found: {pdb_file}")
        # For local files, we don't create a temp_dir
        return pdb_file, None
        
    else:
        raise ValueError(f"Unknown type: {type_value}")
    
    # Clean downloaded file
    if type_value != 'file':
        clean_pdb_file(pdb_file)
    
    return pdb_file, temp_dir


def run_dssp(pdb_file):
    """Runs DSSP and returns path to dssp file"""
    dssp_file = pdb_file.replace('.pdb', '.dssp')
    
    # Make DSSP executable
    if os.path.exists(DSSP_BIN):
        os.chmod(DSSP_BIN, os.stat(DSSP_BIN).st_mode | stat.S_IEXEC)
    
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = os.getcwd() + ":" + env.get("LD_LIBRARY_PATH", "")
    
    result = subprocess.run(
        [DSSP_BIN, pdb_file, dssp_file],
        capture_output=True,
        text=True,
        env=env
    )
    
    if result.returncode != 0:
        raise RuntimeError(f"DSSP Error: {result.stderr}")
    
    return dssp_file


def get_topology_string(analyzer, result):
    """Returns topology string same as in web application"""
    if not result:
        return "", "", ""
    
    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout
    
    analyzer.print_linear_topology_from_result(result)
    
    sys.stdout = old_stdout
    linear_topology = new_stdout.getvalue()
    
    lines = linear_topology.strip().split('\n')
    topology_line = ""
    
    for line in lines:
        if ('↑' in line or '↓' in line) and '[' in line and ']' in line:
            topology_line = line.strip()
            break
    
    if not topology_line:
        return "", "", ""
    
    # Remove distance annotations
    topology_line = re.sub(r'\s*\([\d.]+ Å\)', '', topology_line)
    full_topology = topology_line
    
    # Strands only (without Hd/Hu helices)
    strands_only = topology_line
    strands_only = re.sub(r'H[ud][_a-zA-Z0-9]*\[[\d-]+\]\s*—\s*', '', strands_only)
    strands_only = re.sub(r'^[Hh][ud][_a-zA-Z0-9]*\[[\d-]+\]\s*—\s*', '', strands_only)
    strands_only = re.sub(r'\s*—\s*[Hh][ud][_a-zA-Z0-9]*\[[\d-]+\]$', '', strands_only)
    
    # Extract directions in descending strand number order (S7 → S6 → ... → S-1)
    s_numbers = re.findall(r'S(-?\d+)\([↑↓]\)', full_topology)
    s_numbers_sorted = sorted(set(int(n) for n in s_numbers), reverse=True)
    
    directions = []
    for num in s_numbers_sorted:
        s_name = f"S{num}"
        match = re.search(f'{s_name}\([↑↓]\)', full_topology)
        if match:
            if '(↑)' in match.group():
                directions.append('↑')
            elif '(↓)' in match.group():
                directions.append('↓')
    
    directions_str = ''.join(directions)
    
    return full_topology, strands_only, directions_str


def analyze_structure(pdb_file):
    """Analyzes a single PDB structure"""
    temp_dirs = []
    
    try:
        # Run DSSP
        dssp_file = run_dssp(pdb_file)
        
        # Create analyzer
        analyzer = MTaseAnalyzer()
        if not analyzer.load_dssp(dssp_file):
            return None
        
        analyzer.find_all_strands()
        analyzer.build_sheet_adjacency()
        
        motifs = analyzer.find_all_motifs()
        motifs = analyzer.filter_motifs_by_topology(motifs)
        
        if not motifs:
            return None
        
        results = []
        for motif_data in motifs:
            chain = motif_data['chain']
            motif_text = motif_data['text']
            motif_res = motif_data['res']
            motif_position = f"{motif_res}-{motif_res + len(motif_text) - 1}"
            
            result = analyzer.analyze_topology(motif_data=motif_data)
            if not result:
                continue
            
            full_topology, strands_only, directions = get_topology_string(analyzer, result)
            
            if not full_topology:
                continue
            
            # Classify topology
            classification = classify_topology(full_topology, motif_text)
            
            results.append({
                'chain': chain,
                'found_motif': motif_text,
                'found_motif_position': motif_position,
                'full_secondary_elements': full_topology,
                'strands_secondary_elements': strands_only,
                'strand_directions': directions,
                'class': classification['class'],
                'class_confidence': classification['confidence'],
                'class_reasons': classification['reasons'],
                'strand_order': classification['strand_order'],
                'has_s0': classification['has_s0'],
                'has_s-1': classification['has_s-1'],
                'gap_s6_s7': classification['gap_s6_s7'],
                'has_n_helix': classification['has_n_helix'],
                'has_c_helix': classification['has_c_helix']
            })
        
        return results
    
    except Exception as e:
        print(f"  ❌ Error: {str(e)}")
        return None


def main():
    print("=" * 70)
    print("MTase Batch Analyzer")
    print("=" * 70)
    print("Input CSV format: ID,Type")
    print("  Type: pdb, alphafold, file")
    print("=" * 70)
    
    # Parse command line arguments
    input_file = 'input.csv'
    output_file = 'output.csv'
    
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    # Load input CSV
    try:
        df = pd.read_csv(input_file)
        print(f"\n✅ Loaded {len(df)} entries from {input_file}")
    except FileNotFoundError:
        print(f"\n❌ File {input_file} not found!")
        print("   Please create input.csv with columns: ID,Type")
        sys.exit(1)
    
    # Validate columns
    if 'ID' not in df.columns or 'Type' not in df.columns:
        print("\n❌ CSV must have columns: ID, Type")
        print("   Example:")
        print("   ID,Type")
        print("   3S1S,pdb")
        print("   A0A7R8ZSU6,alphafold")
        print("   /path/to/file.pdb,file")
        sys.exit(1)
    
    # Results container
    all_results = []
    
    print("\n📊 Analyzing structures...")
    print("-" * 70)
    
    for idx, row in df.iterrows():
        id_value = str(row['ID']).strip()
        type_value = str(row['Type']).strip().lower()
        
        print(f"\n🔬 {idx+1}/{len(df)}: {id_value} ({type_value})")
        
        pdb_file = None
        temp_dir = None
        
        try:
            # Get structure
            pdb_file, temp_dir = get_structure(id_value, type_value)
            
            # Analyze
            results = analyze_structure(pdb_file)
            
            if results:
                for res in results:
                    all_results.append({
                        'source_id': id_value,
                        'source_type': type_value,
                        'chain': res['chain'],
                        'found_motif': res['found_motif'],
                        'found_motif_position': res['found_motif_position'],
                        'full_secondary_elements': res['full_secondary_elements'],
                        'strands_secondary_elements': res['strands_secondary_elements'],
                        'strand_directions': res['strand_directions'],
                        'class': res['class'],
                        'class_confidence': res['class_confidence'],
                        'class_reasons': res['class_reasons'],
                        'strand_order': res['strand_order'],
                        'has_s0': res['has_s0'],
                        'has_s-1': res['has_s-1'],
                        'gap_s6_s7': res['gap_s6_s7'],
                        'has_n_helix': res['has_n_helix'],
                        'has_c_helix': res['has_c_helix']
                    })
                print(f"  ✅ Found {len(results)} chain(s)")
            else:
                all_results.append({
                    'source_id': id_value,
                    'source_type': type_value,
                    'chain': None,
                    'found_motif': None,
                    'found_motif_position': None,
                    'full_secondary_elements': None,
                    'strands_secondary_elements': None,
                    'strand_directions': None,
                    'class': None,
                    'class_confidence': None,
                    'class_reasons': None,
                    'strand_order': None,
                    'has_s0': None,
                    'has_s-1': None,
                    'gap_s6_s7': None,
                    'has_n_helix': None,
                    'has_c_helix': None
                })
                print(f"  ⚠️ No motifs found")
            
            # Clean up temporary directory
            if temp_dir and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
                
        except Exception as e:
            print(f"  ❌ Error: {str(e)}")
            all_results.append({
                'source_id': id_value,
                'source_type': type_value,
                'chain': None,
                'found_motif': None,
                'found_motif_position': None,
                'full_secondary_elements': None,
                'strands_secondary_elements': None,
                'strand_directions': None,
                'class': None,
                'class_confidence': None,
                'class_reasons': None,
                'strand_order': None,
                'has_s0': None,
                'has_s-1': None,
                'gap_s6_s7': None,
                'has_n_helix': None,
                'has_c_helix': None,
                'error': str(e)
            })
            
            # Clean up on error
            if temp_dir and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
    
    # Create output DataFrame
    df_out = pd.DataFrame(all_results)
    
    # Save results
    df_out.to_csv(output_file, index=False, encoding='utf-8-sig')
    
    # Print summary
    print("\n" + "=" * 70)
    print(f"✅ Done! Saved to {output_file}")
    print(f"📊 Processed: {len(df)} entries")
    
    chains_with_motifs = len([r for r in all_results if r['found_motif']])
    print(f"📊 Chains with motifs: {chains_with_motifs}")
    
    # Class distribution
    class_counts = df_out['class'].value_counts()
    if len(class_counts) > 0:
        print(f"\n📊 Class distribution:")
        for cls, count in class_counts.items():
            if cls:
                print(f"   Class {cls}: {count}")
    
    print("=" * 70)


if __name__ == "__main__":
    main()