#!/usr/bin/env python3
"""
Скрипт для автоматического обновления таблицы MTases
Анализирует структуры и добавляет колонки:
- chain: цепь анализа
- found_motif: найденный каталитический мотив
- found_motif_position: позиция мотива
- full_secondary_elements: полная топология (Hd/Hu + S со стрелками)
- strands_secondary_elements: только тяжи (S без стрелок)
- strand_directions: направления тяжей в порядке убывания номеров (S7→S6→...→S-1)
"""

import pandas as pd
import sys
import io
import re
from analyzer import MTaseAnalyzer
from utils.helpers import download_structure

# Минимальная длина тяжа для надежного определения направления
MIN_STRAND_LENGTH = 3

def get_topology_string(analyzer, result):
    """
    Возвращает топологию так же, как в веб-приложении
    Использует print_linear_topology_from_result для получения правильного порядка N→C
    """
    if not result:
        return "", "", ""
    
    # Захватываем вывод print_linear_topology_from_result
    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout
    
    analyzer.print_linear_topology_from_result(result)
    
    sys.stdout = old_stdout
    linear_topology = new_stdout.getvalue()
    
    # Извлекаем строку с элементами (содержит ↑ ↓ и [])
    lines = linear_topology.strip().split('\n')
    topology_line = ""
    
    for line in lines:
        if ('↑' in line or '↓' in line) and '[' in line and ']' in line:
            topology_line = line.strip()
            break
    
    if not topology_line:
        return "", "", ""
    
    # Убираем расстояния вида (15.8 Å)
    topology_line = re.sub(r'\s*\([\d.]+ Å\)', '', topology_line)
    
    # Полная топология (со спиралями и стрелками)
    full_topology = topology_line
    
    # Создаем версию только с тяжами (без Hd/Hu)
    strands_only = topology_line
    # Убираем Hd/Hu спирали
    strands_only = re.sub(r'H[ud][_a-zA-Z0-9]*\[[\d-]+\]\s*—\s*', '', strands_only)
    strands_only = re.sub(r'^[Hh][ud][_a-zA-Z0-9]*\[[\d-]+\]\s*—\s*', '', strands_only)
    strands_only = re.sub(r'\s*—\s*[Hh][ud][_a-zA-Z0-9]*\[[\d-]+\]$', '', strands_only)
    # Убираем BIG LOOP из версии с тяжами (опционально)
    strands_only = re.sub(r'—\s*🔴 BIG LOOP 🔴\s*—', ' — ', strands_only)
    
    # Направления в порядке убывания номеров тяжей (S7 → S6 → S5 → ... → S-1)
    # Получаем все номера тяжей из full_topology
    s_numbers = re.findall(r'S(-?\d+)\([↑↓]\)', full_topology)
    # Сортируем по убыванию (от большего к меньшему)
    s_numbers_sorted = sorted(set(int(n) for n in s_numbers), reverse=True)
    
    directions = []
    for num in s_numbers_sorted:
        s_name = f"S{num}"
        # Ищем этот тяж в full_topology
        match = re.search(f'{s_name}\([↑↓]\)', full_topology)
        if match:
            if '(↑)' in match.group():
                directions.append('↑')
            elif '(↓)' in match.group():
                directions.append('↓')
    
    directions_str = ''.join(directions)
    
    return full_topology, strands_only, directions_str

def analyze_all_chains(pdb_id):
    """
    Анализирует ВСЕ цепи в структуре и возвращает список результатов
    """
    results = []
    
    try:
        analyzer = MTaseAnalyzer()
        
        # Загружаем структуру
        result_files = download_structure(pdb_id, source='pdb')
        
        if not result_files or not analyzer.load_dssp(result_files['dssp']):
            print(f"  ❌ Failed to load {pdb_id}")
            return results
        
        # Находим вторичные структуры
        analyzer.find_all_strands()
        analyzer.build_sheet_adjacency()
        
        # Находим и фильтруем мотивы
        motifs = analyzer.find_all_motifs()
        motifs = analyzer.filter_motifs_by_topology(motifs)
        
        if not motifs:
            print(f"  ⚠️ No motifs found in {pdb_id}")
            return results
        
        print(f"  🔍 Found {len(motifs)} motifs in {pdb_id}")
        
        # Анализируем каждый мотив (каждую цепь)
        for motif_data in motifs:
            chain = motif_data['chain']
            motif_text = motif_data['text']
            motif_res = motif_data['res']
            motif_position = f"{motif_res}-{motif_res + len(motif_text) - 1}"
            
            # Анализируем топологию
            result = analyzer.analyze_topology(motif_data=motif_data)
            if not result:
                print(f"    ⚠️ No topology for chain {chain}")
                continue
            
            # Проверяем длину тяжей и предупреждаем о коротких
            short_strands = []
            for idx in result['full_path']:
                strand = result['strands'][idx]
                if len(strand) < MIN_STRAND_LENGTH:
                    name = result['strand_names'].get(idx, f"S{idx}")
                    short_strands.append(f"{name}({len(strand)} residues)")
            
            if short_strands:
                print(f"    ⚠️ Chain {chain}: short strands - {', '.join(short_strands)}")
            
            # Получаем топологию через функцию, которая использует print_linear_topology_from_result
            full_topology, strands_only, directions = get_topology_string(analyzer, result)
            
            if not full_topology:
                print(f"    ⚠️ Could not parse topology for chain {chain}")
                continue
            
            results.append({
                'chain': chain,
                'found_motif': motif_text,
                'found_motif_position': motif_position,
                'full_secondary_elements': full_topology,
                'strands_secondary_elements': strands_only,
                'strand_directions': directions
            })
            
            print(f"    ✅ Chain {chain}: motif {motif_text} at {motif_position}, directions: {directions}")
        
        return results
    
    except Exception as e:
        print(f"  ❌ Error analyzing {pdb_id}: {str(e)}")
        return results

def main():
    print("=" * 70)
    print("MTase Table Topology Updater")
    print("=" * 70)
    print("Добавляет колонки:")
    print("  - chain")
    print("  - found_motif")
    print("  - found_motif_position")
    print("  - full_secondary_elements")
    print("  - strands_secondary_elements")
    print("  - strand_directions")
    print("=" * 70)
    
    input_file = 'mtases_table.csv'
    output_file = 'mtases_table_updated.csv'
    
    try:
        df = pd.read_csv(input_file)
        print(f"\n✅ Loaded {len(df)} entries from {input_file}")
    except FileNotFoundError:
        print(f"\n❌ File {input_file} not found!")
        print("Please make sure mtases_table.csv is in the current directory")
        sys.exit(1)
    
    # Создаем новый список для расширенной таблицы
    new_rows = []
    
    print("\n📊 Analyzing structures...")
    print("-" * 70)
    
    for idx, row in df.iterrows():
        pdb_id = row['Repr. PDB code']
        if pd.notna(pdb_id) and str(pdb_id).strip():
            print(f"\n🔬 {idx+1}/{len(df)}: {pdb_id}")
            
            # Анализируем все цепи
            chain_results = analyze_all_chains(str(pdb_id).strip())
            
            if chain_results:
                # Для каждой цепи создаем отдельную строку
                for chain_result in chain_results:
                    new_row = row.to_dict()
                    new_row['chain'] = chain_result['chain']
                    new_row['found_motif'] = chain_result['found_motif']
                    new_row['found_motif_position'] = chain_result['found_motif_position']
                    new_row['full_secondary_elements'] = chain_result['full_secondary_elements']
                    new_row['strands_secondary_elements'] = chain_result['strands_secondary_elements']
                    new_row['strand_directions'] = chain_result['strand_directions']
                    new_rows.append(new_row)
            else:
                # Если нет мотивов, добавляем строку с пустыми значениями
                new_row = row.to_dict()
                new_row['chain'] = None
                new_row['found_motif'] = None
                new_row['found_motif_position'] = None
                new_row['full_secondary_elements'] = None
                new_row['strands_secondary_elements'] = None
                new_row['strand_directions'] = None
                new_rows.append(new_row)
                print(f"  ⚠️ No motifs found in any chain")
        else:
            new_row = row.to_dict()
            new_row['chain'] = None
            new_row['found_motif'] = None
            new_row['found_motif_position'] = None
            new_row['full_secondary_elements'] = None
            new_row['strands_secondary_elements'] = None
            new_row['strand_directions'] = None
            new_rows.append(new_row)
    
    # Создаем новый DataFrame
    df_new = pd.DataFrame(new_rows)
    
    # Сохраняем с UTF-8 BOM для совместимости с Excel
    df_new.to_csv(output_file, index=False, encoding='utf-8-sig')
    
    print("\n" + "=" * 70)
    print(f"✅ Done! Saved to {output_file}")
    print(f"📊 Original rows: {len(df)}")
    print(f"📊 New rows (after expanding chains): {len(df_new)}")
    print(f"\n📋 Added columns:")
    print("   - chain")
    print("   - found_motif")
    print("   - found_motif_position")
    print("   - full_secondary_elements")
    print("   - strands_secondary_elements")
    print("   - strand_directions")
    print("=" * 70)

if __name__ == "__main__":
    main()
