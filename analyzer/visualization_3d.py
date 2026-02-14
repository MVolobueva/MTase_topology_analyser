import numpy as np
import os
from scipy.spatial import distance_matrix
from .core import MTaseAnalyzer
import py3Dmol

def visualize_3d_structure(self, result, pdb_id=None, chain=None, pdb_file=None):
    """3D визуализация структуры с учетом цепи"""
    print(f"🔍 3D DEBUG - pdb_id: {pdb_id}, chain: {chain}, pdb_file: {pdb_file}")
    print(f"🔍 3D DEBUG - result keys: {result.keys() if result else 'None'}")
    
    if not result:
        print("Ошибка: нет результата анализа")
        return None

    if not self.motif_info:
        print("Ошибка: каталитический мотив не найден")
        return None

    motif_text = result.get('motif_text', 'Motif')
    motif_res_num = result.get('motif_res', 0)
    motif_chain = result.get('chain', 'A')


    if chain is None:
        chain = motif_chain
    elif chain != motif_chain:
        print(f"⚠️ Запрошена цепь {chain}, но мотив в цепи {motif_chain}")
        chain = motif_chain

    print(f"\n{'='*60}")
    print(f"3D ВИЗУАЛИЗАЦИЯ ТОПОЛОГИИ")
    print(f"Цепь анализа: '{chain}'")
    print(f"{'='*60}")

    # -----------------------------------------------------------------
    # 1. СБОР ЭЛЕМЕНТОВ ДЛЯ РАСКРАСКИ
    # -----------------------------------------------------------------
    elements = []
    full_path = result['full_path']
    s4_idx = result['s4_idx']
    s4_pos = full_path.index(s4_idx)
    path_map = result['path_map']
    strand_names = result.get('strand_names', {})

    # Добавляем тяжи из result['strands']
    for i, idx in enumerate(full_path):
        if idx in strand_names:
            s_name = strand_names[idx]
        else:
            s_name = f"S{4 - (i - s4_pos)}"
        
        s_range = result['strands'][idx]
        first_res_key = s_range[0]
        
        if first_res_key in self.res_data:
            strand_chain = self.res_data[first_res_key]['chain']
            if strand_chain == chain:
                s_start = self._get_res_num(s_range[0])
                s_end = self._get_res_num(s_range[-1])
                elements.append({
                    'name': s_name,
                    'start': s_start,
                    'end': s_end,
                    'type': 'strand',
                    'chain': chain
                })
                print(f"  Тяж {s_name}: {s_start}-{s_end} (цепь {chain})")

    # Добавляем спирали из result['helices']
    for h_keys in result['helices']:
        if len(h_keys) < self.MIN_HELIX_LENGTH:
            continue
            
        h_start = self._get_res_num(h_keys[0])
        h_end = self._get_res_num(h_keys[-1])

        # Проверяем цепь
        if h_keys[0] in self.res_data:
            helix_chain = self.res_data[h_keys[0]]['chain']
            if helix_chain != chain:
                continue

        if h_start in result['helix_sides']:
            side = result['helix_sides'][h_start]
            num_part = self._get_helix_number(h_start, h_end, path_map)
            name = self._get_helix_name(side, h_start, num_part)

            elements.append({
                'name': name,
                'start': h_start,
                'end': h_end,
                'type': 'helix',
                'side': side,
                'chain': chain
            })
            print(f"  Спираль {name}: {h_start}-{h_end} (цепь {chain}, {side})")

    elements.sort(key=lambda x: x['start'])

    if not elements:
        print("Ошибка: не найдены элементы для раскраски")
        return None

    # -----------------------------------------------------------------
    # 2. ЦВЕТА
    # -----------------------------------------------------------------
    strand_color = self.COLORS['strand']
    hu_color = self.COLORS['Hu']
    hd_color = self.COLORS['Hd']

    # -----------------------------------------------------------------
    # 3. СОЗДАНИЕ VIEWER (С ПОДДЕРЖКОЙ PDB ФАЙЛОВ)
    # -----------------------------------------------------------------
    try:
        if pdb_file and os.path.exists(pdb_file):
            # Загружаем локальный PDB файл
            with open(pdb_file, 'r') as f:
                pdb_data = f.read()
            view = py3Dmol.view(data=pdb_data, format='pdb')
            print(f"  Загружен локальный PDB файл: {pdb_file}")
        elif pdb_id:
            # Загружаем по PDB ID
            view = py3Dmol.view(query=f'pdb:{pdb_id}')
            print(f"  Загружена структура PDB: {pdb_id}")
        else:
            # Пустой viewer
            view = py3Dmol.view()
            print("  Используется пустой viewer")
    except Exception as e:
        print(f"  Ошибка загрузки: {e}")
        return None

    # -----------------------------------------------------------------
    # 4. БАЗОВЫЙ СТИЛЬ - ВЕСЬ БЕЛОК СВЕТЛО-СЕРЫЙ
    # -----------------------------------------------------------------
    # Сначала весь белок - светло-серый полупрозрачный
    view.setStyle({'cartoon': {'color': '#cccccc', 'opacity': 0.2}})
    
    # Потом текущая цепь - чуть ярче (НО НЕ ПЕРЕЗАПИСЫВАЕМ!)
    if chain:
        view.addStyle({'chain': chain}, {'cartoon': {'color': '#e0e0e0', 'opacity': 0.6}})
    
    print(f"  Базовый стиль: весь белок светло-серый, цепь {chain} выделена")

    # -----------------------------------------------------------------
    # 5. РАСКРАСКА ЭЛЕМЕНТОВ ПО ТОПОЛОГИИ
    # -----------------------------------------------------------------
    for elem in elements:
        selector = {'chain': chain, 'resi': f"{elem['start']}-{elem['end']}"}

        if elem['type'] == 'strand':
            color = strand_color
            view.addStyle(selector, {
                'cartoon': {
                    'color': color,
                    'arrows': True,
                    'opacity': 1.0,
                    'thickness': 0.8
                }
            })
        else:
            color = hu_color if elem['side'] == 'Hu' else hd_color
            view.addStyle(selector, {
                'cartoon': {
                    'color': color,
                    'opacity': 1.0,
                    'thickness': 0.8
                }
            })

        # Добавляем метку с именем
        center_res = (elem['start'] + elem['end']) // 2
        view.addLabel(
            elem['name'],
            {
                'fontSize': 12,
                'fontColor': 'black',
                'backgroundColor': 'white',
                'backgroundOpacity': 0.8,
                'borderColor': color,
                'borderWidth': 1
            },
            {'chain': chain, 'resi': center_res}
        )

    # -----------------------------------------------------------------
    # 6. ВЫДЕЛЕНИЕ КАТАЛИТИЧЕСКОГО МОТИВА
    # -----------------------------------------------------------------
    view.addStyle(
        {'chain': chain, 'resi': motif_res_num},
        {
            'stick': {
                'colorscheme': 'yellowCarbon',
                'radius': 0.3,
                'singleBonds': True
            }
        }
    )
    
    view.addLabel(
        f"{motif_text} ({motif_res_num})",
        {
            'fontSize': 11,
            'fontColor': 'black',
            'backgroundColor': 'yellow',
            'backgroundOpacity': 0.8,
            'borderColor': 'orange',
            'borderWidth': 1
        },
        {'chain': chain, 'resi': motif_res_num}
    )

    # -----------------------------------------------------------------
    # 7. МЕТКИ N И C-КОНЦОВ
    # -----------------------------------------------------------------
    if elements:
        n_elem = elements[0]
        view.addLabel(
            'N',
            {
                'fontSize': 14,
                'fontColor': '#2c3e50',
                'backgroundColor': 'yellow',
                'backgroundOpacity': 0.8,
                'borderColor': '#2c3e50',
                'borderWidth': 1
            },
            {'chain': chain, 'resi': n_elem['start']}
        )

        c_elem = elements[-1]
        view.addLabel(
            'C',
            {
                'fontSize': 14,
                'fontColor': '#2c3e50',
                'backgroundColor': 'yellow',
                'backgroundOpacity': 0.8,
                'borderColor': '#2c3e50',
                'borderWidth': 1
            },
            {'chain': chain, 'resi': c_elem['end']}
        )

    # -----------------------------------------------------------------
    # 8. НАСТРОЙКА ПРОСМОТРА
    # -----------------------------------------------------------------
    if chain:
        view.zoomTo({'chain': chain})

    # Поворачиваем камеру согласно твоей системе координат
    # north, east, up - из твоей coord_system
    if result.get('coord_system'):
        cs = result['coord_system']
        # Устанавливаем ориентацию камеры
        view.rotate(180, 'z')
        view.rotate(0, 'x')
        view.rotate(45, 'y')  # пример - настрой под свои оси
        #view.rotate(45, 'x')

    view.setBackgroundColor('white')
    view.render()
    
    print(f"\n✅ 3D визуализация готова для цепи {chain}")
    print(f"{'='*60}")
    
    return view

# Прикрепляем метод к классу
MTaseAnalyzer.visualize_3d_structure = visualize_3d_structure