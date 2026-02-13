import numpy as np
from .core import MTaseAnalyzer

def visualize_3d_structure(self, result, pdb_id=None, chain=None):
    """3D визуализация с раскраской по топологии"""
    if not result:
        print("Ошибка: нет результата анализа")
        return None

    try:
        import py3Dmol
    except ImportError:
        print("Ошибка: требуется установка py3Dmol: pip install py3Dmol")
        return None

    if not self.motif_info:
        print("Ошибка: каталитический мотив не найден")
        return None

    motif_res = self.motif_info['res']
    motif_chain = self.motif_info.get('chain', 'A')
    motif_key = f"{motif_chain}:{motif_res}"
    motif_text = self.motif_info['text']

    if motif_key in self.res_data:
        motif_res_num = self.res_data[motif_key]['res_num']
        motif_chain = self.res_data[motif_key]['chain']
    else:
        motif_res_num = motif_res

    if chain is None:
        chain = motif_chain
    elif chain != motif_chain:
        print(f"⚠️ Запрошена цепь {chain}, но мотив в цепи {motif_chain}")
        chain = motif_chain

    print(f"\n{'='*60}")
    print(f"3D ВИЗУАЛИЗАЦИЯ ТОПОЛОГИИ")
    print(f"PDB: {pdb_id or 'локальный файл'}")
    print(f"Цепь анализа: '{chain}'")
    print(f"{'='*60}")

    # -----------------------------------------------------------------
    # 1. СБОР ЭЛЕМЕНТОВ ДЛЯ РАСКРАСКИ - ТОЛЬКО ИЗ result!
    # -----------------------------------------------------------------
    elements = []
    full_path = result['full_path']
    s4_idx = result['s4_idx']
    s4_pos = full_path.index(s4_idx)
    path_map = result['path_map']
    strand_names = result['strand_names']
    
    # ✅ Используем result['strands'] - это уже отфильтрованные тяжи!
    working_strands = result['strands']
    
    # ✅ Используем result['helices'] - это уже отфильтрованные спирали!
    working_helices = result['helices']

    # Добавляем тяжи
    for i, idx in enumerate(full_path):
        if idx in strand_names:
            s_name = strand_names[idx]
        else:
            s_name = f"S{4 - (i - s4_pos)}"
        
        s_range = working_strands[idx]
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

    # Добавляем спирали - БЕРЕМ ИЗ result['helices']!
    unique_helices_3d = {}
    for h_keys in working_helices:  # ✅ ИСПРАВЛЕНО!
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
            helix_key = f"{h_start}-{h_end}"
            if helix_key not in unique_helices_3d:
                side = result['helix_sides'][h_start]
                num_part = self._get_helix_number(h_start, h_end, path_map)
                name = self._get_helix_name(side, h_start, num_part)
                unique_helices_3d[helix_key] = (name, side, h_start, h_end)

    for helix_key, (name, side, h_start, h_end) in unique_helices_3d.items():
        elements.append({
            'name': name,
            'start': h_start,
            'end': h_end,
            'type': 'helix',
            'side': side,
            'chain': chain
        })
        print(f"  Спираль {name}: {h_start}-{h_end} (цепь {chain}, {side})")

    # Сортируем по номеру остатка
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
    # 3. СОЗДАНИЕ VIEWER
    # -----------------------------------------------------------------
    try:
        if pdb_id:
            view = py3Dmol.view(query=f'pdb:{pdb_id}')
            print(f"  Загружена структура PDB: {pdb_id}")
        else:
            view = py3Dmol.view()
            print("  Используется локальный PDB файл")
    except Exception as e:
        print(f"  Ошибка загрузки PDB: {e}")
        return None

    # -----------------------------------------------------------------
    # 4. БАЗОВЫЙ СТИЛЬ
    # -----------------------------------------------------------------
    view.setStyle({'cartoon': {'color': '#e0e0e0', 'opacity': 0.2}})
    view.setStyle({'chain': chain}, {'cartoon': {'color': '#e0e0e0', 'opacity': 0.6}})

    # -----------------------------------------------------------------
    # 5. РАСКРАСКА ЭЛЕМЕНТОВ
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
                'borderWidth': 2
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
            },
            'sphere': {
                'colorscheme': 'yellowCarbon',
                'radius': 0.5,
                'opacity': 0.6
            }
        }
    )
    
    view.addLabel(
        f"{motif_text} ({motif_res_num})",
        {
            'fontSize': 14,
            'fontColor': 'black',
            'backgroundColor': 'yellow',
            'backgroundOpacity': 0.9,
            'borderColor': 'orange',
            'borderWidth': 2
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
                'fontSize': 16,
                'fontColor': '#2c3e50',
                'backgroundColor': 'yellow',
                'backgroundOpacity': 0.9,
                'borderColor': '#2c3e50',
                'borderWidth': 2
            },
            {'chain': chain, 'resi': n_elem['start']}
        )

        c_elem = elements[-1]
        view.addLabel(
            'C',
            {
                'fontSize': 16,
                'fontColor': '#2c3e50',
                'backgroundColor': 'yellow',
                'backgroundOpacity': 0.9,
                'borderColor': '#2c3e50',
                'borderWidth': 2
            },
            {'chain': chain, 'resi': c_elem['end']}
        )

    # -----------------------------------------------------------------
    # 8. НАСТРОЙКА ПРОСМОТРА
    # -----------------------------------------------------------------
    view.zoomTo({'chain': chain})
    
    # Добавляем поверхность для контекста
    view.addSurface(py3Dmol.VDW, {'opacity': 0.1, 'color': 'gray'})
    
    # Настройки отображения
    view.setBackgroundColor('white')
    view.render()
    
    print(f"\n✅ 3D визуализация готова для цепи {chain}")
    print(f"{'='*60}")
    view.show()
    
    return view

MTaseAnalyzer.visualize_3d_structure = visualize_3d_structure