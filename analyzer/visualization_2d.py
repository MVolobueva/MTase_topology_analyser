import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from scipy.spatial import distance_matrix
from .core import MTaseAnalyzer
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

def visualize_topology_interactive(self, result):
    """ИНТЕРАКТИВНАЯ 2D визуализация с Plotly"""
    if not result:
        return None
    
    full_path = result['full_path']
    path_map = result['path_map']
    strand_names = result['strand_names']
    s4_idx = result['s4_idx']
    s4_pos = full_path.index(s4_idx)
    v4 = result['v4']
    
    elements = []
    
    # Сбор элементов - ТЯЖИ
    for i, idx in enumerate(full_path):
        if idx in strand_names:
            s_name = strand_names[idx]
        else:
            s_name = f"S{4 - (i - s4_pos)}"
        s_range = result['strands'][idx]
        s_start = self._get_res_num(s_range[0])
        s_end = self._get_res_num(s_range[-1])
        vi = self._get_strand_vector(s_range)
        direction = "UP" if np.dot(vi, v4) > 0 else "DOWN"
        elements.append({
            'type': 'strand',
            'name': s_name,
            'start': s_start,
            'end': s_end,
            'direction': direction,
            'index': i,
            'strand_idx': idx
        })
    
    # Сбор элементов - ВСЕ СПИРАЛИ (КАЖДАЯ УНИКАЛЬНАЯ!)
    unique_helices_2d = {}
    helix_count = 0
    for h_keys in result['helices']:
        if len(h_keys) < self.MIN_HELIX_LENGTH:
            continue
        h_start = self._get_res_num(h_keys[0])
        h_end = self._get_res_num(h_keys[-1])
        
        if h_start in result['helix_sides']:
            side = result['helix_sides'][h_start]
            dist = result['helix_distances'].get(h_start, 0)
            num_part = self._get_helix_number(h_start, h_end, path_map)
            name = self._get_helix_name(side, h_start, num_part)
            if num_part:
                display_name = f"{side}{num_part}_{h_start}"  # Hd1_153
            else:
                display_name = f"{side}_{h_start}"            # Hd_153
            contacts = []
            for i, idx in enumerate(full_path):
                s_coords = np.array([self.res_data[key]['coords'] for key in result['strands'][idx]])
                h_c = np.array([self.res_data[r]['coords'] for r in h_keys])
                if np.min(distance_matrix(s_coords, h_c)) < self.HELIX_RADIUS:
                    contacts.append(idx)
            
            # УНИКАЛЬНЫЙ КЛЮЧ: имя + диапазон
            helix_key = f"{name}_{h_start}_{h_end}"
            if helix_key not in unique_helices_2d:
                unique_helices_2d[helix_key] = (display_name, side, h_start, h_end, contacts, dist)
                helix_count += 1

    print(f"Найдено {helix_count} спиралей для отображения")
    
    for helix_key, (display_name, side, h_start, h_end, contacts, dist) in unique_helices_2d.items():
        elements.append({
            'type': 'helix',
            'name': display_name,
            'start': h_start,
            'end': h_end,
            'side': side,
            'strand_contacts': contacts,
            'distance': dist
        })
    
    elements.sort(key=lambda x: x['start'])
    
    print(f"Всего элементов для отображения: {len(elements)}")
    
    # =================================================================
    # ПОЗИЦИОНИРОВАНИЕ ПО ЦЕНТРУ МАСС
    # =================================================================
    pos_map = {}
    strand_x_positions = {}
    strands = [e for e in elements if e['type'] == 'strand']
    
    # Получаем центры тяжей
    strand_centers = {}
    for strand in strands:
        idx = strand['strand_idx']
        s_range = result['strands'][idx]
        coords = np.array([self.res_data[key]['coords'] for key in s_range])
        strand_centers[idx] = coords.mean(axis=0)
    
    strands_sorted = sorted(strands, key=lambda x: x['start'])
    
    # Вычисляем проекцию на ось цепи
    if len(strands_sorted) > 1:
        first_center = strand_centers[strands_sorted[0]['strand_idx']]
        last_center = strand_centers[strands_sorted[-1]['strand_idx']]
        chain_vector = last_center - first_center
        chain_length = np.linalg.norm(chain_vector)
        if chain_length > 0:
            chain_direction = chain_vector / chain_length
        else:
            chain_direction = np.array([1, 0, 0])
    else:
        chain_direction = np.array([1, 0, 0])
        first_center = strand_centers[strands_sorted[0]['strand_idx']] if strands_sorted else np.array([0, 0, 0])
    
    # Проецируем центры на направление цепи
    projections = []
    for strand in strands_sorted:
        center = strand_centers[strand['strand_idx']]
        vec = center - first_center
        proj = np.dot(vec, chain_direction)
        projections.append(proj)
    
    # Нормируем для 2D
    min_proj = min(projections) if projections else 0
    max_proj = max(projections) if projections else 1
    scale_factor = 15.0
    
    # НЕ ЗЕРКАЛИМ ЗДЕСЬ!
    for i, strand in enumerate(strands_sorted):
        if max_proj > min_proj:
            norm_proj = (projections[i] - min_proj) / (max_proj - min_proj) * scale_factor
            # norm_proj = scale_factor - norm_proj  # НЕ ЗЕРКАЛИМ!
        else:
            norm_proj = (len(strands_sorted) - 1 - i) * 1.5
            
        original_dir = strand['direction']
        pos_map[strand['name']] = {
            'x': norm_proj,
            'y': 1,
            'original_direction': original_dir,
            'start': strand['start'],
            'end': strand['end'],
            'type': 'strand',
            'strand_idx': strand['strand_idx']
        }
        strand_x_positions[strand['strand_idx']] = norm_proj
    
    # =================================================================
    # ФУНКЦИЯ ДЛЯ ПРОВЕРКИ ПЕРЕСЕЧЕНИЙ
    # =================================================================
    def elements_intersect(pos1, pos2, threshold=1.0):
        """Проверяет, пересекаются ли два элемента"""
        dist = np.sqrt((pos1['x'] - pos2['x'])**2 + (pos1['y'] - pos2['y'])**2)
        return dist < threshold
    
    # =================================================================
    # РАЗНОСИМ ТЯЖИ, СОХРАНЯЯ ПРАВИЛЬНЫЙ ПОРЯДОК (ОТ БОЛЬШЕГО НОМЕРА К МЕНЬШЕМУ)
    # =================================================================
    # Создаем список тяжей
    strand_positions = []
    for name, data in pos_map.items():
        if data['type'] == 'strand':
            # Извлекаем номер из имени (S7, S6, S5, S4, S3, S2, S1, S0, S-1)
            if name.startswith('S-'):
                num = -int(name[2:])  # S-1 -> -1
            else:
                num = int(name[1:])    # S7 -> 7, S0 -> 0
            
            strand_positions.append({
                'name': name,
                'x': data['x'],
                'y': data['y'],
                'num': num,
                'data': data
            })
    
    # СОРТИРУЕМ ПО УБЫВАНИЮ ЧИСЛОВОГО НОМЕРА (ОТ БОЛЬШЕГО К МЕНЬШЕМУ)
    # S7 (7) -> S6 (6) -> S5 (5) -> S4 (4) -> S3 (3) -> S2 (2) -> S1 (1) -> S0 (0) -> S-1 (-1)
    strand_positions.sort(key=lambda x: x['num'], reverse=True)
    
    print(f"Правильный порядок тяжей (слева направо): {[s['name'] for s in strand_positions]}")
    
    # Разносим пересекающиеся тяжи, НО НЕ МЕНЯЕМ ПОРЯДОК!
    moved = True
    max_iterations = 50
    iteration = 0
    min_distance = 1.5
    
    while moved and iteration < max_iterations:
        moved = False
        iteration += 1
        
        # Проверяем только СОСЕДНИЕ тяжи в правильном порядке
        for i in range(len(strand_positions) - 1):
            current = strand_positions[i]
            next_strand = strand_positions[i + 1]
            
            dist = next_strand['x'] - current['x']
            
            if dist < min_distance:
                moved = True
                overlap = min_distance - dist
                current['x'] -= overlap / 2
                next_strand['x'] += overlap / 2
    
    # Обновляем позиции тяжей
    for sp in strand_positions:
        pos_map[sp['name']]['x'] = sp['x']
        strand_x_positions[pos_map[sp['name']]['strand_idx']] = sp['x']
    
    # =================================================================
    # ЗЕРКАЛЬНОЕ ОТРАЖЕНИЕ - ПОСЛЕ РАЗНЕСЕНИЯ!
    # =================================================================
    all_x = [data['x'] for data in pos_map.values() if data['type'] == 'strand']
    if all_x:
        min_x = min(all_x)
        max_x = max(all_x)
        
        # Зеркально отражаем (чтобы S1 был слева)
        for name, data in pos_map.items():
            if data['type'] == 'strand':
                data['x'] = max_x + min_x - data['x']
        
        # Обновляем strand_x_positions после отражения
        for name, data in pos_map.items():
            if data['type'] == 'strand':
                strand_x_positions[data['strand_idx']] = data['x']
    
    # =================================================================
    # ПОЗИЦИОНИРОВАНИЕ СПИРАЛЕЙ С ИСКУССТВЕННЫМ РАЗВЕДЕНИЕМ
    # =================================================================
    helices = [e for e in elements if e['type'] == 'helix']
    print(f"Позиционирование {len(helices)} спиралей")
    
    # Сначала размещаем все спирали на базовых позициях
    helix_positions = []
    for helix in helices:
        # Находим базовую X координату
        min_dist = float('inf')
        base_x = 0
        for strand_name, s_data in pos_map.items():
            if s_data['type'] == 'strand':
                dist = abs(helix['start'] - s_data['start'])
                if dist < min_dist:
                    min_dist = dist
                    base_x = s_data['x']
        
        # Если есть контакты, усредняем
        if 'strand_contacts' in helix and helix['strand_contacts']:
            contact_xs = [strand_x_positions[idx] for idx in helix['strand_contacts']
                        if idx in strand_x_positions]
            if contact_xs:
                base_x = np.mean(contact_xs)
        
        base_y = 2.2 if helix['side'] == 'Hu' else -0.2
        
        helix_positions.append({
            'name': helix['name'],
            'helix': helix,
            'x': base_x,
            'y': base_y,
            'side': helix['side']
        })
    
    # Сортируем по X
    helix_positions.sort(key=lambda x: x['x'])
    
    # Разносим пересекающиеся спирали
    moved = True
    iteration = 0
    
    while moved and iteration < max_iterations:
        moved = False
        iteration += 1
        
        # Проверяем пересечения спиралей между собой
        for i in range(len(helix_positions)):
            for j in range(i+1, len(helix_positions)):
                pos_i = {'x': helix_positions[i]['x'], 'y': helix_positions[i]['y']}
                pos_j = {'x': helix_positions[j]['x'], 'y': helix_positions[j]['y']}
                
                if elements_intersect(pos_i, pos_j, 1.2):
                    moved = True
                    # Раздвигаем по горизонтали
                    if i < j:
                        helix_positions[i]['x'] = min(helix_positions[i]['x'], 
                                                      helix_positions[j]['x'] - 1.5)
                        helix_positions[j]['x'] = max(helix_positions[j]['x'], 
                                                      helix_positions[i]['x'] + 1.5)
                    else:
                        helix_positions[j]['x'] = min(helix_positions[j]['x'], 
                                                      helix_positions[i]['x'] - 1.5)
                        helix_positions[i]['x'] = max(helix_positions[i]['x'], 
                                                      helix_positions[j]['x'] + 1.5)
        
        # Проверяем пересечения спиралей с тяжами
        for h in helix_positions:
            for s in strand_positions:
                pos_h = {'x': h['x'], 'y': h['y']}
                pos_s = {'x': s['x'], 'y': s['y']}
                
                if elements_intersect(pos_h, pos_s, 1.2):
                    moved = True
                    # Сдвигаем спираль по вертикали
                    if h['y'] > 1:  # Hu спираль
                        h['y'] += 0.3
                    else:  # Hd спираль
                        h['y'] -= 0.3
    
    # Сохраняем позиции спиралей
    for hp in helix_positions:
        pos_map[hp['name']] = {
            'x': hp['x'],
            'y': hp['y'],
            'start': hp['helix']['start'],
            'end': hp['helix']['end'],
            'distance': hp['helix'].get('distance', 0),
            'side': hp['side'],
            'type': 'helix'
        }
    
    # =================================================================
    # ОПРЕДЕЛЯЕМ N И C КОНЦЫ
    # =================================================================
    if elements:
        n_element = elements[0]  # первый по последовательности
        c_element = elements[-1]  # последний по последовательности
    else:
        n_element = None
        c_element = None
    
    # =================================================================
    # ПОСТРОЕНИЕ УМНЫХ СОЕДИНЕНИЙ
    # =================================================================
    sequence_order = []
    for elem in sorted(elements, key=lambda x: x['start']):
        sequence_order.append(elem['name'])
    
    print(f"Порядок следования: {sequence_order}")
    
    # Функция для подсчета пересечений линии с элементами
    def count_line_intersections(x1, y1, x2, y2, excluded_names):
        count = 0
        for t in np.linspace(0, 1, 30):
            x = x1 + t * (x2 - x1)
            y = y1 + t * (y2 - y1)
            for name, data in pos_map.items():
                if name in excluded_names:
                    continue
                dist = np.sqrt((x - data['x'])**2 + (y - data['y'])**2)
                if dist < 0.8:
                    count += 1
                    break
        return count
    
    # Функция для подсчета пересечений дуги с элементами
    def count_arc_intersections(x1, y1, x2, y2, cy, excluded_names):
        count = 0
        for t in np.linspace(0, 1, 40):
            x = x1 + t * (x2 - x1)
            y = y1 + t * (y2 - y1) + 4 * (cy - (y1+y2)/2) * t * (1 - t)
            for name, data in pos_map.items():
                if name in excluded_names:
                    continue
                dist = np.sqrt((x - data['x'])**2 + (y - data['y'])**2)
                if dist < 0.8:
                    count += 1
                    break
        return count
    
    # Собираем информацию о расстояниях между элементами
    element_starts = {}
    for elem in elements:
        element_starts[elem['name']] = elem['start']
    
    # Создаем соединения для ВСЕХ последовательных элементов
    connections = []
    
    for i in range(len(sequence_order)-1):
        name1, name2 = sequence_order[i], sequence_order[i+1]
        if name1 in pos_map and name2 in pos_map:
            x1, y1 = pos_map[name1]['x'], pos_map[name1]['y']
            x2, y2 = pos_map[name2]['x'], pos_map[name2]['y']
            
            # Вычисляем расстояние в аминокислотах
            start1 = element_starts[name1]
            start2 = element_starts[name2]
            aa_distance = abs(start2 - start1)
            
            # Определяем тип соединения
            is_large_insertion = aa_distance > 50
            base_arc_height = 1.0 if is_large_insertion else 0.6
            
            # Сначала пробуем прямую линию
            line_intersections = count_line_intersections(x1, y1, x2, y2, [name1, name2])
            
            if line_intersections == 0:
                # Прямая линия идеальна
                connections.append({
                    'start': name1,
                    'end': name2,
                    'type': 'line',
                    'params': (x1, y1, x2, y2),
                    'is_large': is_large_insertion
                })
                continue
            
            # Пробуем дугу вверх с увеличивающейся высотой
            best_up_intersections = float('inf')
            best_up_cy = None
            
            for height_mult in [1.0, 1.5, 2.0, 2.5, 3.0]:
                test_height = base_arc_height * height_mult
                test_cy = max(y1, y2) + test_height
                intersections = count_arc_intersections(x1, y1, x2, y2, test_cy, [name1, name2])
                
                if intersections < best_up_intersections:
                    best_up_intersections = intersections
                    best_up_cy = test_cy
                
                if intersections == 0:
                    break
            
            # Пробуем дугу вниз с увеличивающейся глубиной
            best_down_intersections = float('inf')
            best_down_cy = None
            
            for height_mult in [1.0, 1.5, 2.0, 2.5, 3.0]:
                test_height = base_arc_height * height_mult
                test_cy = min(y1, y2) - test_height
                intersections = count_arc_intersections(x1, y1, x2, y2, test_cy, [name1, name2])
                
                if intersections < best_down_intersections:
                    best_down_intersections = intersections
                    best_down_cy = test_cy
                
                if intersections == 0:
                    break
            
            # Выбираем лучший вариант
            if best_up_intersections == 0:
                connections.append({
                    'start': name1,
                    'end': name2,
                    'type': 'arc',
                    'params': (x1, y1, x2, y2, best_up_cy),
                    'is_large': is_large_insertion
                })
            elif best_down_intersections == 0:
                connections.append({
                    'start': name1,
                    'end': name2,
                    'type': 'arc',
                    'params': (x1, y1, x2, y2, best_down_cy),
                    'is_large': is_large_insertion
                })
            else:
                # Выбираем вариант с наименьшим количеством пересечений
                if best_up_intersections <= best_down_intersections:
                    connections.append({
                        'start': name1,
                        'end': name2,
                        'type': 'arc',
                        'params': (x1, y1, x2, y2, best_up_cy),
                        'is_large': is_large_insertion
                    })
                else:
                    connections.append({
                        'start': name1,
                        'end': name2,
                        'type': 'arc',
                        'params': (x1, y1, x2, y2, best_down_cy),
                        'is_large': is_large_insertion
                    })
    
    print(f"Создано {len(connections)} соединений")
    
    # =================================================================
    # СОЗДАНИЕ PLOTLY ФИГУРЫ
    # =================================================================
    fig = go.Figure()
    
    # Добавляем ВСЕ тяжи
    for name, data in pos_map.items():
        if data['type'] == 'strand':
            symbol = 'triangle-up' if data['original_direction'] == 'UP' else 'triangle-down'
            
            # Определяем цвет (золотой для N/C)
            if n_element and n_element['name'] == name:
                marker_color = '#FFD700'  # золотой для N-конца
                text = 'N'
            elif c_element and c_element['name'] == name:
                marker_color = '#FFA500'  # оранжевый для C-конца
                text = 'C'
            else:
                marker_color = '#1a5276'
                text = ''
            
            fig.add_trace(go.Scatter(
                x=[data['x']],
                y=[data['y']],
                mode='markers+text',
                marker=dict(
                    symbol=symbol,
                    size=45,
                    color=marker_color,
                    line=dict(color='white', width=2)
                ),
                text=[text],
                textposition="middle center",
                textfont=dict(size=14, color='white', family='Arial Black'),
                name=name,
                showlegend=False,
                hovertemplate=f"<b>{name}</b><br>Range: {data['start']}-{data['end']}<br>Direction: {data['original_direction']}<extra></extra>"
            ))
    
    # Добавляем ВСЕ спирали
    for name, data in pos_map.items():
        if data['type'] == 'helix':
            # Определяем цвет
            if n_element and n_element['name'] == name:
                marker_color = '#FFD700'  # золотой для N-конца
                text = 'N'
            elif c_element and c_element['name'] == name:
                marker_color = '#FFA500'  # оранжевый для C-конца
                text = 'C'
            else:
                if data['side'] == 'Hu':
                    marker_color = '#27ae60'  # зеленый
                else:
                    marker_color = '#e74c3c'  # красный
                text = ''
            
            fig.add_trace(go.Scatter(
                x=[data['x']],
                y=[data['y']],
                mode='markers+text',
                marker=dict(
                    symbol='circle',
                    size=40,
                    color=marker_color,
                    line=dict(color='white', width=2)
                ),
                text=[text],
                textposition="middle center",
                textfont=dict(size=14, color='white', family='Arial Black'),
                name=name,
                showlegend=False,
                hovertemplate=f"<b>{name}</b><br>Range: {data['start']}-{data['end']}<br>Distance: {data.get('distance', 0):.1f}Å<extra></extra>"
            ))
    
    # Добавляем соединения
    for conn in connections:
        if conn['type'] == 'line':
            x1, y1, x2, y2 = conn['params']
            line_color = '#ff9900' if conn['is_large'] else '#8e44ad'
            line_width = 4 if conn['is_large'] else 3
            line_dash = 'dot' if conn['is_large'] else 'solid'
            
            fig.add_trace(go.Scatter(
                x=[x1, x2],
                y=[y1, y2],
                mode='lines',
                line=dict(color=line_color, width=line_width, dash=line_dash),
                showlegend=False,
                hoverinfo='none'
            ))
            
        elif conn['type'] == 'arc':
            x1, y1, x2, y2, cy = conn['params']
            
            t = np.linspace(0, 1, 30)
            x = x1 + t * (x2 - x1)
            y = y1 + t * (y2 - y1) + 0.5 * (cy - (y1+y2)/2) * np.sin(t * np.pi)
            
            line_color = '#ff9900' if conn['is_large'] else '#8e44ad'
            line_width = 4 if conn['is_large'] else 3
            line_dash = 'dot' if conn['is_large'] else 'solid'
            
            fig.add_trace(go.Scatter(
                x=x, y=y,
                mode='lines',
                line=dict(color=line_color, width=line_width, dash=line_dash),
                showlegend=False,
                hoverinfo='none'
            ))
    
    # =================================================================
    # ДОБАВЛЯЕМ ЛЕГЕНДУ
    # =================================================================
    # Добавляем невидимые трейсы для легенды
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='#8e44ad', width=3),
        name='N→C chain (normal)',
        showlegend=True
    ))
    
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='#ff9900', width=4, dash='dot'),
        name='N→C chain (large insertion >50aa)',
        showlegend=True
    ))
    
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='triangle-up', size=15, color='#1a5276'),
        name='Beta-strand UP',
        showlegend=True
    ))
    
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='triangle-down', size=15, color='#1a5276'),
        name='Beta-strand DOWN',
        showlegend=True
    ))
    
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='circle', size=15, color='#27ae60'),
        name='Helix Hu (Up)',
        showlegend=True
    ))
    
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='circle', size=15, color='#e74c3c'),
        name='Helix Hd (Down)',
        showlegend=True
    ))
    
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='triangle-up', size=15, color='#FFD700', line=dict(color='black', width=1)),
        name='N-terminus',
        showlegend=True
    ))
    
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='triangle-down', size=15, color='#FFA500', line=dict(color='black', width=1)),
        name='C-terminus',
        showlegend=True
    ))
    
    # Настройка макета
    all_x = [data['x'] for data in pos_map.values()]
    max_x = max(all_x) if all_x else 10
    min_x = min(all_x) if all_x else 0
    
    fig.update_layout(
        title=f"2D Protein Topology: Chain {result['chain']}<br>Motif: {self.motif_info['text']} ({self.motif_info['res']}) | S4: {path_map['S4'][0]}-{path_map['S4'][1]}",
        title_font_size=18,
        width=1400,
        height=900,
        showlegend=True,
        hovermode='closest',
        plot_bgcolor='white',
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            visible=False,
            range=[min_x - 2, max_x + 2]
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            visible=False,
            range=[-3, 6]
        ),
        legend=dict(
            groupclick="toggleitem",
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99,
            bgcolor='rgba(255,255,255,0.8)',
            font=dict(size=12)
        )
    )
    
    # Добавляем очень бледные опорные линии
    fig.add_hline(y=1, line_dash="dot", line_color="lightgray", opacity=0.1)
    fig.add_hline(y=2.2, line_dash="dot", line_color="lightgray", opacity=0.1)
    fig.add_hline(y=-0.2, line_dash="dot", line_color="lightgray", opacity=0.1)
    
    print("Фигура создана успешно")
    return fig


def _merge_display_helices(self, elements, result):
    """Объединение соседних спиралей для отображения"""
    return elements


# Прикрепляем методы к классу
MTaseAnalyzer.visualize_topology_interactive = visualize_topology_interactive
MTaseAnalyzer._merge_display_helices = _merge_display_helices