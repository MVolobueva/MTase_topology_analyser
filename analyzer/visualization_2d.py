import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import numpy as np
from .core import MTaseAnalyzer

def visualize_topology_from_analysis(self, result):
    """2D визуализация топологии"""
    if not result:
        return

    full_path = result['full_path']
    path_map = result['path_map']
    strand_names = result['strand_names']
    s4_idx = result['s4_idx']
    s4_pos = full_path.index(s4_idx)
    motif_res = self.motif_info['res']
    motif_chain = result['chain']
    motif_text = self.motif_info['text']

    elements = []
    v4 = result['v4']

    # -----------------------------------------------------------------
    # 1. ДОБАВЛЯЕМ ТЯЖИ
    # -----------------------------------------------------------------
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

    # -----------------------------------------------------------------
    # 2. ДОБАВЛЯЕМ СПИРАЛИ
    # -----------------------------------------------------------------
    unique_helices_2d = {}
    for h_keys in result['helices']:
        if len(h_keys) < self.MIN_HELIX_LENGTH:
            continue
            
        h_start = self._get_res_num(h_keys[0])
        h_end = self._get_res_num(h_keys[-1])
        helix_key = f"{h_start}-{h_end}"

        if h_start in result['helix_sides'] and helix_key not in unique_helices_2d:
            side = result['helix_sides'][h_start]
            dist = result['helix_distances'][h_start]
            num_part = self._get_helix_number(h_start, h_end, path_map)
            name = self._get_helix_name(side, h_start, num_part)

            # Находим все тяжи в контакте со спиралью
            contacts = []
            for i, idx in enumerate(full_path):
                s_coords = np.array([self.res_data[key]['coords'] 
                                    for key in result['strands'][idx]])
                h_c = np.array([self.res_data[r]['coords'] for r in h_keys])
                if np.min(distance_matrix(s_coords, h_c)) < self.HELIX_RADIUS:
                    contacts.append(idx)

            unique_helices_2d[helix_key] = (name, side, h_start, h_end, contacts, dist)

    for helix_key, (name, side, h_start, h_end, contacts, dist) in unique_helices_2d.items():
        elements.append({
            'type': 'helix',
            'name': name,
            'start': h_start,
            'end': h_end,
            'side': side,
            'strand_contacts': contacts,
            'distance': dist
        })

    # -----------------------------------------------------------------
    # 3. ОБЪЕДИНЕНИЕ СОСЕДНИХ СПИРАЛЕЙ
    # -----------------------------------------------------------------
    merged_elements = []
    elements.sort(key=lambda x: x['start'])
    i = 0
    while i < len(elements):
        if elements[i]['type'] == 'helix':
            j = i
            while j + 1 < len(elements) and elements[j+1]['type'] == 'helix':
                if elements[j]['name'] == elements[j+1]['name']:
                    # Проверяем 3D расстояние между спиралями
                    dist_3d = float('inf')
                    for h_keys in result['helices']:
                        h1_start = self._get_res_num(h_keys[0])
                        h1_end = self._get_res_num(h_keys[-1])
                        if h1_start == elements[j]['start'] and h1_end == elements[j]['end']:
                            for h2_keys in result['helices']:
                                h2_start = self._get_res_num(h2_keys[0])
                                h2_end = self._get_res_num(h2_keys[-1])
                                if h2_start == elements[j+1]['start'] and h2_end == elements[j+1]['end']:
                                    coord1 = self.res_data[h_keys[-1]]['coords']
                                    coord2 = self.res_data[h2_keys[0]]['coords']
                                    dist_3d = np.linalg.norm(coord1 - coord2)
                                    break
                            break
                    if dist_3d <= 5.0:
                        j += 1
                    else:
                        break
                else:
                    break
                    
            if j > i:
                merged_helix = elements[i].copy()
                merged_helix['start'] = min(e['start'] for e in elements[i:j+1])
                merged_helix['end'] = max(e['end'] for e in elements[i:j+1])
                merged_helix['merged'] = True
                all_contacts = []
                all_distances = []
                for e in elements[i:j+1]:
                    if 'strand_contacts' in e:
                        all_contacts.extend(e['strand_contacts'])
                    if 'distance' in e:
                        all_distances.append(e['distance'])
                merged_helix['strand_contacts'] = list(set(all_contacts))
                merged_helix['distance'] = min(all_distances) if all_distances else 0
                merged_elements.append(merged_helix)
                i = j + 1
            else:
                merged_elements.append(elements[i])
                i += 1
        else:
            merged_elements.append(elements[i])
            i += 1
    elements = merged_elements
    elements.sort(key=lambda x: x['start'])

    # -----------------------------------------------------------------
    # 4. ПОЗИЦИОНИРОВАНИЕ ТЯЖЕЙ
    # -----------------------------------------------------------------
    pos_map = {}
    strand_x_positions = {}
    strands = [e for e in elements if e['type'] == 'strand']

    def get_strand_sort_key(strand):
        name = strand['name']
        if name.startswith('S-'):
            return -int(name[2:])
        else:
            return int(name[1:])

    strands_sorted = sorted(strands, key=get_strand_sort_key)

    x_pos = 0
    for strand in strands_sorted:
        original_dir = strand['direction']
        graph_dir = "DOWN" if original_dir == "UP" else "UP"
        pos_map[strand['name']] = {
            'x': x_pos,
            'y': 1,
            'direction': graph_dir,
            'original_direction': original_dir,
            'start': strand['start'],
            'end': strand['end']
        }
        strand_x_positions[strand['strand_idx']] = x_pos
        x_pos += 1

    # -----------------------------------------------------------------
    # 5. ПОЗИЦИОНИРОВАНИЕ СПИРАЛЕЙ
    # -----------------------------------------------------------------
    helices = [e for e in elements if e['type'] == 'helix']
    helices.sort(key=lambda h: h['start'])
    helix_groups = {}

    for helix in helices:
        if 'strand_contacts' in helix and helix['strand_contacts']:
            contact_xs = [strand_x_positions[idx] for idx in helix['strand_contacts']
                        if idx in strand_x_positions]
            base_x = np.mean(contact_xs) if contact_xs else 0
        else:
            base_x = 0

        y = 2 if helix['side'] == 'Hu' else 0
        group_key = (round(base_x * 2) / 2, y)

        if group_key not in helix_groups:
            helix_groups[group_key] = []
        helix_groups[group_key].append((helix, base_x))

    for (group_x, group_y), helix_list in helix_groups.items():
        if len(helix_list) == 1:
            helix, base_x = helix_list[0]
            final_x = base_x
            pos_map[helix['name']] = {
                'x': final_x,
                'y': group_y,
                'direction': None,
                'start': helix['start'],
                'end': helix['end'],
                'distance': helix.get('distance', 0)
            }
        else:
            total_width = 0.8 * (len(helix_list) - 1)
            start_x = group_x - total_width / 2
            for idx, (helix, base_x) in enumerate(helix_list):
                final_x = start_x + idx * 0.8
                pos_map[helix['name']] = {
                    'x': final_x,
                    'y': group_y,
                    'direction': None,
                    'start': helix['start'],
                    'end': helix['end'],
                    'distance': helix.get('distance', 0)
                }

    # -----------------------------------------------------------------
    # 6. СВЯЗИ МЕЖДУ СПИРАЛЯМИ И ТЯЖАМИ
    # -----------------------------------------------------------------
    connections = []
    for elem in elements:
        if elem['type'] == 'helix' and 'strand_contacts' in elem:
            for strand_idx in elem['strand_contacts']:
                for strand_elem in elements:
                    if strand_elem['type'] == 'strand' and strand_elem['strand_idx'] == strand_idx:
                        connections.append((elem['name'], strand_elem['name']))
                        break

    # -----------------------------------------------------------------
    # 7. ПОРЯДОК В ЦЕПИ
    # -----------------------------------------------------------------
    sequence_order = []
    for elem in sorted(elements, key=lambda x: x['start']):
        sequence_order.append(elem['name'])

    # -----------------------------------------------------------------
    # 8. РИСОВАНИЕ
    # -----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(16, 7))
    
    strand_color = self.COLORS['strand']
    hu_color = self.COLORS['Hu']
    hd_color = self.COLORS['Hd']
    connection_color = '#7f8c8d'
    main_path_color = '#8e44ad'
    text_positions = {}

    # Рисуем тяжи
    for name, data in pos_map.items():
        x, y = data['x'], data['y']
        if name.startswith("S"):
            marker = '^' if data['direction'] == "UP" else 'v'
            ax.plot(x, y, marker, markersize=40, color=strand_color,
                    mfc='white', mew=3.5, zorder=5, alpha=0.9)
            
            dir_symbol = "↑" if data['original_direction'] == "UP" else "↓"
            label = f"{name}{dir_symbol}"
            
            text_offset = 0.7 if y >= 1 else -0.7
            key = (round(x, 1), round(y + text_offset, 1))
            if key in text_positions:
                text_offset *= 1.3
            text_positions[key] = name
            
            ax.text(x, y - text_offset, label, ha='center', fontsize=14, fontweight='bold',
                    color=strand_color, bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                                edgecolor=strand_color, alpha=0.8))
            
            ax.text(x, y - text_offset - 0.4, f"{data['start']}-{data['end']}", 
                    ha='center', fontsize=9, color='black', style='italic', alpha=0.8)

    # Рисуем спирали
    for name, data in pos_map.items():
        if not name.startswith("S"):
            x, y = data['x'], data['y']
            color = hu_color if "Hu" in name else hd_color
            
            ax.plot(x, y, 'o', markersize=35, color=color,
                    mfc=color, mew=2.5, zorder=5, alpha=0.9)
            
            dist_text = f"{data.get('distance', 0):.1f}Å"
            display_name = name
            
            text_offset = -0.9 if "Hd" in name else 0.9
            key = (round(x, 1), round(y + text_offset, 1))
            
            attempts = 0
            while key in text_positions and attempts < 3:
                text_offset *= 1.2
                key = (round(x, 1), round(y + text_offset, 1))
                attempts += 1
            text_positions[key] = name
            
            ax.text(x, y + text_offset, display_name, ha='center', fontsize=12, fontweight='bold',
                    color='white', bbox=dict(boxstyle="round,pad=0.3", facecolor=color,
                                            edgecolor='white', alpha=0.9))
            
            ax.text(x, y + text_offset - 0.4, dist_text, ha='center', fontsize=8,
                    color='white', fontweight='bold', alpha=0.9,
                    bbox=dict(boxstyle="round,pad=0.1", facecolor=color, alpha=0.7))
            
            ax.text(x, y + text_offset - 0.8, f"{data['start']}-{data['end']}", 
                    ha='center', fontsize=8, color='black', style='italic', alpha=0.8)

    # Рисуем связи спираль-тяж
    for start, end in connections:
        if start in pos_map and end in pos_map:
            x_vals = [pos_map[start]['x'], pos_map[end]['x']]
            y_vals = [pos_map[start]['y'], pos_map[end]['y']]
            if x_vals[0] > x_vals[1]:
                x_vals = [x_vals[1], x_vals[0]]
                y_vals = [y_vals[1], y_vals[0]]
            
            ax.plot(x_vals, y_vals, color=connection_color,
                    linestyle='--', linewidth=2.0, alpha=0.7, zorder=3)

    # Рисуем путь N→C
    for i in range(len(sequence_order)-1):
        name1, name2 = sequence_order[i], sequence_order[i+1]
        if name1 in pos_map and name2 in pos_map:
            x1, y1 = pos_map[name1]['x'], pos_map[name1]['y']
            x2, y2 = pos_map[name2]['x'], pos_map[name2]['y']
            
            dx = x2 - x1
            arc_height = 0.8
            
            # Проверяем, не пересекает ли дуга другие элементы
            for other_name, other_data in pos_map.items():
                if other_name not in [name1, name2]:
                    ox, oy = other_data['x'], other_data['y']
                    if min(x1, x2) < ox < max(x1, x2) and abs(oy - (y1 + y2)/2) < 0.5:
                        arc_height = 1.2
                        break
            
            if abs(dx) > 1.5:
                control_x = (x1 + x2) / 2
                control_y = max(y1, y2) + arc_height
                
                verts = [(x1, y1), (control_x, control_y), (x2, y2)]
                codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
                path = Path(verts, codes)
                patch = patches.PathPatch(path, facecolor='none', edgecolor=main_path_color,
                                        linewidth=3.0, alpha=0.9, zorder=4)
                ax.add_patch(patch)
            else:
                ax.plot([x1, x2], [y1, y2], color=main_path_color,
                        linewidth=3.0, alpha=0.9, zorder=4)

    # Метки N и C
    if sequence_order:
        first_elem = sequence_order[0]
        if first_elem in pos_map:
            x_n, y_n = pos_map[first_elem]['x'], pos_map[first_elem]['y']
            ax.text(x_n - 0.4, y_n + 0.4, 'N', fontsize=16, fontweight='bold', color='#2c3e50',
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="yellow", alpha=1.0, edgecolor='none'),
                    zorder=20, ha='center', va='center')
            ax.text(x_n - 0.4, y_n + 0.1, f"{pos_map[first_elem]['start']}",
                    fontsize=8, color='#2c3e50', ha='center', va='center', alpha=0.7)
        
        last_elem = sequence_order[-1]
        if last_elem in pos_map:
            x_c, y_c = pos_map[last_elem]['x'], pos_map[last_elem]['y']
            ax.text(x_c + 0.4, y_c + 0.4, 'C', fontsize=16, fontweight='bold', color='#2c3e50',
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="yellow", alpha=1.0, edgecolor='none'),
                    zorder=20, ha='center', va='center')
            ax.text(x_c + 0.4, y_c + 0.1, f"{pos_map[last_elem]['end']}",
                    fontsize=8, color='#2c3e50', ha='center', va='center', alpha=0.7)

    # Легенда
    legend_elements = [
        plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='white',
                    markeredgecolor=strand_color, markersize=15, markeredgewidth=2,
                    label='Beta-strand (UP)'),
        plt.Line2D([0], [0], marker='v', color='w', markerfacecolor='white',
                    markeredgecolor=strand_color, markersize=15, markeredgewidth=2,
                    label='Beta-strand (DOWN)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=hu_color,
                    markersize=15, label='Helix (Hu) - Up'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=hd_color,
                    markersize=15, label='Helix (Hd) - Down'),
        plt.Line2D([0], [0], color=main_path_color, lw=3, label='N→C chain'),
        plt.Line2D([0], [0], color=connection_color, linestyle='--',
                    lw=2, label='Spatial contacts')
    ]

    ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.9)
    
    # Настройки графика
    ax.set_xlim(-1.5, x_pos + 0.5)
    ax.set_ylim(-1.5, 4.5)
    ax.set_axis_off()
    
    # Заголовок
    title_text = f"2D Protein Topology: Chain {motif_chain}\n"
    title_text += f"Motif: {motif_text} ({motif_res}) | S4: {path_map['S4'][0]}-{path_map['S4'][1]}"
    plt.title(title_text, fontsize=16, fontweight='bold', pad=20)
    
    # Опорные линии
    ax.axhline(y=1, color='gray', linestyle=':', alpha=0.3, zorder=0)
    ax.axhline(y=2, color='gray', linestyle=':', alpha=0.3, zorder=0)
    ax.axhline(y=0, color='gray', linestyle=':', alpha=0.3, zorder=0)
    
    plt.tight_layout()
    plt.show()

MTaseAnalyzer.visualize_topology_from_analysis = visualize_topology_from_analysis