from .core import MTaseAnalyzer

def analyze_topology(self, motif_data=None):
    """АНАЛИЗ ТОПОЛОГИИ - ТВОЙ ПРОВЕРЕННЫЙ АЛГОРИТМ"""
    if motif_data:
        self.motif_info = motif_data
    if not self.motif_info:
        print("Ошибка: каталитический мотив не найден")
        return None
    s3_idx = None
    s5_idx = None
    # Очищаем словари спиралей
    self.helix_sides = {}
    self.helix_distances = {}
    self.helix_nearest_strand = {}

    s4_global_idx = self.motif_info['s4_idx']
    motif_key = self.motif_info.get('key', '')
    motif_chain = self.motif_info.get('chain', 'A')
    s4_local_idx = self.motif_info.get('s4_local_idx', None)
    motif_res = self.motif_info['res']
    motif_text = self.motif_info['text']

    print(f"\n{'='*60}")
    print(f"🔬 АНАЛИЗ ТОПОЛОГИИ ДЛЯ ЦЕПИ: {motif_chain}")
    print(f"   Мотив: {motif_text} ({motif_res})")
    print(f"{'='*60}")

    self.current_chain = motif_chain

    # =================================================================
    # ФИЛЬТРАЦИЯ - ТОЛЬКО ЦЕПЬ МОТИВА!
    # =================================================================

    # 1. Фильтруем тяжи
    chain_strands = []
    chain_strand_indices = []
    chain_global_to_local = {}

    for i, strand in enumerate(self.strands):
        if strand and self._get_chain(strand[0]) == motif_chain:
            local_idx = len(chain_strands)
            chain_strands.append(strand)
            chain_strand_indices.append(i)
            chain_global_to_local[i] = local_idx

    print(f"\n📊 ТЯЖИ ЦЕПИ {motif_chain}: {len(chain_strands)}")
    
    if len(chain_strands) == 0:
        print(f"❌ ОШИБКА: В цепи {motif_chain} не найдено тяжей!")
        return None

    # 2. Определяем локальный индекс S4
    if s4_local_idx is not None:
        s4_idx = s4_local_idx
        print(f"   Используем сохраненный локальный индекс S4: {s4_idx}")
    else:
        if s4_global_idx in chain_global_to_local:
            s4_idx = chain_global_to_local[s4_global_idx]
            print(f"   S4 глобальный индекс {s4_global_idx} -> локальный {s4_idx}")
        else:
            print(f"❌ ОШИБКА: S4 (глобальный индекс {s4_global_idx}) не найден в цепи {motif_chain}!")
            print(f"   Доступные глобальные индексы: {chain_strand_indices}")
            return None

    # 3. Сохраняем оригиналы
    original_strands = self.strands
    original_adj = self.adj
    original_helices = self.helices

    # 4. Заменяем self.strands на отфильтрованные
    self.strands = chain_strands

    # 5. Создаем новый adj только для этой цепи
    self.adj = collections.defaultdict(set)
    for i in range(len(self.strands)):
        for j in range(i + 1, len(self.strands)):
            s1_c = np.array([self.res_data[r]['coords'] for r in self.strands[i]])
            s2_c = np.array([self.res_data[r]['coords'] for r in self.strands[j]])
            if np.min(distance_matrix(s1_c, s2_c)) < self.CONTACT_DIST:
                self.adj[i].add(j)
                self.adj[j].add(i)

    # 6. Фильтруем спирали
    self.helices = [h for h in self.helices if h and self._get_chain(h[0]) == motif_chain]
    print(f"   Спиралей в цепи: {len(self.helices)}")


    # =================================================================
    # ОСНОВНОЙ АНАЛИЗ
    # =================================================================

    s4_strand = self.strands[s4_idx]
    s4_start = self._get_res_num(s4_strand[0])
    s4_end = self._get_res_num(s4_strand[-1])
    v4 = self._get_strand_vector(s4_strand)

    v_set = {s4_idx}
    side_down = self._expand_path(s4_idx, v_set)
    side_up = self._expand_path(s4_idx, v_set)
    full_path = side_up[::-1][:-1] + side_down
    s4_pos = full_path.index(s4_idx)

    strand_names = {s4_idx: "S4"}
    path_map = {"S4": (s4_start, s4_end)}
    processed = {s4_idx}

    # -----------------------------------------------------------------
    # РАЗДЕЛЕНИЕ СОСЕДЕЙ НА ЛЕВЫХ И ПРАВЫХ
    # -----------------------------------------------------------------
    left_neighbors = []
    right_neighbors = []

    for neighbor in self.adj[s4_idx]:
        if neighbor in processed:
            continue
        strand = self.strands[neighbor]
        strand_end = self._get_res_num(strand[-1])
        strand_start = self._get_res_num(strand[0])
        if strand_end < s4_start:
            left_neighbors.append((neighbor, strand_end, strand_start))
        elif strand_start > s4_end:
            right_neighbors.append((neighbor, strand_start, strand_end))

    # ------------------------------------------------------------
    # ШАГ 1: ОБРАБАТЫВАЕМ C-КОНЕЦ (правые соседи)
    # ------------------------------------------------------------
    if right_neighbors:
        print(f"\n📌 ОБРАБОТКА C-КОНЦА: {len(right_neighbors)} тяжа(ей)")
        right_neighbors.sort(key=lambda x: x[1])

        # ПЕРВЫЙ = S3
        s3_idx = right_neighbors[0][0]
        strand_names[s3_idx] = "S3"
        path_map["S3"] = (
            self._get_res_num(self.strands[s3_idx][0]),
            self._get_res_num(self.strands[s3_idx][-1])
        )
        processed.add(s3_idx)
        print(f"  → S3: тяж {s3_idx} ({path_map['S3'][0]}-{path_map['S3'][1]})")

        # ВТОРОЙ = S5 (если есть)
        if len(right_neighbors) > 1:
            s5_idx = right_neighbors[1][0]
            strand_names[s5_idx] = "S5"
            path_map["S5"] = (
                self._get_res_num(self.strands[s5_idx][0]),
                self._get_res_num(self.strands[s5_idx][-1])
            )
            processed.add(s5_idx)
            print(f"  → S5: тяж {s5_idx} ({path_map['S5'][0]}-{path_map['S5'][1]})")

            # Идем ВЛЕВО от S5: S6, S7, S8...
            current_idx = s5_idx
            current_number = 6
            while True:
                next_idx = None
                for neighbor in self.adj[current_idx]:
                    if neighbor not in processed:
                        next_idx = neighbor
                        break
                if next_idx is None:
                    break
                if next_idx in full_path:
                    name = f"S{current_number}"
                    strand_names[next_idx] = name
                    path_map[name] = (
                        self._get_res_num(self.strands[next_idx][0]),
                        self._get_res_num(self.strands[next_idx][-1])
                    )
                    processed.add(next_idx)
                    print(f"    → {name}: тяж {next_idx} ({path_map[name][0]}-{path_map[name][1]})")
                    current_idx = next_idx
                    current_number += 1
                else:
                    break

        # Идем ВПРАВО от S3: S2, S1, S0...
        current_idx = s3_idx
        current_number = 2
        while True:
            next_idx = None
            for neighbor in self.adj[current_idx]:
                if neighbor not in processed:
                    next_idx = neighbor
                    break
            if next_idx is None:
                break
            if next_idx in full_path:
                name = f"S{current_number}"
                strand_names[next_idx] = name
                path_map[name] = (
                    self._get_res_num(self.strands[next_idx][0]),
                    self._get_res_num(self.strands[next_idx][-1])
                )
                processed.add(next_idx)
                print(f"    → {name}: тяж {next_idx} ({path_map[name][0]}-{path_map[name][1]})")
                current_idx = next_idx
                current_number -= 1
            else:
                break

    # ------------------------------------------------------------
    # ШАГ 2: ОБРАБАТЫВАЕМ N-КОНЕЦ (левые соседи)
    # ------------------------------------------------------------
    if left_neighbors:
        print(f"\n📌 ОБРАБОТКА N-КОНЦА: {len(left_neighbors)} тяжа(ей)")
        left_neighbors.sort(key=lambda x: x[1], reverse=True)

        # ПЕРВЫЙ = S5
        s5_idx = left_neighbors[0][0]
        strand_names[s5_idx] = "S5"
        path_map["S5"] = (
            self._get_res_num(self.strands[s5_idx][0]),
            self._get_res_num(self.strands[s5_idx][-1])
        )
        processed.add(s5_idx)
        print(f"  → S5: тяж {s5_idx} ({path_map['S5'][0]}-{path_map['S5'][1]})")

        # Идем ВЛЕВО от S5: S6, S7, S8...
        current_idx = s5_idx
        current_number = 6
        while True:
            next_idx = None
            for neighbor in self.adj[current_idx]:
                if neighbor not in processed:
                    next_idx = neighbor
                    break
            if next_idx is None:
                break
            if next_idx in full_path:
                name = f"S{current_number}"
                strand_names[next_idx] = name
                path_map[name] = (
                    self._get_res_num(self.strands[next_idx][0]),
                    self._get_res_num(self.strands[next_idx][-1])
                )
                processed.add(next_idx)
                print(f"    → {name}: тяж {next_idx} ({path_map[name][0]}-{path_map[name][1]})")
                current_idx = next_idx
                current_number += 1
            else:
                break

        # ВТОРОЙ = S3 (если есть)
        if len(left_neighbors) > 1:
            s3_idx = left_neighbors[1][0]
            strand_names[s3_idx] = "S3"
            path_map["S3"] = (
                self._get_res_num(self.strands[s3_idx][0]),
                self._get_res_num(self.strands[s3_idx][-1])
            )
            processed.add(s3_idx)
            print(f"  → S3: тяж {s3_idx} ({path_map['S3'][0]}-{path_map['S3'][1]})")

            # Идем ВПРАВО от S3: S2, S1, S0...
            current_idx = s3_idx
            current_number = 2
            while True:
                next_idx = None
                for neighbor in self.adj[current_idx]:
                    if neighbor not in processed:
                        next_idx = neighbor
                        break
                if next_idx is None:
                    break
                if next_idx in full_path:
                    name = f"S{current_number}"
                    strand_names[next_idx] = name
                    path_map[name] = (
                        self._get_res_num(self.strands[next_idx][0]),
                        self._get_res_num(self.strands[next_idx][-1])
                    )
                    processed.add(next_idx)
                    print(f"    → {name}: тяж {next_idx} ({path_map[name][0]}-{path_map[name][1]})")
                    current_idx = next_idx
                    current_number -= 1
                else:
                    break

    # ------------------------------------------------------------
    # ОСТАЛЬНЫЕ ТЯЖИ В ПУТИ
    # ------------------------------------------------------------
    for i, idx in enumerate(full_path):
        if idx not in strand_names and idx != s4_idx:
            s_name = f"S{4 - (i - s4_pos)}"
            strand_names[idx] = s_name
            path_map[s_name] = (
                self._get_res_num(self.strands[idx][0]),
                self._get_res_num(self.strands[idx][-1])
            )
            processed.add(idx)

    # =================================================================
    # ПРОВЕРКА S3 И СОЗДАНИЕ СИСТЕМЫ КООРДИНАТ
    # =================================================================
    if s3_idx is None:
        print(f"❌ ОШИБКА: S3 не найден!")
        self.strands = original_strands
        self.adj = original_adj
        self.helices = original_helices
        return None

    # Создаем систему координат
    self._setup_coordinate_system(s4_idx, s3_idx)

    # =================================================================
    # ОПРЕДЕЛЕНИЕ СТОРОН СПИРАЛЕЙ
    # =================================================================
    print(f"\n🔍 ОПРЕДЕЛЕНИЕ СТОРОН СПИРАЛЕЙ ПО СИСТЕМЕ КООРДИНАТ:")

    for h_keys in self.helices:
        if len(h_keys) < self.MIN_HELIX_LENGTH:
            continue

        h_c = np.array([self.res_data[r]['coords'] for r in h_keys])
        h_center = h_c.mean(axis=0)
        h_start = self._get_res_num(h_keys[0])
        h_end = self._get_res_num(h_keys[-1])

        min_dist = float('inf')
        nearest_idx = None
        nearest_center = None
        for idx in full_path:
            s_coords = np.array([self.res_data[key]['coords'] for key in self.strands[idx]])
            s_center = s_coords.mean(axis=0)
            dist = np.linalg.norm(h_center - s_center)
            if dist < min_dist:
                min_dist = dist
                nearest_idx = idx
                nearest_center = s_center

        if min_dist < self.HELIX_RADIUS:
            nearest_name = strand_names.get(nearest_idx, f"S?({nearest_idx})")
            
            side, proj = self._get_helix_side_by_coords(h_center, nearest_center)

            self.helix_sides[h_start] = side
            self.helix_distances[h_start] = min_dist
            self.helix_nearest_strand[h_start] = nearest_name

            print(f"  Спираль {h_start:3d}-{h_end:3d}: проекция {proj:7.3f} -> {side:2s} | ближ.тяж {nearest_name} ({min_dist:.2f} Å)")
        else:
            print(f"  Спираль {h_start:3d}-{h_end:3d}: пропущена (расстояние {min_dist:.2f} Å > {self.HELIX_RADIUS} Å)")

    # =================================================================
    # ВЫВОД ТАБЛИЦЫ
    # =================================================================
    all_strands = []
    for idx, name in strand_names.items():
        strand = self.strands[idx]
        s_start = self._get_res_num(strand[0])
        s_end = self._get_res_num(strand[-1])
        all_strands.append({'idx': idx, 'name': name, 'start': s_start, 'end': s_end, 'strand': strand})

    def sort_key_table(x):
        name = x['name']
        if name.startswith('S'):
            num_str = name[1:].split('_')[0].split('(')[0]
            try:
                num = int(num_str)
                return -num
            except:
                return 0
        return 0
    all_strands.sort(key=sort_key_table)

    print(f"\n✅ MOTIF: {motif_text} ({motif_res}) | S4: {s4_start}-{s4_end}")
    print("-" * 120)
    print(f"{'Strand':<7} | {'Range':<10} | {'Dir':<7} | {'Bond':<10} | {'Hu (Up)':<30} | {'Hd (Down)'}")
    print("-" * 120)

    for item in all_strands:
        idx = item['idx']
        s_name = item['name']
        s_start = item['start']
        s_end = item['end']
        strand_keys = item['strand']

        vi = self._get_strand_vector(strand_keys)
        dir_str = "UP" if np.dot(vi, v4) > 0 else "DOWN"

        if idx == full_path[0]:
            bond = "Edge"
        else:
            try:
                pos = full_path.index(idx)
                if pos > 0:
                    prev_idx = full_path[pos-1]
                    prev_vi = self._get_strand_vector(self.strands[prev_idx])
                    bond = "Para" if np.dot(vi, prev_vi) > 0 else "Anti"
                else:
                    bond = "Edge"
            except:
                bond = "Edge"

        hu_list, hd_list = [], []
        s_coords = np.array([self.res_data[key]['coords'] for key in strand_keys])

        for h_keys in self.helices:
            if len(h_keys) < self.MIN_HELIX_LENGTH:
                continue

            h_c = np.array([self.res_data[r]['coords'] for r in h_keys])
            if np.min(distance_matrix(s_coords, h_c)) < self.HELIX_RADIUS:
                h_start = self._get_res_num(h_keys[0])
                h_end = self._get_res_num(h_keys[-1])

                if h_start in self.helix_sides:
                    side = self.helix_sides[h_start]
                    dist = self.helix_distances[h_start]
                    num_part = self._get_helix_number(h_start, h_end, path_map)
                    name = self._get_helix_name(side, h_start, num_part)
                    label = f"{name} ({h_start}-{h_end}) [{dist:.1f} Å]"

                    if side == "Hu":
                        if label not in hu_list:
                            hu_list.append(label)
                    else:
                        if label not in hd_list:
                            hd_list.append(label)

        hu_str = ', '.join(hu_list) if hu_list else '--'
        hd_str = ', '.join(hd_list) if hd_list else '--'
        print(f"{s_name:<7} | {s_start:>4}-{s_end:<5} | {dir_str:<7} | {bond:<10} | {hu_str:<30} | {hd_str}")

    # =================================================================
    # СОХРАНЯЕМ РЕЗУЛЬТАТ И ВОССТАНАВЛИВАЕМ ОРИГИНАЛЫ
    # =================================================================
    result = {
        'full_path': full_path,
        'path_map': path_map,
        'strands': self.strands,
        'helices': self.helices,
        'strand_names': strand_names,
        's4_idx': s4_idx,
        's4_start': s4_start,
        's4_end': s4_end,
        'v4': v4,
        'chain': motif_chain,
        's4_global_idx': s4_global_idx,
        's3_idx': s3_idx,
        'helix_sides': self.helix_sides,
        'helix_distances': self.helix_distances,
        'helix_nearest_strand': self.helix_nearest_strand,
        'coord_system': self.coord_system
    }

    # Восстанавливаем оригиналы
    self.strands = original_strands
    self.adj = original_adj
    self.helices = original_helices
    self.current_chain = None
    self.coord_system = None
    self.strand_names = None

    return result
# =====================================================================
# ЛИНЕЙНАЯ ТОПОЛОГИЯ
# =====================================================================
def print_linear_topology_from_result(self, result):
    if not result:
        return
    full_path = result['full_path']
    path_map = result['path_map']
    strand_names = result['strand_names']
    s4_idx = result['s4_idx']
    v4 = result['v4']
    s4_pos = full_path.index(s4_idx)
    all_elements = []
    unique_helices = {}
    motif_chain = result['chain']

    for h_keys in self.helices:
        if self._get_chain(h_keys[0]) != motif_chain:
            continue
        
        if len(h_keys) >= self.MIN_HELIX_LENGTH:
            h_start = self._get_res_num(h_keys[0])
            h_end = self._get_res_num(h_keys[-1])

            if h_start in result['helix_sides']:
                helix_key = f"{h_start}-{h_end}"
                if helix_key not in unique_helices:
                    side = result['helix_sides'][h_start]
                    dist = result['helix_distances'][h_start]
                    num_part = self._get_helix_number(h_start, h_end, path_map)
                    name = self._get_helix_name(side, h_start, num_part)
                    unique_helices[helix_key] = (h_start, f"{name}[{h_start}-{h_end}] ({dist:.1f} Å)")

    for h_start, helix_name in sorted(unique_helices.values(), key=lambda x: x[0]):
        all_elements.append((h_start, helix_name))

    for i, idx in enumerate(full_path):
        if idx in strand_names:
            s_name = strand_names[idx]
        else:
            s_name = f"S{4 - (i - s4_pos)}"
        strand_keys = result['strands'][idx]
        s_start = self._get_res_num(strand_keys[0])
        s_end = self._get_res_num(strand_keys[-1])
        vi = self._get_strand_vector(strand_keys)
        direction = "↑" if np.dot(vi, v4) > 0 else "↓"
        strand_name = f"{s_name}({direction})[{s_start}-{s_end}]"
        all_elements.append((s_start, strand_name))

    all_elements.sort(key=lambda x: x[0])
    print(f"\n{'='*80}")
    print("ЛИНЕЙНАЯ ТОПОЛОГИЯ (N -> C):")
    print(" — ".join([name for _, name in all_elements]))
    print('='*80)

MTaseAnalyzer.analyze_topology = analyze_topology
MTaseAnalyzer.print_linear_topology_from_result = print_linear_topology_from_result