import numpy as np
import os
import collections
import re
from scipy.spatial import distance_matrix

class MTaseAnalyzer:
    def __init__(self, contact_dist=5.2, helix_radius=20, max_loop=5, min_helix_length=4):
        self.CONTACT_DIST = contact_dist
        self.HELIX_RADIUS = helix_radius
        self.MAX_LOOP = max_loop
        self.MOTIF_PATTERNS = [r"[SND]P[PL][YFW]", r"P[CSP]"]
        self.MIN_HELIX_LENGTH = min_helix_length

        self.res_data = {}
        self.full_seq = ""
        self.res_map = []
        self.strands = []
        self.helices = []
        self.adj = collections.defaultdict(set)
        self.motif_info = None
        self.chain_data = {}
        self.helix_sides = {}
        self.helix_distances = {}
        self.helix_nearest_strand = {}
        self.current_chain = None
        self.coord_system = None
        self.strand_names = None

        self.COLORS = {
            'Hu': '#27ae60',
            'Hd': '#e74c3c',
            'strand': '#1a5276'
        }

    def _get_res_num(self, key):
        if ':' in key:
            return int(key.split(':')[1])
        try:
            return int(key)
        except:
            return 0

    def _get_chain(self, key):
        if ':' in key:
            return key.split(':')[0]
        return 'A'

    def load_dssp(self, file_path):
        if not os.path.exists(file_path):
            print(f"Ошибка: Файл {file_path} не найден")
            return False

        self.res_data, self.full_seq, self.res_map = {}, "", []
        self.chain_data = {}

        with open(file_path, 'r') as f:
            lines = f.readlines()

        header_idx = next(i for i, l in enumerate(lines) if "  #  RESIDUE" in l) + 1

        for line in lines[header_idx:]:
            try:
                rn = int(line[5:10].strip())
                chain_id = line[11].strip()
                if not chain_id:
                    chain_id = 'A'

                unique_id = f"{chain_id}:{rn}"

                self.res_data[unique_id] = {
                    'struct': line[16],
                    'coords': np.array([float(line[115:122]), float(line[123:130]), float(line[131:138])]),
                    'aa': line[13] if len(line) > 13 else 'X',
                    'chain': chain_id,
                    'res_num': rn,
                    'unique_id': unique_id
                }
                self.full_seq += self.res_data[unique_id]['aa']
                self.res_map.append(unique_id)

                if chain_id not in self.chain_data:
                    self.chain_data[chain_id] = []
                self.chain_data[chain_id].append(unique_id)

            except:
                continue

        print(f"\n{'='*60}")
        print("ИНФОРМАЦИЯ О ЦЕПЯХ В DSSP:")
        for chain_id, residues in self.chain_data.items():
            res_nums = [self._get_res_num(r) for r in residues]
            print(f"Цепь '{chain_id}': {len(residues)} остатков ({min(res_nums)}-{max(res_nums)})")

        return True

    def find_all_motifs(self):
        """Найти все мотивы во всех цепях"""
        motifs = []

        for pattern in self.MOTIF_PATTERNS:
            for m in re.finditer(pattern, self.full_seq):
                seq_pos = m.start()
                if seq_pos < len(self.res_map):
                    motif_key = self.res_map[seq_pos]
                    if motif_key in self.res_data:
                        chain = self.res_data[motif_key]['chain']
                        motif_res_num = self.res_data[motif_key]['res_num']

                        potential = []
                        for i, strand in enumerate(self.strands):
                            if strand:
                                last_key = strand[-1]
                                if last_key in self.res_data:
                                    last_chain = self.res_data[last_key]['chain']
                                    last_num = self.res_data[last_key]['res_num']
                                    if last_chain == chain and last_num < motif_res_num:
                                        potential.append((i, last_num))

                        if potential:
                            idx, last_num = max(potential, key=lambda x: x[1])
                            if motif_res_num - last_num <= self.MAX_LOOP + 1:
                                motifs.append({
                                    'text': m.group(),
                                    'res': motif_res_num,
                                    'key': motif_key,
                                    'chain': chain,
                                    's4_idx': idx,
                                    's4_start': self._get_res_num(self.strands[idx][0]),
                                    's4_end': last_num
                                })
                                print(f"  ✅ Цепь {chain}: мотив {m.group()} ({motif_res_num}), S4 {last_num} (глобальный индекс {idx})")

        print(f"\n📊 Всего найдено мотивов: {len(motifs)}")
        return motifs

    def _merge_helices(self, helices):
        """Объединение перекрывающихся спиралей по 3D расстоянию"""
        if not helices:
            return helices

        chains = {}
        for h in helices:
            chain = self._get_chain(h[0])
            if chain not in chains:
                chains[chain] = []
            chains[chain].append(h)

        merged = []
        for chain, chain_helices in chains.items():
            chain_helices.sort(key=lambda x: self._get_res_num(x[0]))

            current = chain_helices[0].copy()
            for next_helix in chain_helices[1:]:
                curr_end_key = current[-1]
                next_start_key = next_helix[0]

                curr_end_coord = self.res_data[curr_end_key]['coords']
                next_start_coord = self.res_data[next_start_key]['coords']

                dist_3d = np.linalg.norm(curr_end_coord - next_start_coord)

                if dist_3d <= 5.0:
                    combined = list(set(current + next_helix))
                    combined.sort(key=lambda x: self._get_res_num(x))
                    current = combined
                else:
                    merged.append(current)
                    current = next_helix.copy()
            merged.append(current)

        unique = []
        seen = set()
        for h in merged:
            h_start = self._get_res_num(h[0])
            h_end = self._get_res_num(h[-1])
            h_chain = self._get_chain(h[0])
            key = f"{h_chain}:{h_start}-{h_end}"
            if key not in seen:
                seen.add(key)
                unique.append(h)

        return unique

    def find_all_strands(self):
        """Нахождение всех бета-тяжей и альфа-спиралей"""
        self.strands, self.helices = [], []

        for chain_id in self.chain_data:
            chain_keys = sorted(self.chain_data[chain_id], key=lambda x: self._get_res_num(x))
            curr_s, curr_h = [], []

            for key in chain_keys:
                st = self.res_data[key]['struct']

                if st == 'E':
                    if not curr_s:
                        curr_s.append(key)
                    else:
                        last_num = self._get_res_num(curr_s[-1])
                        curr_num = self._get_res_num(key)
                        if curr_num == last_num + 1:
                            curr_s.append(key)
                        else:
                            if curr_s:
                                self.strands.append(curr_s)
                                curr_s = []
                            curr_s = [key]
                elif st in ['H', 'G', 'I']:
                    if not curr_h:
                        curr_h.append(key)
                    else:
                        last_num = self._get_res_num(curr_h[-1])
                        curr_num = self._get_res_num(key)
                        if curr_num - last_num <= 5:
                            curr_h.append(key)
                        else:
                            if curr_h:
                                self.helices.append(curr_h)
                            curr_h = [key]

            if curr_s:
                self.strands.append(curr_s)
            if curr_h:
                self.helices.append(curr_h)

        self.helices = self._merge_helices(self.helices)
        print(f"Найдено спиралей после объединения: {len(self.helices)}")
        return self.strands, self.helices

    def _get_strand_vector(self, strand_keys):
        start_coord = self.res_data[strand_keys[0]]['coords']
        end_coord = self.res_data[strand_keys[-1]]['coords']
        return end_coord - start_coord

    def _get_strand_center(self, strand_keys):
        coords = []
        for key in strand_keys:
            if key in self.res_data:
                coords.append(self.res_data[key]['coords'])
        return np.array(coords).mean(axis=0)

    def build_sheet_adjacency(self):
        self.adj = collections.defaultdict(set)
        for i in range(len(self.strands)):
            for j in range(i + 1, len(self.strands)):
                s1_c = np.array([self.res_data[r]['coords'] for r in self.strands[i]])
                s2_c = np.array([self.res_data[r]['coords'] for r in self.strands[j]])
                if np.min(distance_matrix(s1_c, s2_c)) < self.CONTACT_DIST:
                    self.adj[i].add(j)
                    self.adj[j].add(i)
        return self.adj

    def _expand_path(self, start_node, visited):
        p_list = [start_node]
        while True:
            ns = [x for x in self.adj[p_list[-1]] if x not in visited]
            if not ns:
                break
            nxt = min(ns, key=lambda x: abs(
                self._get_res_num(self.strands[x][0]) - self._get_res_num(self.strands[p_list[-1]][0])))
            p_list.append(nxt)
            visited.add(nxt)
        return p_list

