import numpy as np

from .core import MTaseAnalyzer

def _get_helix_number(self, h_start, h_end, path_map):
    def get_r(name):
        return path_map.get(name)
    if get_r('S7') and get_r('S4') and get_r('S7')[1] < h_start < get_r('S4')[0]:
        return "3"
    if get_r('S4') and h_start > get_r('S4')[1]:
        if not get_r('S3') or h_end < get_r('S3')[0]:
            return "2"
    if get_r('S3') and get_r('S1') and get_r('S3')[1] < h_start < get_r('S1')[0]:
        return "1"
    if get_r('S4') and get_r('S5') and get_r('S4')[1] < h_start < get_r('S5')[0]:
        return "1"
    if get_r('S5') and get_r('S6') and get_r('S5')[1] < h_start < get_r('S6')[0]:
        return "2"
    if get_r('S6') and get_r('S7') and get_r('S6')[1] < h_start < get_r('S7')[0]:
        return "3"
    return ""

def _determine_helix_side(self, proj):
    return "Hu" if proj > 0 else "Hd"

def _get_helix_name(self, side, h_start, num_part):
    if num_part:
        return f"{side}{num_part}"
    else:
        return f"{side}_{h_start}"

def _find_hbond_between_strands(self, strand1, strand2):
    """Поиск водородной связи между двумя тяжами по координатам Cα"""
    min_dist = float('inf')
    best_pair = None
    best_vector = None
    
    for key1 in strand1:
        coord1 = self.res_data[key1]['coords']
        res1 = self._get_res_num(key1)
        
        for key2 in strand2:
            coord2 = self.res_data[key2]['coords']
            res2 = self._get_res_num(key2)
            
            dist = np.linalg.norm(coord1 - coord2)
            if dist < 7.0 and dist < min_dist:
                min_dist = dist
                best_pair = (res1, res2, dist)
                best_vector = coord2 - coord1
                
    return best_vector, best_pair

def _setup_coordinate_system(self, s4_idx, s3_idx):
    """
    СОЗДАНИЕ СИСТЕМЫ КООРДИНАТ ДЛЯ ТЕКУЩЕЙ ЦЕПИ
    """
    print(f"\n{'='*60}")
    print(f"🌍 СОЗДАНИЕ СИСТЕМЫ КООРДИНАТ ДЛЯ ЦЕПИ {self.current_chain}")
    print(f"{'='*60}")
    
    if s4_idx >= len(self.strands) or s3_idx >= len(self.strands):
        raise ValueError(f"Индекс тяжа вне диапазона")
    
    s4_strand = self.strands[s4_idx]
    s3_strand = self.strands[s3_idx]
    
    s4_center = self._get_strand_center(s4_strand)
    s3_center = self._get_strand_center(s3_strand)
    
    print(f"\n📌 ЦЕНТР МИРА (S4 цепи {self.current_chain}):")
    print(f"   [{s4_center[0]:>8.2f}, {s4_center[1]:>8.2f}, {s4_center[2]:>8.2f}]")
    
    s4_first = s4_strand[0]
    s4_last = s4_strand[-1]
    
    s4_first_ca = self.res_data[s4_first]['coords']
    s4_last_ca = self.res_data[s4_last]['coords']
    
    north = s4_last_ca - s4_first_ca
    north_norm = north / np.linalg.norm(north)
    
    print(f"\n🧭 СЕВЕР (+Y): N→C вдоль S4")
    print(f"   Вектор: [{north_norm[0]:>8.3f}, {north_norm[1]:>8.3f}, {north_norm[2]:>8.3f}]")
    
    hbond_vector, hbond_pair = self._find_hbond_between_strands(s4_strand, s3_strand)
    
    if hbond_vector is not None:
        east = hbond_vector / np.linalg.norm(hbond_vector)
        print(f"\n🧭 ВОСТОК (+X): ОТ S4 К S3 (по водородной связи)")
        print(f"   Вектор: [{east[0]:>8.3f}, {east[1]:>8.3f}, {east[2]:>8.3f}]")
    else:
        east = s3_center - s4_center
        east = east / np.linalg.norm(east)
        print(f"\n🧭 ВОСТОК (+X): ОТ S4 К S3 (геометрический центр)")
    
    up = np.cross(north_norm, east)
    up_norm = up / np.linalg.norm(up)
    up_norm = -up_norm
    
    print(f"\n⬆️ ВВЕРХ (+Z): СЕВЕР × ВОСТОК")
    print(f"   Вектор: [{up_norm[0]:>8.3f}, {up_norm[1]:>8.3f}, {up_norm[2]:>8.3f}]")
    
    self.coord_system = {
        'chain': self.current_chain,
        's4_center': s4_center,
        's4_idx': s4_idx,
        's3_idx': s3_idx,
        'north': north_norm,
        'east': east,
        'up': up_norm,
        'hbond_used': hbond_pair
    }
    
    print(f"\n✅ СИСТЕМА КООРДИНАТ СОЗДАНА ДЛЯ ЦЕПИ {self.current_chain}")
    print(f"{'='*60}")
    
    return self.coord_system

def _get_helix_side_by_coords(self, helix_center, nearest_strand_center):
    """Определение стороны спирали ОТНОСИТЕЛЬНО БЛИЖАЙШЕГО ТЯЖА"""
    if self.coord_system is None:
        raise ValueError(f"Система координат не создана для цепи {self.current_chain}")
    
    # ✅ Вектор ОТ БЛИЖАЙШЕГО ТЯЖА К СПИРАЛИ!
    v = helix_center - nearest_strand_center
    proj = np.dot(v, self.coord_system['up'])
    return "Hu" if proj > 0 else "Hd", proj

# Прикрепляем методы
MTaseAnalyzer._find_hbond_between_strands = _find_hbond_between_strands
MTaseAnalyzer._setup_coordinate_system = _setup_coordinate_system
MTaseAnalyzer._get_helix_side_by_coords = _get_helix_side_by_coords
MTaseAnalyzer._get_helix_number = _get_helix_number
MTaseAnalyzer._determine_helix_side = _determine_helix_side
MTaseAnalyzer._get_helix_name = _get_helix_name