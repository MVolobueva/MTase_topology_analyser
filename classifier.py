#!/usr/bin/env python3
"""
Module for determining structural class of MTase based on beta-sheet topology
Classes: A, B, C, D, E, F
Pure topology-based - NO MOTIF REQUIRED!
"""

import re

# Minimum number of strands required for canonical MTase catalytic domain
MIN_STRAND_COUNT = 6


def is_s0_at_n_terminus(strand_list):
    """Check if S0 appears before S5 (N-terminus)"""
    if 'S0' not in strand_list or 'S5' not in strand_list:
        return False
    s0_idx = strand_list.index('S0')
    s5_idx = strand_list.index('S5')
    return s0_idx < s5_idx


def is_s0_at_c_terminus(strand_list):
    """Check if S0 appears after S2 (C-terminus)"""
    if 'S0' not in strand_list or 'S2' not in strand_list:
        return False
    s0_idx = strand_list.index('S0')
    s2_idx = strand_list.index('S2')
    return s0_idx > s2_idx


def is_permuted_order(strand_list):
    """Checks if strand order is permuted (S7 appears before S6)"""
    s6_index = strand_list.index('S6') if 'S6' in strand_list else -1
    s7_index = strand_list.index('S7') if 'S7' in strand_list else -1
    return (s7_index < s6_index) if (s6_index != -1 and s7_index != -1) else False


def classify_by_topology(strand_list, has_s0, has_s_minus_1, gap_s6_s7, has_n_helix, has_c_helix):
    """
    Pure topology-based classification - NO MOTIF REQUIRED!
    
    Parameters:
    - strand_list: order of strands in N→C (e.g., ['S5','S6','S7','S4','S3','S1','S2'])
    - has_s0: presence of S0 strand
    - has_s_minus_1: presence of S-1 strand
    - gap_s6_s7: distance between S6 and S7 (in aa), None if S6 or S7 missing
    - has_n_helix: presence of N-terminal helix (before first strand)
    - has_c_helix: presence of C-terminal helix (after last strand)
    
    Returns:
    - class_name: A, B, C, D, E, F, or Unknown
    - confidence: high/medium/low
    - reasons: list of reasons for classification
    """
    reasons = []
    
    # ========== MINIMUM STRAND CHECK ==========
    # Canonical MTase requires at least 6 strands in the sheet
    if len(strand_list) < MIN_STRAND_COUNT:
        reasons.append(f"Too few strands ({len(strand_list)} < {MIN_STRAND_COUNT}) - incomplete catalytic domain")
        return 'Unknown', 'low', reasons
    
    # ========== CLASS A (S0 + S-1 both present at C-terminus) ==========
    if has_s0 and has_s_minus_1 and is_s0_at_c_terminus(strand_list):
        reasons.append("S0 and S-1 present at C-terminus (C-terminal extension)")
        return 'A', 'high', reasons
    
    # ========== CLASS B (permuted order) ==========
    # Check if S7 appears early (within first 3 positions) and S5/S6 appear at the end
    # Handle cases where S0 may appear before S7
    s7_index = -1
    for i, s in enumerate(strand_list):
        if s == 'S7':
            s7_index = i
            break
    
    s5_s6_at_end = ('S5' in strand_list[-2:] or 'S6' in strand_list[-2:])
    
    if s7_index >= 0 and s7_index <= 3 and s5_s6_at_end:
        reasons.append(f"Permuted order: S7 at position {s7_index}, S5/S6 at C-terminus")
        return 'B', 'high', reasons
    
    # ========== HYBRID D/E (S0 at N-terminus + large gap S6-S7) ==========
    # Special case: shares features with both Class D (large gap) and Class E (S0 at N-terminus)
    if (has_s0 and not has_s_minus_1 and 
        is_s0_at_n_terminus(strand_list) and
        gap_s6_s7 is not None and gap_s6_s7 >= 100):
        reasons.append(f"S0 at N-terminus with large gap S6-S7 ({gap_s6_s7} aa) - hybrid D/E feature")
        return 'D', 'medium', reasons
    
    # ========== CLASS E (S0 at N-terminus, no S-1, small gap S6-S7) ==========
    if (has_s0 and not has_s_minus_1 and 
        is_s0_at_n_terminus(strand_list) and
        gap_s6_s7 is not None and gap_s6_s7 < 100 and
        has_n_helix):
        reasons.append(f"S0 at N-terminus, small gap S6-S7 ({gap_s6_s7} aa), no S-1")
        return 'E', 'high', reasons
    
    # ========== CLASS D (large gap S6-S7 = TRD insertion) ==========
    if (not has_s0 and not has_s_minus_1 and 
        gap_s6_s7 is not None and gap_s6_s7 >= 100 and
        has_n_helix and not has_c_helix):
        reasons.append(f"Large gap S6-S7 ({gap_s6_s7} aa) - TRD insertion between Hd3 and S7")
        return 'D', 'high', reasons
    
    # ========== CLASS C (C-helix only, no N-helix, no S0/S-1) ==========
    if (not has_s0 and not has_s_minus_1 and 
        has_c_helix and not has_n_helix):
        if gap_s6_s7 is not None:
            reasons.append(f"Small gap S6-S7 ({gap_s6_s7} aa), C-helix present")
        else:
            reasons.append("Reduced beta-sheet, C-helix present")
        return 'C', 'high', reasons
    
    # ========== CLASS F (N-helix present, no S0/S-1) ==========
    if (not has_s0 and not has_s_minus_1 and has_n_helix):
        reasons.append("N-terminal helix present")
        return 'F', 'high', reasons
    
    # ========== Partial Class A (S0 at C-terminus, S-1 missing) ==========
    # This gets MEDIUM confidence because S-1 is missing (incomplete C-terminus)
    if has_s0 and not has_s_minus_1 and is_s0_at_c_terminus(strand_list):
        reasons.append("S0 present at C-terminus but S-1 missing (incomplete C-terminal extension)")
        return 'A', 'medium', reasons
    
    # ========== UNKNOWN ==========
    reasons.append("No clear class assignment")
    return 'Unknown', 'low', reasons


def parse_strand_order(full_topology):
    """Extracts strand order from full topology string"""
    matches = re.findall(r'S(-?\d+)\([↑↓]\)', full_topology)
    return [f"S{m}" for m in matches]


def has_strand(strand_list, name):
    """Checks if a strand exists in the list"""
    return name in strand_list


def get_gap_between_strands(full_topology, s1, s2):
    """Calculates distance between two strands in amino acids"""
    pattern_s1 = rf'{s1}\([↑↓]\)\[(\d+)-(\d+)\]'
    pattern_s2 = rf'{s2}\([↑↓]\)\[(\d+)-(\d+)\]'
    
    match_s1 = re.search(pattern_s1, full_topology)
    match_s2 = re.search(pattern_s2, full_topology)
    
    if match_s1 and match_s2:
        end_s1 = int(match_s1.group(2))
        start_s2 = int(match_s2.group(1))
        return start_s2 - end_s1
    
    return None


def has_n_terminal_helix(full_topology):
    """Checks for presence of N-terminal helix before the first strand"""
    first_s = re.search(r'S-?\d+', full_topology)
    if first_s:
        before_first = full_topology[:first_s.start()]
        if re.search(r'Hd_\d+\[\d+-\d+\]', before_first):
            return True
    return False


def has_c_terminal_helix(full_topology):
    """Checks for presence of C-terminal helix after the last strand"""
    last_s = re.search(r'S-?\d+\([↑↓]\)\[\d+-\d+\](?!.*S-?\d+)', full_topology)
    if last_s:
        after_last = full_topology[last_s.end():]
        if re.search(r'H[ud]_\d+\[\d+-\d+\]', after_last):
            return True
    return False


def classify_topology(full_topology, motif=None):
    """
    Main classification function - uses pure topology, motif is ignored!
    
    Parameters:
    - full_topology: string with full secondary structure elements
    - motif: (optional) ignored, kept for backward compatibility
    
    Returns:
    - dictionary with classification results
    """
    strand_list = parse_strand_order(full_topology)
    
    has_s0 = has_strand(strand_list, 'S0')
    has_s_minus_1 = has_strand(strand_list, 'S-1')
    has_s6 = has_strand(strand_list, 'S6')
    has_s7 = has_strand(strand_list, 'S7')
    
    # Calculate gap only if in normal order (S6 before S7)
    gap_s6_s7 = None
    if has_s6 and has_s7 and not is_permuted_order(strand_list):
        gap_s6_s7 = get_gap_between_strands(full_topology, 'S6', 'S7')
    
    has_n_helix = has_n_terminal_helix(full_topology)
    has_c_helix = has_c_terminal_helix(full_topology)
    
    class_name, confidence, reasons = classify_by_topology(
        strand_list, has_s0, has_s_minus_1, gap_s6_s7, has_n_helix, has_c_helix
    )
    
    # Format gap for output
    gap_display = gap_s6_s7 if gap_s6_s7 is not None else 'N/A'
    if has_s6 and has_s7 and is_permuted_order(strand_list):
        gap_display = 'N/A (permuted)'
    
    return {
        'class': class_name,
        'confidence': confidence,
        'reasons': '; '.join(reasons),
        'strand_order': ' → '.join(strand_list),
        'has_s0': has_s0,
        'has_s-1': has_s_minus_1,
        'gap_s6_s7': gap_display,
        'has_n_helix': has_n_helix,
        'has_c_helix': has_c_helix
    }


# Test function
if __name__ == "__main__":
    test_cases = [
        # Class A example (3S1S)
        ("Hd_191[191-211] — Hd_301[301-315] — S5(↑)[324-327] — Hd2[334-341] — S6(↑)[351-354] — Hd3[358-369] — S7(↑)[383-385] — Hu3[389-397] — S4(↑)[400-405] — Hu2[443-454] — S3(↑)[460-466] — Hu1[468-484] — S1(↑)[490-495] — S2(↓)[509-515] — S0(↑)[522-528] — Hd_532[532-546] — S-1(↓)[558-564] — Hu_565[565-570]", None),
        
        # Class B example (6PBD)
        ("S7(↑)[6-10] — Hu3[13-19] — S4(↑)[25-30] — Hu2[56-80] — S3(↑)[81-91] — Hu1[96-105] — S1(↑)[109-118] — S2(↓)[134-141] — S0(↑)[168-171] — Hd1[195-205] — S5(↑)[211-214] — Hd2[221-228] — S6(↑)[232-237] — Hd_240[240-251]", None),
        
        # Class D example (2G1P)
        ("Hd_15[15-24] — S5(↑)[30-33] — Hd2[40-43] — S6(↑)[49-54] — Hd3[57-121] — Hd3[145-157] — S7(↑)[158-162] — Hu3[165-169] — S4(↑)[176-180] — Hu2[203-218] — S3(↑)[223-228] — Hu1[231-236] — S1(↑)[241-245] — S2(↓)[263-268]", None),
        
        # Class E example (Q8RNV6)
        ("Hd_31[31-35] — S0(↑)[41-43] — Hd_63[63-73] — S5(↑)[79-82] — Hd2[89-96] — S6(↑)[100-105] — Hd3[108-122] — S7(↑)[131-135] — Hu3[138-143] — S4(↑)[150-155] — Hu2[193-215] — S3(↑)[216-226] — Hu1[238-248] — S1(↑)[252-260] — S2(↓)[284-291]", None),
        
        # Hybrid D/E example (C1EIH6)
        ("Hd_53[53-58] — S0(↓)[64-68] — Hd_71[71-87] — Hd_112[112-121] — S5(↑)[126-129] — Hd2[136-143] — S6(↑)[147-152] — Hd3[155-165] — S7(↑)[353-356] — S4(↑)[376-381] — Hu2[391-402] — Hu2[440-445] — Hu2[466-494] — S3(↑)[495-506] — Hu1[513-523] — S1(↑)[527-535] — S2(↓)[551-558]", None),
        
        # Class C example (3G7U)
        ("S5(↑)[2-6] — S6(↑)[23-28] — S7(↑)[46-48] — S4(↑)[72-75] — S3(↑)[112-117] — S1(↑)[147-150] — S2(↓)[161-169]", None),
    ]
    
    for full_topology, motif in test_cases:
        result = classify_topology(full_topology, motif)
        print(f"\nClass: {result['class']} (confidence: {result['confidence']})")
        print(f"Strand order: {result['strand_order']}")
        print(f"Reasons: {result['reasons']}")