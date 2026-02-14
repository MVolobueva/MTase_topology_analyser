import streamlit as st
import pandas as pd

def show_results_table(results):
    """Отображает таблицу с найденными мотивами"""
    chains_info = []
    for key, data in results.items():
        motif = data['motif']
        result = data['result']
        chains_info.append({
            'Chain': data['display_chain'],
            'Motif': data['display_motif'],
            'Strands': len(result['full_path']),
            'S4': f"{result['s4_start']}-{result['s4_end']}",
            'key': key
        })
    
    df = pd.DataFrame(chains_info)
    st.dataframe(df[['Chain', 'Motif', 'Strands', 'S4']], 
                width='stretch', 
                hide_index=True)
    
    return chains_info
    