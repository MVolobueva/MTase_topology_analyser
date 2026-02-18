import streamlit as st
import plotly.graph_objects as go
import tempfile
import os
import re
import io
import sys

def show_linear_topology(result, analyzer):
    """Показывает линейную топологию с BIG LOOP для больших разрывов"""
    if not result:
        st.warning("No topology data available")
        return
    
    # Получаем сырую линейную топологию из анализатора
    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout
    
    analyzer.print_linear_topology_from_result(result)
    
    sys.stdout = old_stdout
    linear_topology = new_stdout.getvalue()
    
    # Извлекаем только строку с элементами
    lines = linear_topology.strip().split('\n')
    topology_line = ""
    
    for line in lines:
        # Ищем строку, которая содержит элементы (со скобками и стрелками)
        if ('↑' in line or '↓' in line) and '[' in line and ']' in line:
            topology_line = line.strip()
            break
    
    if not topology_line:
        st.warning("Could not parse topology")
        return
    
    # Убираем все расстояния вида (15.8 Å)
    topology_line = re.sub(r'\s*\([\d.]+ Å\)', '', topology_line)
    
    # Разбираем на элементы для добавления BIG LOOP
    elements = topology_line.split(' — ')
    
    # Проверяем разрывы и добавляем BIG LOOP
    processed_elements = []
    for i, elem in enumerate(elements):
        processed_elements.append(elem)
        
        # Извлекаем конечный номер из текущего элемента
        current_end_match = re.search(r'\[(\d+)-(\d+)\]', elem)
        if i < len(elements) - 1 and current_end_match:
            current_end = int(current_end_match.group(2))
            
            # Извлекаем начальный номер из следующего элемента
            next_elem = elements[i + 1]
            next_start_match = re.search(r'\[(\d+)-(\d+)\]', next_elem)
            if next_start_match:
                next_start = int(next_start_match.group(1))
                
                # Если разрыв больше 50 а.к.
                if next_start - current_end > 50:
                    processed_elements.append("🔴 **BIG LOOP** 🔴")
    
    # Собираем обратно
    final_topology = ' — '.join(processed_elements)
    
    # Отображаем
    st.info(f"**{final_topology}**")

def show_2d_topology(result, analyzer):
    """Показывает 2D топологию (статическая версия)"""
    fig = analyzer.visualize_topology_from_analysis(result)
    if fig:
        st.pyplot(fig)
    else:
        st.warning("Could not generate 2D topology")


def show_3d_topology(result, analyzer, pdb_id=None, pdb_file=None):
    """Показывает 3D структуру"""
    view = analyzer.visualize_3d_structure(
        result, 
        pdb_id=pdb_id,
        chain=result['chain'],
        pdb_file=pdb_file
    )
    
    if view:
        with tempfile.NamedTemporaryFile(suffix='.html', delete=False) as f:
            view.write_html(f.name)
            with open(f.name, 'r') as f_html:
                html_content = f_html.read()
            st.components.v1.html(html_content, height=600)
        os.unlink(f.name)
    else:
        st.warning("Could not generate 3D visualization")