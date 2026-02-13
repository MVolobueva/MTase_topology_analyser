import streamlit as st
import io
import sys
import tempfile
import os

def show_linear_topology(result, analyzer):
    """Показывает линейную топологию"""
    old_stdout = sys.stdout
    new_stdout = io.StringIO()
    sys.stdout = new_stdout
    
    analyzer.print_linear_topology_from_result(result)
    
    sys.stdout = old_stdout
    linear_topology = new_stdout.getvalue()
    
    st.text(linear_topology)
    
    st.download_button(
        "📥 Download",
        linear_topology,
        file_name=f"linear_topology.txt"
    )

def show_2d_topology(result, analyzer):
    """Показывает 2D топологию"""
    fig = analyzer.visualize_topology_from_analysis(result)
    if fig:
        st.pyplot(fig)
    else:
        st.warning("Could not generate 2D topology")

def show_3d_topology(result, analyzer, pdb_id):
    """Показывает 3D структуру"""
    view = analyzer.visualize_3d_structure(
        result, 
        pdb_id=pdb_id,
        chain=result['chain']
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