import streamlit as st

def show():
    st.header("📚 **Documentation**")
    
    with st.expander("ℹ️ About", expanded=True):
        st.markdown("""
        ### MTase Topology Analyzer
        
        This tool analyzes the topology of MTase catalytic domains.
        
        **Interpretation:**
        - **S4** - catalytic strand containing the motif
        - **S3, S2, S1** - strands on C-terminal side
        - **S5, S6, S7** - strands on N-terminal side
        - **Hu** - helices above the sheet (green)
        - **Hd** - helices below the sheet (red)
        """)