import streamlit as st

def render_sidebar():
    with st.sidebar:
        st.header("📥 **Input**")
        
        # Предустановленные MTases
        st.subheader("Pre-studied MTases")
        preset_mtases = {
            "None": None,
            "M.TaqI (PDB: 1G38)": "1G38",
            "M.HhaI (PDB: 5MHT)": "5MHT",
            "M.EcoRI (PDB: 1QPS)": "1QPS",
            "M.EcoRV (PDB: 4XQK)": "4XQK",
            "M.Bse634I (PDB: 1KRA)": "1KRA",
            "M.Thermophilus (PDB: 7M6B)": "7M6B",
            "Custom structure": "custom"
        }
        
        selected_preset = st.selectbox(
            "Select MTase:",
            list(preset_mtases.keys()),
            index=0,
            key="preset"
        )
        
        st.markdown("---")
        
        # Выбор источника структуры
        st.subheader("Structure Input")
        input_type = st.radio(
            "Choose input method:",
            ["PDB ID", "AlphaFold ID", "Upload PDB file"],
            key="input_type"
        )
        
        if input_type == "PDB ID":
            pdb_id = st.text_input(
                "Enter PDB ID:",
                value=preset_mtases[selected_preset] if selected_preset != "None" and preset_mtases[selected_preset] != "custom" else "",
                placeholder="4XQK",
                key="pdb_id"
            ).upper()
        elif input_type == "AlphaFold ID":
            uniprot_id = st.text_input(
                "Enter UniProt ID:",
                placeholder="P04392",
                key="uniprot_id"
            ).upper()
        else:
            uploaded_file = st.file_uploader(
                "Upload PDB file:",
                type=['pdb', 'ent', 'cif'],
                key="uploaded_file"
            )
        
        st.markdown("---")
        
        # Параметры анализа
        st.subheader("Analysis Parameters")
        
        motif_choice = st.selectbox(
            "Catalytic motifs to search:",
            ["Default (DPPY, NPPY, PC, PS)", "Custom motifs"],
            key="motif_choice"
        )
        
        if motif_choice == "Custom motifs":
            custom_motifs = st.text_area(
                "Enter motifs:",
                placeholder="[SND]PP[YFW]\nP[CS]",
                key="custom_motifs"
            )
        
        st.markdown("---")
        
        # Кнопка запуска
        run_button = st.button(
            "🔬 **Run Analysis**", 
            type="primary", 
            use_container_width=True,
            key="run_button"
        )