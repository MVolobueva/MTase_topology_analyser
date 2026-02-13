import streamlit as st
from analyzer import MTaseAnalyzer
from utils.helpers import download_structure, parse_uploaded_file
from components.results_table import show_results_table
from components.visualizations import show_linear_topology, show_2d_topology, show_3d_topology

def show():
    # Проверяем, что нужно запустить анализ
    run_button = st.session_state.get('run_button', False)
    
    # Определяем источник структуры
    structure_source = None
    if st.session_state.get('input_type') == "PDB ID" and st.session_state.get('pdb_id'):
        structure_source = st.session_state.pdb_id
        source_type = 'pdb'
    elif st.session_state.get('input_type') == "AlphaFold ID" and st.session_state.get('uniprot_id'):
        structure_source = st.session_state.uniprot_id
        source_type = 'alphafold'
    elif st.session_state.get('input_type') == "Upload PDB file" and st.session_state.get('uploaded_file'):
        structure_source = st.session_state.uploaded_file
        source_type = 'upload'
    
    if run_button and structure_source:
        with st.spinner("🔄 Loading structure and running analysis..."):
            try:
                # Создаем анализатор
                analyzer = MTaseAnalyzer()
                
                # Загружаем структуру
                if source_type == 'pdb':
                    dssp_file = download_structure(structure_source, source='pdb')
                elif source_type == 'alphafold':
                    dssp_file = download_structure(structure_source, source='alphafold')
                else:  # upload
                    dssp_file = parse_uploaded_file(structure_source)
                
                if not dssp_file:
                    st.error(f"Failed to download {source_type} structure")
                    return
                
                # Загружаем DSSP
                if not analyzer.load_dssp(dssp_file):
                    st.error("Failed to load DSSP file. Make sure DSSP is installed.")
                    return
                
                # Находим вторичные структуры
                analyzer.find_all_strands()
                analyzer.build_sheet_adjacency()
                
                # Находим мотивы
                if st.session_state.get('motif_choice') == "Default (DPPY, NPPY, PC, PS)":
                    motifs = analyzer.find_all_motifs()
                else:
                    original_patterns = analyzer.MOTIF_PATTERNS
                    custom = st.session_state.get('custom_motifs', '').strip().split('\n')
                    analyzer.MOTIF_PATTERNS = [m.strip() for m in custom if m.strip()]
                    motifs = analyzer.find_all_motifs()
                    analyzer.MOTIF_PATTERNS = original_patterns
                
                # Анализируем каждый мотив
                results = {}
                for motif in motifs:
                    chain = motif['chain']
                    motif_res = motif['res']
                    motif_text = motif['text']
                    key = f"{chain}_{motif_text}_{motif_res}"
                    
                    result = analyzer.analyze_topology(motif_data=motif)
                    if result:
                        results[key] = {
                            'motif': motif,
                            'result': result,
                            'display_chain': chain,
                            'display_motif': f"{motif_text} ({motif_res})"
                        }
                
                st.session_state.analyzer = analyzer
                st.session_state.results = results
                st.success(f"✅ Analysis complete! Found {len(motifs)} motifs in {len(results)} chains")
                
            except Exception as e:
                st.error(f"Error during analysis: {str(e)}")
                import traceback
                st.code(traceback.format_exc())
    
    # Отображаем результаты
    if st.session_state.get('results'):
        st.header("📊 **Analysis Results**")
        
        chains_info = show_results_table(st.session_state.results)
        
        if chains_info:
            selected_idx = st.selectbox(
                "Choose chain:",
                range(len(chains_info)),
                format_func=lambda i: chains_info[i]['Motif'],
                key="chain_selector"
            )
            
            selected_key = chains_info[selected_idx]['key']
            data = st.session_state.results[selected_key]
            
            st.markdown("---")
            st.header(f"🧬 **Chain {data['display_chain']}**")
            st.subheader(f"Motif: {data['display_motif']} | S4: {data['result']['s4_start']}-{data['result']['s4_end']}")
            
            tab1, tab2, tab3 = st.tabs(["📝 Linear Topology", "🖼️ 2D Topology", "🔬 3D Structure"])
            
            with tab1:
                show_linear_topology(data['result'], st.session_state.analyzer)
            with tab2:
                show_2d_topology(data['result'], st.session_state.analyzer)
            with tab3:
                pdb_for_view = st.session_state.pdb_id if st.session_state.get('input_type') == "PDB ID" else None
                show_3d_topology(data['result'], st.session_state.analyzer, pdb_for_view)