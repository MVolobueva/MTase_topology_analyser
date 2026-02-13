import streamlit as st
import os
import tempfile
import numpy as np
from analyzer import MTaseAnalyzer
from utils.helpers import download_structure, parse_uploaded_file

# Конфигурация страницы
st.set_page_config(
    page_title="MTase Topology Analyzer",
    page_icon="🧬",
    layout="wide"
)

# Заголовок
st.title("🧬 **MTase Topology Analyzer**")
st.markdown("---")

# Инициализация session state
if 'analyzer' not in st.session_state:
    st.session_state.analyzer = None
if 'motifs' not in st.session_state:
    st.session_state.motifs = []
if 'results' not in st.session_state:
    st.session_state.results = {}
if 'selected_chain' not in st.session_state:
    st.session_state.selected_chain = None

# =====================================================================
# Боковая панель - Ввод данных
# =====================================================================
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
        index=0
    )
    
    st.markdown("---")
    
    # Выбор источника структуры
    st.subheader("Structure Input")
    input_type = st.radio(
        "Choose input method:",
        ["PDB ID", "AlphaFold ID", "Upload PDB file"]
    )
    
    structure_source = None
    pdb_id = None
    
    if input_type == "PDB ID":
        pdb_id = st.text_input(
            "Enter PDB ID (e.g., 4XQK):",
            value=preset_mtases[selected_preset] if selected_preset != "None" and preset_mtases[selected_preset] != "custom" else "",
            placeholder="4XQK"
        ).upper()
        if pdb_id:
            structure_source = f"PDB:{pdb_id}"
            
    elif input_type == "AlphaFold ID":
        uniprot_id = st.text_input(
            "Enter UniProt ID (e.g., P04392):",
            placeholder="P04392"
        ).upper()
        if uniprot_id:
            structure_source = f"AF:{uniprot_id}"
            
    else:  # Upload PDB file
        uploaded_file = st.file_uploader(
            "Upload PDB file:",
            type=['pdb', 'ent', 'cif']
        )
        if uploaded_file:
            structure_source = "upload"
    
    st.markdown("---")
    
    # Параметры анализа
    st.subheader("Analysis Parameters")
    
    # Каталитические мотивы
    motif_options = {
        "Default (DPPY, NPPY, PC, PS)": "default",
        "Custom motifs": "custom"
    }
    
    motif_choice = st.selectbox(
        "Catalytic motifs to search:",
        list(motif_options.keys())
    )
    
    custom_motifs = []
    if motif_choice == "Custom motifs":
        custom_motifs = st.text_area(
            "Enter motifs (one per line, regex format):",
            placeholder="[SND]PP[YFW]\nP[CS]"
        ).strip().split('\n')
        custom_motifs = [m.strip() for m in custom_motifs if m.strip()]
    
    st.markdown("---")
    
    # Кнопка запуска
    run_button = st.button("🔬 **Run Analysis**", type="primary", use_container_width=True)

# =====================================================================
# Основная область - Результаты
# =====================================================================
main_col = st.columns([2, 1])[0] if st.session_state.results else st.container()

with main_col:
    if run_button and structure_source:
        with st.spinner("🔄 Loading structure and running analysis..."):
            try:
                # Создаем анализатор
                analyzer = MTaseAnalyzer()
                
                # Загружаем структуру
                if input_type == "PDB ID" and pdb_id:
                    dssp_file = download_structure(pdb_id, source='pdb')
                elif input_type == "AlphaFold ID" and uniprot_id:
                    dssp_file = download_structure(uniprot_id, source='alphafold')
                elif input_type == "Upload PDB file" and uploaded_file:
                    dssp_file = parse_uploaded_file(uploaded_file)
                else:
                    st.error("Please provide structure input")
                    st.stop()
                
                # Загружаем DSSP
                if not analyzer.load_dssp(dssp_file):
                    st.error("Failed to load DSSP file. Make sure DSSP is installed.")
                    st.stop()
                
                # Находим вторичные структуры
                analyzer.find_all_strands()
                analyzer.build_sheet_adjacency()
                
                # Находим мотивы
                if motif_choice == "Default (DPPY, NPPY, PC, PS)":
                    motifs = analyzer.find_all_motifs()  # использует дефолтные паттерны
                else:
                    # Временно заменяем паттерны
                    original_patterns = analyzer.MOTIF_PATTERNS
                    analyzer.MOTIF_PATTERNS = custom_motifs
                    motifs = analyzer.find_all_motifs()
                    analyzer.MOTIF_PATTERNS = original_patterns
                
                st.session_state.analyzer = analyzer
                st.session_state.motifs = motifs
                
                # Очищаем предыдущие результаты
                st.session_state.results = {}
                st.session_state.selected_chain = None
                
                # Анализируем каждый мотив
                for motif in motifs:
                    chain = motif['chain']
                    result = analyzer.analyze_topology(motif_data=motif)
                    if result:
                        st.session_state.results[chain] = {
                            'motif': motif,
                            'result': result
                        }
                
                st.success(f"✅ Analysis complete! Found {len(motifs)} motifs in {len(st.session_state.results)} chains")
                
            except Exception as e:
                st.error(f"Error during analysis: {str(e)}")
                st.stop()
    
    # =================================================================
    # Отображение результатов
    # =================================================================
    if st.session_state.results:
        st.header("📊 **Analysis Results**")
        
        # Сводка по цепям
        chains_info = []
        for chain, data in st.session_state.results.items():
            motif = data['motif']
            result = data['result']
            n_strands = len(result['full_path'])
            chains_info.append({
                'Chain': chain,
                'Motif': f"{motif['text']} ({motif['res']})",
                'Strands in sheet': n_strands,
                'S4': f"{result['s4_start']}-{result['s4_end']}"
            })
        
        # Таблица со сводкой
        st.subheader("Detected catalytic motifs")
        import pandas as pd
        df = pd.DataFrame(chains_info)
        st.dataframe(df, use_container_width=True, hide_index=True)
        
        # Выбор цепи для детального просмотра
        st.subheader("🔍 **Select chain for detailed view**")
        
        chain_options = [f"Chain {info['Chain']}: {info['Motif']}" for info in chains_info]
        selected_idx = st.selectbox(
            "Choose chain:",
            range(len(chain_options)),
            format_func=lambda i: chain_options[i],
            key="chain_selector"
        )
        
        selected_chain = chains_info[selected_idx]['Chain']
        
        if selected_chain != st.session_state.selected_chain:
            st.session_state.selected_chain = selected_chain
            st.rerun()
        
        # =================================================================
        # Детальный просмотр выбранной цепи
        # =================================================================
        if st.session_state.selected_chain:
            chain_data = st.session_state.results[st.session_state.selected_chain]
            motif = chain_data['motif']
            result = chain_data['result']
            
            st.markdown("---")
            st.header(f"🧬 **Chain {st.session_state.selected_chain}**")
            st.subheader(f"Motif: {motif['text']} ({motif['res']}) | S4: {result['s4_start']}-{result['s4_end']}")
            
            # Вкладки для разных видов
            tab1, tab2, tab3 = st.tabs(["📝 Linear Topology", "🖼️ 2D Topology", "🔬 3D Structure"])
            
            with tab1:
                st.subheader("Linear Topology (N → C)")
                
                # Получаем линейную топологию
                import io
                import sys
                
                old_stdout = sys.stdout
                new_stdout = io.StringIO()
                sys.stdout = new_stdout
                
                st.session_state.analyzer.print_linear_topology_from_result(result)
                
                sys.stdout = old_stdout
                linear_topology = new_stdout.getvalue()
                
                st.text(linear_topology)
            
            with tab2:
                st.subheader("2D Topology Diagram")
                
                # Создаем 2D визуализацию
                fig = st.session_state.analyzer.visualize_topology_from_analysis(result)
                if fig:
                    st.pyplot(fig)
                else:
                    st.warning("Could not generate 2D topology")
            
            with tab3:
                st.subheader("3D Structure Visualization")
                
                pdb_for_view = pdb_id if input_type == "PDB ID" else None
                
                # Получаем 3D визуализацию
                view = st.session_state.analyzer.visualize_3d_structure(
                    result, 
                    pdb_id=pdb_for_view,
                    chain=motif['chain']  # ← ВОТ ТАК ПРАВИЛЬНО!
                )
                
                if view:
                    # Сохраняем и отображаем HTML
                    import tempfile
                    with tempfile.NamedTemporaryFile(suffix='.html', delete=False) as f:
                        view.write_html(f.name)
                        with open(f.name, 'r') as f_html:
                            html_content = f_html.read()
                        st.components.v1.html(html_content, height=600)
                    os.unlink(f.name)
                else:
                    st.warning("Could not generate 3D visualization")
            
            # Дополнительная информация
            with st.expander("📋 Detailed strand information"):
                # Здесь можно вывести таблицу из analyze_topology
                # Она уже выводится в консоль, нужно перехватить
                pass

# =====================================================================
# Информация и инструкция
# =====================================================================
with st.expander("ℹ️ **About & Instructions**"):
    st.markdown("""
    ### **MTase Topology Analyzer**
    
    This tool analyzes the topology of MTase catalytic domains based on DSSP data.
    
    **How to use:**
    1. **Select a pre-studied MTase** from the dropdown, or
    2. **Input your own structure** via PDB ID, AlphaFold ID, or upload a PDB file
    3. **Choose catalytic motifs** to search (default or custom)
    4. Click **Run Analysis**
    5. **Select a chain** from the results table
    6. **Explore** linear topology, 2D diagram, and 3D structure
    
    **Interpretation:**
    - **S4** - catalytic strand containing the motif
    - **S3, S2, S1** - strands on C-terminal side
    - **S5, S6, S7** - strands on N-terminal side
    - **Hu** - helices above the sheet (green)
    - **Hd** - helices below the sheet (red)
    
    **References:**
    - Based on the algorithm for MTase topology analysis
    - Uses DSSP for secondary structure assignment
    """)

# Footer
st.markdown("---")
st.markdown("🧬 **MTase Topology Analyzer** | Created for research purposes")