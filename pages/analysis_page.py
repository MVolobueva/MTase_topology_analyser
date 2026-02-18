import streamlit as st
from analyzer import MTaseAnalyzer
from utils.helpers import download_structure, parse_uploaded_file
from components.visualizations import show_linear_topology, show_2d_topology, show_3d_topology
import pandas as pd
import os

def show():
    st.title("MTase Topology Analysis")
    
    # Проверяем, что нужно запустить анализ
    run_button = st.session_state.get('run_button', False)
    
    # Определяем источник структуры
    structure_source = None
    source_type = None
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
        # Прогресс-бар для отслеживания шагов
        progress_bar = st.progress(0, text="Starting analysis...")
        
        with st.spinner("🔄 Loading structure and running analysis..."):
            try:
                # Шаг 1: Создаем анализатор
                progress_bar.progress(10, text="Initializing analyzer...")
                analyzer = MTaseAnalyzer()
                
                # Шаг 2: Загружаем структуру
                progress_bar.progress(20, text="Downloading structure...")
                if source_type == 'pdb':
                    result_files = download_structure(structure_source, source='pdb')
                    dssp_file = result_files['dssp']
                    st.session_state.current_pdb_file = result_files['pdb']
                    st.session_state.structure_source = 'pdb'
                    
                elif source_type == 'alphafold':
                    result_files = download_structure(structure_source, source='alphafold')
                    if result_files is None:
                        st.error(f"Failed to download AlphaFold structure for {structure_source}")
                        return
                    dssp_file = result_files['dssp']
                    st.session_state.current_pdb_file = result_files['pdb']
                    st.session_state.structure_source = 'alphafold'
                    
                elif source_type == 'upload':
                    result_files = parse_uploaded_file(structure_source)
                    if result_files is None:
                        st.error("Failed to parse uploaded PDB file")
                        return
                    dssp_file = result_files['dssp']
                    st.session_state.current_pdb_file = result_files['pdb']
                    st.session_state.structure_source = 'upload'
                
                if not dssp_file:
                    st.error("Failed to load structure")
                    return
                
                # Шаг 3: Загружаем DSSP
                progress_bar.progress(40, text="Loading DSSP data...")
                if not analyzer.load_dssp(dssp_file):
                    st.error("Failed to load DSSP file. Make sure DSSP is installed.")
                    return
                
                # Кнопка для скачивания DSSP файла
                with open(dssp_file, 'r') as f:
                    dssp_content = f.read()
                st.download_button(
                    label="📥 Download DSSP file",
                    data=dssp_content,
                    file_name=f"{structure_source}.dssp",
                    mime="text/plain"
                )
                
                # Шаг 4: Находим вторичные структуры
                progress_bar.progress(60, text="Identifying secondary structures...")
                analyzer.find_all_strands()
                analyzer.build_sheet_adjacency()
                
                # Шаг 5: Находим мотивы
                progress_bar.progress(80, text="Searching catalytic motifs...")
                if st.session_state.get('motif_choice') == "Default (DPPY, NPPY, PC, PS)":
                    motifs = analyzer.find_all_motifs()
                else:
                    # Получаем пользовательские мотивы
                    custom_motifs_text = st.session_state.get('custom_motifs', '')
                    if custom_motifs_text.strip():
                        custom_patterns = [m.strip() for m in custom_motifs_text.split('\n') if m.strip()]
                        motifs = analyzer.find_all_motifs(custom_patterns=custom_patterns)
                    else:
                        motifs = analyzer.find_all_motifs()
                        st.warning("No custom motifs entered, using default patterns")
                
                # Фильтруем мотивы
                motifs = analyzer.filter_motifs_by_topology(motifs)
                
                # Шаг 6: Анализируем каждый мотив
                progress_bar.progress(90, text="Analyzing topology...")
                results = {}
                for motif in motifs:
                    chain = motif['chain']
                    motif_res = motif['res']
                    motif_text = motif['text']
                    key = f"{chain}_{motif_text}_{motif_res}"
                    
                    result = analyzer.analyze_topology(motif_data=motif)
                    
                    if result:
                        result['motif'] = motif
                        result['motif_text'] = motif_text
                        result['motif_res'] = motif_res
                        
                        # Получаем последовательность тяжей для отображения
                        strand_sequence = []
                        # Сортируем full_path по порядку в цепи
                        sorted_indices = sorted(result['full_path'], key=lambda idx: result['strands'][idx][0] if result['strands'][idx] else 0)
                        
                        for idx in sorted_indices:
                            # Ищем имя тяжа по его индексу
                            for name, i in result['strand_names'].items():
                                if i == idx:
                                    strand_sequence.append(name)
                                    break
                        result['strand_sequence'] = ' → '.join(strand_sequence)
                        
                        # Получаем координаты листа (первый и последний тяж)
                        if result['full_path']:
                            first_idx = result['full_path'][0]
                            last_idx = result['full_path'][-1]
                            first_strand = result['strands'][first_idx]
                            last_strand = result['strands'][last_idx]
                            result['sheet_start'] = analyzer._get_res_num(first_strand[0])
                            result['sheet_end'] = analyzer._get_res_num(last_strand[-1])
                        else:
                            result['sheet_start'] = 0
                            result['sheet_end'] = 0
                        
                        results[key] = {
                            'motif': motif,
                            'result': result,
                            'display_chain': chain,
                            'display_motif': f"{motif_text}",
                            'strand_sequence': result['strand_sequence'],
                            'sheet_range': f"{result['sheet_start']}-{result['sheet_end']}",
                            'n_strands': len(result['full_path'])
                        }
                
                progress_bar.progress(100, text="Analysis complete!")
                st.session_state.analyzer = analyzer
                st.session_state.results = results
                st.success(f"✅ Analysis complete! Found {len(motifs)} motifs in {len(results)} chains")
                
                # Убираем прогресс-бар после завершения
                progress_bar.empty()
                
            except Exception as e:
                st.error(f"Error during analysis: {str(e)}")
                import traceback
                st.code(traceback.format_exc())
                progress_bar.empty()
    
    # Отображаем результаты
    if st.session_state.get('results'):
        st.header("📊 **Analysis Results**")
        
        # Создаем улучшенную таблицу
    
        # Создаем улучшенную таблицу
        table_data = []
        for key, data in st.session_state.results.items():
            result = data['result']
            motif = data['motif']
            
            # Диапазон мотива
            motif_text = motif['text']
            motif_length = len(motif_text)
            motif_start = motif['res']
            motif_end = motif_start + motif_length - 1
            motif_coords = f"{motif_start}-{motif_end}"
            
            # ✅ ПОЛУЧАЕМ ПОСЛЕДОВАТЕЛЬНОСТЬ ТЯЖЕЙ В ПРАВИЛЬНОМ ПОРЯДКЕ
            strand_dict = {}
            for idx in result['full_path']:
                strand_name = result['strand_names'][idx]
                # Получаем первый остаток тяжа
                first_residue = result['strands'][idx][0]
                # Преобразуем в число
                start_num = st.session_state.analyzer._get_res_num(first_residue)
                strand_dict[start_num] = strand_name
            
            # Сортируем по номеру старта
            strand_sequence = [strand_dict[start] for start in sorted(strand_dict.keys())]
            strand_sequence_str = ' → '.join(strand_sequence)
            
            # Для отладки
            print(f"Старты: {strand_dict}")
            print(f"Последовательность: {strand_sequence_str}")
            
            # Получаем координаты листа
            if result['full_path']:
                first_idx = result['full_path'][0]
                last_idx = result['full_path'][-1]
                first_strand = result['strands'][first_idx]
                last_strand = result['strands'][last_idx]
                sheet_start = st.session_state.analyzer._get_res_num(first_strand[0])
                sheet_end = st.session_state.analyzer._get_res_num(last_strand[-1])
                sheet_range = f"{sheet_start}-{sheet_end}"
            else:
                sheet_range = "N/A"
            
            table_data.append({
                'Chain': data['display_chain'],
                'Motif': motif_text,
                'Motif Position': motif_coords,
                'Strands in sheet': len(result['full_path']),
                'Sheet sequence': strand_sequence_str,
                'Sheet range': sheet_range
            })    
        df = pd.DataFrame(table_data)
        st.dataframe(df, use_container_width=True, hide_index=True)
        
        # Выбор мотива
        if table_data:
            motif_options = [f"{row['Chain']}: {row['Motif']} (sheet: {row['Sheet sequence']})" for row in table_data]
            selected_idx = st.selectbox(
                "Select motif for detailed view:",
                range(len(motif_options)),
                format_func=lambda i: motif_options[i],
                key="motif_selector_detail"
            )
            
            selected_key = list(st.session_state.results.keys())[selected_idx]
            data = st.session_state.results[selected_key]
            
            st.markdown("---")
            st.header(f"🧬 **Chain {data['display_chain']}**")
            
            # Дополнительная информация о листе
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Motif", data['display_motif'])
            with col2:
                st.metric("Strands in sheet", data['n_strands'])
            with col3:
                st.metric("Sheet range", data['sheet_range'])
            
            # Вкладки - только 3!
            tab1, tab2, tab3 = st.tabs([
                "Linear Topology", 
                "2D Topology", 
                "3D Structure"
            ])

            with tab1:
                st.subheader("Linear Topology")
                
                # Две колонки внутри первой вкладки
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**N→C Sequence of Strands**")
                    
                    # Получаем последовательность тяжей с координатами
                    strand_positions = []
                    for idx in result['full_path']:
                        strand_name = result['strand_names'][idx]
                        # Получаем первый и последний остаток тяжа
                        first_res = result['strands'][idx][0]
                        last_res = result['strands'][idx][-1]
                        start_num = st.session_state.analyzer._get_res_num(first_res)
                        end_num = st.session_state.analyzer._get_res_num(last_res)
                        strand_positions.append((start_num, f"{strand_name}[{start_num}-{end_num}]"))
                    
                    # Сортируем по start_num
                    strand_positions.sort(key=lambda x: x[0])
                    
                    # Формируем строку с проверкой на большие разрывы
                    sequence_parts = []
                    for i, (start, part) in enumerate(strand_positions):
                        sequence_parts.append(part)
                        
                        # Проверяем разрыв до следующего
                        if i < len(strand_positions) - 1:
                            next_start = strand_positions[i+1][0]
                            if next_start - start > 50:
                                sequence_parts.append("🔴 **BIG LOOP** 🔴")
                    
                    strand_sequence_str = ' — '.join(sequence_parts)
                    st.info(f"**{strand_sequence_str}**")
                    st.caption("Order of beta-strands from N to C terminus. 🔴 BIG LOOP indicates insertion >50 amino acids between consecutive strands.")
                
                with col2:
                    st.markdown("**Linear Topology with All Elements**")
                    show_linear_topology(data['result'], st.session_state.analyzer)

            with tab2:
                st.subheader("2D Topology Diagram")
                fig = st.session_state.analyzer.visualize_topology_interactive(data['result'])
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning("Could not generate interactive 2D topology")
            
            with tab3:
                st.subheader("3D Structure Visualization")
                pdb_for_view = st.session_state.pdb_id if st.session_state.get('input_type') == "PDB ID" else None
                show_3d_topology(
                    data['result'], 
                    st.session_state.analyzer, 
                    pdb_id=pdb_for_view,
                    pdb_file=st.session_state.get('current_pdb_file')
                )