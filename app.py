import streamlit as st
from components.sidebar import render_sidebar
from pages.analysis_page import show as show_analysis  # ✅ правильный импорт
from pages.documentation_page import show as show_docs

# Конфигурация страницы
st.set_page_config(
    page_title="MTase Topology Analyzer",
    page_icon="🧬",
    layout="wide"
)

# Заголовок
st.title("🧬 **MTase Topology Analyzer**")
st.markdown("---")

# Боковая панель
render_sidebar()

# Навигация по страницам
page = st.sidebar.radio("Navigate", ["Analysis", "Documentation"])

if page == "Analysis":
    show_analysis()  # ✅ вызываем правильную функцию
else:
    show_docs()