import streamlit as st
from components.sidebar import render_sidebar
from pages import analysis_page, documentation_page

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
    analysis_page.show()
else:
    documentation_page.show()