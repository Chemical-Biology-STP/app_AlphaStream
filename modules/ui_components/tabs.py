import streamlit as st
from modules.ui_components.retrieve_output_ui import render_retrieve_output_tab
from modules.ui_components.submit_job_ui import render_submit_job_tab
from modules.ui_components.view_jobs_ui import render_view_jobs_tab


def render_tabs(email):
    """
    Renders the tabs in the Streamlit app and delegates UI components to their respective modules.
    """
    tab1, tab2, tab3 = st.tabs(
        ["Submit Job", "View Submitted Jobs", "Retrieve Job Output"]
    )

    with tab1:
        render_submit_job_tab(email)

    with tab2:
        render_view_jobs_tab(email)

    with tab3:
        render_retrieve_output_tab(email)
