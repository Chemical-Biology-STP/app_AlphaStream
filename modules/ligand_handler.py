import streamlit as st
from rdkit import Chem

from modules.input_validation import parse_ids


def collect_ligands(unique_ids):
    """
    Collects ligand information.
    """
    st.subheader("💊 Define Ligands")
    st.info(
        "💡 **Multiple copies?** Separate IDs with commas (e.g., A, B, C). This implies a **ligand with multiple copies**."
    )

    num_ligands = st.number_input(
        "Number of Ligands", min_value=0, value=0, step=1
    )
    ligands = []
    ligands_valid = True

    for i in range(num_ligands):
        with st.expander(f"Ligand {i + 1} Details", expanded=True):
            ligand_id_input = st.text_input(
                "Ligand ID(s)", key=f"ligand_id_{i}"
            )
            smiles = st.text_input(
                "SMILES (optional)", key=f"ligand_smiles_{i}"
            ).strip()

            ids, error = parse_ids(ligand_id_input, unique_ids, "Ligand")
            if error:
                st.error(error)
                ligands_valid = False
            elif smiles and not Chem.MolFromSmiles(smiles):
                st.error("❌ Invalid SMILES format. Please check your input.")
                ligands_valid = False
            else:
                ligands.append(
                    {
                        "ligand": {
                            "id": ids if len(ids) > 1 else ids[0],
                            "smiles": smiles,
                        }
                    }
                )

    return ligands, ligands_valid
