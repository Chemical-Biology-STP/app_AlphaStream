import json
import re

import streamlit as st
from modules.input_validation import parse_ids, sanitize_protein_sequence

# ✅ Load CCD.json (Modify path if necessary)
#with open("CCD.json", "r") as f:
with open("app_AlphaStream/CCD.json", "r") as f:
    CCD_DATA = json.load(f)

# ✅ Base URL for CCD modification information
CCD_BASE_URL = "https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/"


def collect_protein_sequences(unique_ids):
    """
    Collects and validates protein sequences along with PTMs.
    """
    st.subheader("🧬 Define Protein Sequences")
    st.info(
        "💡 **Multiple copies?** Separate IDs with commas (e.g., `A, B, C`). This implies a **homomeric protein chain**."
    )

    num_proteins = st.number_input(
        "Number of Protein Entities", min_value=0, value=0, step=1
    )
    protein_sequences = []
    proteins_valid = True

    for i in range(num_proteins):
        with st.expander(f"Protein {i + 1} Details", expanded=True):
            protein_id_input = st.text_input("Protein ID(s)", key=f"protein_id_{i}")
            raw_sequence = st.text_area("Protein Sequence", key=f"protein_seq_{i}")

            ids, error = parse_ids(protein_id_input, unique_ids, "Protein")
            if error:
                st.error(error)
                proteins_valid = False
            elif not raw_sequence:
                st.error("❌ Protein sequence cannot be empty.")
                proteins_valid = False
            else:
                sanitized_sequence = sanitize_protein_sequence(raw_sequence)

                if not sanitized_sequence:
                    st.error(
                        "❌ Invalid sequence! No valid amino acid characters found."
                    )
                    proteins_valid = False
                    continue

                # ✅ Display the collapsible amino acid position table
                display_sequence_positions(sanitized_sequence)

                modifications = collect_protein_modifications(sanitized_sequence, i)

                protein_sequences.append(
                    {
                        "protein": {
                            "id": ids if len(ids) > 1 else ids[0],
                            "sequence": sanitized_sequence,
                            "modifications": (modifications if modifications else []),
                        }
                    }
                )

    return protein_sequences, proteins_valid


def display_sequence_positions(cleaned_sequence):
    """
    Displays the amino acid sequence along with their positions in rows of 10.
    """
    if cleaned_sequence:
        chunk_size = 10  # ✅ Number of amino acids per row
        num_chunks = (
            len(cleaned_sequence) + chunk_size - 1
        ) // chunk_size  # Total number of rows

        table_rows = ""  # Store table HTML rows

        for chunk_index in range(num_chunks):
            start = chunk_index * chunk_size
            end = start + chunk_size

            # Extract a 10-residue chunk
            chunk_positions = cleaned_sequence[start:end]
            chunk_numbers = range(start + 1, min(end + 1, len(cleaned_sequence) + 1))

            # ✅ Format positions and amino acids in separate rows
            positions_html = "".join(
                f"<td style='text-align: center;'>{pos}</td>" for pos in chunk_numbers
            )
            sequence_html = "".join(
                f"<td style='text-align: center; font-weight: bold;'>{aa}</td>"
                for aa in chunk_positions
            )

            table_rows += f"<tr>{positions_html}</tr><tr>{sequence_html}</tr><tr style='height: 10px;'></tr>"  # Add spacing row

        # ✅ Generate full table with row spacing inside an HTML expander
        expander_html = f"""
        <details>
            <summary style="font-size: 18px; font-weight: bold; cursor: pointer;">📍 Show Sequence with Positions</summary>
            <table style="border-collapse: separate; border-spacing: 5px 10px; width: 100%; margin-top: 10px;">
                {table_rows}
            </table>
        </details>
        """

        # ✅ Display the HTML expander in Streamlit
        st.markdown(expander_html, unsafe_allow_html=True)


def collect_protein_modifications(sequence, protein_index):
    """
    Collects PTMs for a given protein sequence, allowing selection based on CCD.json.

    Args:
        sequence (str): The sanitized protein sequence.
        protein_index (int): Index of the protein in the input list.

    Returns:
        list[dict]: A list of PTM dictionaries containing `ptmType` and `ptmPosition`.
    """
    modifications = []
    occupied_positions = set()
    used_ptm_types = set()

    num_modifications = st.number_input(
        f"Number of Modifications for Protein {protein_index + 1}",
        min_value=0,
        value=0,
        step=1,
        key=f"num_modifications_{protein_index}",
    )

    if num_modifications > 0:
        st.markdown("### 🏷 Post-Translational Modifications (PTMs)")

        for j in range(num_modifications):
            st.subheader(f"Modification {j + 1}")

            # Select modification method
            method = st.radio(
                f"Choose method for PTM {j + 1}",
                ("Manual Entry", "Select Amino Acid"),
                key=f"ptm_method_{protein_index}_{j}",
                horizontal=True,
            )

            ptm_type = None

            if method == "Manual Entry":
                ptm_type = (
                    st.text_input(
                        f"PTM Type (CCD Code) {j + 1}",
                        key=f"ptm_type_{protein_index}_{j}",
                    )
                    .strip()
                    .upper()
                )

            else:
                # ✅ Sort amino acids alphabetically
                sorted_amino_acids = sorted(CCD_DATA.keys())

                # Select Amino Acid
                amino_acid = st.selectbox(
                    f"Select Amino Acid for PTM {j + 1}",
                    options=sorted_amino_acids,  # ✅ Now sorted alphabetically
                    key=f"aa_choice_{protein_index}_{j}",
                )

                # Select Modification from CCD.json
                ptm_type = st.selectbox(
                    f"Select Modification for {amino_acid}",
                    options=CCD_DATA.get(amino_acid, []),
                    key=f"ptm_type_select_{protein_index}_{j}",
                )

            if not ptm_type:
                st.error("❌ PTM Type cannot be empty.")
                continue

            # ✅ Generate a clickable link for the selected CCD code
            ptm_url = f"{CCD_BASE_URL}{ptm_type}"
            st.markdown(
                f"[🔗 More Info on {ptm_type}]({ptm_url})",
                unsafe_allow_html=True,
            )

            if ptm_type in used_ptm_types:
                st.error(
                    f"❌ PTM Type '{ptm_type}' has already been used for this protein."
                )
                continue

            available_positions = [
                pos
                for pos in range(1, len(sequence) + 1)
                if pos not in occupied_positions
            ]

            ptm_positions = st.multiselect(
                f"Residue Positions for {ptm_type} (1-based index)",
                options=available_positions,
                key=f"ptm_positions_{protein_index}_{j}",
            )

            if not ptm_positions:
                st.error(f"❌ No positions selected for PTM Type '{ptm_type}'.")
                continue

            for pos in ptm_positions:
                modifications.append({"ptmType": ptm_type, "ptmPosition": pos})
                occupied_positions.add(pos)

            used_ptm_types.add(ptm_type)

    return modifications
