import streamlit as st

from modules.input_validation import parse_ids


def validate_nucleic_acid_sequence(sequence, valid_chars):
    """Ensure sequence contains only valid DNA or RNA characters."""
    return all(char in valid_chars for char in sequence)


def collect_nucleic_acid_sequences(entity_name, valid_chars, unique_ids):
    """
    Collects DNA or RNA sequences from user input.
    """
    st.subheader(f"🧬 Define {entity_name} Sequences")
    st.info(
        f"💡 **Multiple copies?** Separate IDs with commas (e.g., A, B, C). This implies a **homomeric {entity_name} chain**."
    )

    num_chains = st.number_input(
        f"Number of {entity_name} Chains", min_value=0, value=0, step=1
    )
    sequences = []
    valid = True

    for i in range(num_chains):
        with st.expander(f"{entity_name} {i + 1} Details", expanded=True):
            chain_id_input = st.text_input(
                f"{entity_name} ID(s)", key=f"{entity_name.lower()}_id_{i}"
            )
            chain_sequence = (
                st.text_area(
                    f"{entity_name} Sequence",
                    key=f"{entity_name.lower()}_seq_{i}",
                )
                .upper()
                .strip()
            )

            ids, error = parse_ids(chain_id_input, unique_ids, entity_name)
            if error:
                st.error(error)
                valid = False
            elif not validate_nucleic_acid_sequence(
                chain_sequence, valid_chars
            ):
                st.error(
                    f"❌ Invalid {entity_name} sequence. Only {', '.join(valid_chars)} are allowed."
                )
                valid = False
            else:
                sequences.append(
                    {
                        entity_name.lower(): {
                            "id": ids if len(ids) > 1 else ids[0],
                            "sequence": chain_sequence,
                        }
                    }
                )

    return sequences, valid
