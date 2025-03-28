import re

from rdkit import Chem

# Allowed characters in a SMILES string based on your table
ALLOWED_SMILES_CHARS = r"[A-Za-z0-9\-=#:().,+*\[\]~]"


def is_valid_smiles(smiles: str):
    """
    Validates if a given SMILES string is valid using RDKit and allowed characters.

    Steps:
    1. Checks if the input contains any **invalid** characters.
    2. Uses RDKit to verify if the SMILES is correctly formatted.

    Returns:
        (bool, str): True if valid, otherwise False with an error message.
    """
    if not smiles:
        return False, "SMILES cannot be empty."

    # Check for invalid characters
    invalid_chars = re.findall(rf"[^{ALLOWED_SMILES_CHARS}]", smiles)
    if invalid_chars:
        unique_invalid_chars = "".join(set(invalid_chars))
        return (
            False,
            f"Invalid characters detected in SMILES: {unique_invalid_chars}",
        )

    # RDKit validation
    if not Chem.MolFromSmiles(smiles):
        return False, "Invalid SMILES format. Please check your input."

    return True, ""


def sanitize_input(input_text: str):
    """
    Removes non-alphabetic characters from the input text.
    """
    return re.sub(r"[^A-Za-z]", "", input_text)


def is_valid_string(s: str):
    """
    Checks if the string contains only alphabets, numbers, dashes, and underscores.
    If invalid characters are found (including spaces, newlines, tabs, etc.), returns the invalid characters.
    """
    invalid_chars = set(
        re.findall(r"[^a-zA-Z0-9_-]", s)
    )  # Use set to remove duplicates

    if invalid_chars:
        invalid_chars_display = " ".join(
            repr(char) for char in invalid_chars
        )  # Show special characters clearly
        return False, f"Invalid characters: {invalid_chars_display}"
    return True, ""


def is_valid_email(email: str):
    """
    Validates if the email is a Crick email.
    """
    return email.endswith("@crick.ac.uk")


def is_unique_id(entity_id, unique_ids):
    """
    Ensures that the entity ID is unique and follows expected rules.
    - Must be a single uppercase alphabet letter.
    - Must not be a duplicate.
    """
    entity_id = sanitize_input(entity_id).upper()

#    if len(entity_id) != 1 or not entity_id.isalpha():
#        return False, "❌ ID must be a single uppercase alphabet letter."

    if entity_id in unique_ids:
        return (
            False,
            f"❌ Duplicate ID detected: {entity_id}. Please enter a unique ID.",
        )

    unique_ids.add(entity_id)
    return True, ""


def parse_ids(id_input, unique_ids, entity_name):
    """
    Parses user input for IDs, ensuring they are:
    - **Single uppercase alphabetic characters** (A-Z)
    - **Comma-separated**
    - **Unique**

    If duplicates exist, lists all duplicate IDs in the error message.

    Returns:
        tuple: (list of valid IDs, error message)
    """
    ids = [x.strip().upper() for x in id_input.split(",") if x.strip()]

    if not ids:
        return None, f"❌ {entity_name} ID cannot be empty."

    # Detect invalid IDs (must be single uppercase letters)
    invalid_ids = [x for x in ids if not re.fullmatch(r"^[A-Z]+$", x)]
    if invalid_ids:
        return None, f"❌ Invalid {entity_name} IDs: {', '.join(invalid_ids)}"

    # Detect duplicates in the current input
    input_duplicates = {id for id in ids if ids.count(id) > 1}
    if input_duplicates:
        return (
            None,
            f"❌ Duplicate {entity_name} IDs in input: {', '.join(input_duplicates)}",
        )

    # Detect duplicates against existing unique IDs
    existing_duplicates = {id for id in ids if id in unique_ids}
    if existing_duplicates:
        return (
            None,
            f"❌ The following {entity_name} IDs are already used: {', '.join(existing_duplicates)}",
        )

    # Add unique IDs to the tracking set
    unique_ids.update(ids)
    return ids, None


def sanitize_protein_sequence(sequence):
    """
    Extracts only alphabets from a given sequence and converts them to uppercase.

    Args:
        sequence (str): The raw input protein sequence.

    Returns:
        str: Sanitized uppercase protein sequence.
    """
    return re.sub(r"[^A-Za-z]", "", sequence).upper()


def validate_nucleic_acid_sequence(sequence, valid_chars):
    """Ensure sequence contains only valid DNA or RNA characters."""
    return all(char in valid_chars for char in sequence)
