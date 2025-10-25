def calculate_identity(aligned1, aligned2):
    """
    Calcula el porcentaje de identidad entre dos secuencias alineadas
    """
    if len(aligned1) != len(aligned2):
        raise ValueError("Las secuencias alineadas deben tener la misma longitud")
    
    matches = sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != '-')
    total = len(aligned1)
    
    return (matches / total) * 100 if total > 0 else 0


def count_gaps(aligned_seq):
    """
    Cuenta el número de gaps en una secuencia alineada
    """
    return aligned_seq.count('-')


def count_matches(aligned1, aligned2):
    """
    Cuenta el número de matches entre dos secuencias alineadas
    """
    return sum(1 for a, b in zip(aligned1, aligned2) if a == b and a != '-')


def count_mismatches(aligned1, aligned2):
    """
    Cuenta el número de mismatches entre dos secuencias alineadas
    """
    return sum(1 for a, b in zip(aligned1, aligned2) if a != b and a != '-' and b != '-')


def validate_dna(sequence):
    """
    Valida que una secuencia sea ADN válido (solo A, C, G, T)
    """
    valid_chars = set('ACGT')
    return all(char.upper() in valid_chars for char in sequence)


def validate_protein(sequence):
    """
    Valida que una secuencia sea proteína válida (aminoácidos)
    """
    valid_chars = set('ACDEFGHIKLMNPQRSTVWY')
    return all(char.upper() in valid_chars for char in sequence)


def load_sequence_from_file(filepath):
    """
    Carga una secuencia desde un archivo de texto
    """
    try:
        with open(filepath, 'r') as f:
            # Leer y limpiar la secuencia (remover espacios y saltos de línea)
            sequence = ''.join(line.strip() for line in f)
            return sequence.upper()
    except FileNotFoundError:
        raise FileNotFoundError(f"No se encontró el archivo: {filepath}")
