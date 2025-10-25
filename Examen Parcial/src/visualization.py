from .utils import calculate_identity, count_gaps, count_matches, count_mismatches


def print_matrix(matrix, seq1, seq2, title="Matriz de Scoring"):
    """
    Imprime la matriz de forma legible con encabezados
    """
    print(f"\n{'='*60}")
    print(f"{title:^60}")
    print(f"{'='*60}\n")
    
    # Encabezado de columnas
    print("      ", end="")
    print("  -  ", end="")
    for char in seq2:
        print(f"  {char}  ", end="")
    print()
    
    # Primera fila (gaps)
    print("  -  ", end="")
    for j in range(len(seq2) + 1):
        print(f"{matrix[0][j]:4} ", end="")
    print()
    
    # Resto de filas
    for i in range(1, len(seq1) + 1):
        print(f"  {seq1[i-1]}  ", end="")
        for j in range(len(seq2) + 1):
            print(f"{matrix[i][j]:4} ", end="")
        print()
    print()


def print_alignment(seq1, seq2, aligned1, aligned2, score, algorithm=""):
    """
    Imprime el alineamiento de forma visual con matches marcados
    AHORA USA FUNCIONES DE utils.py
    """
    print(f"\n{'='*60}")
    print(f"{algorithm:^60}")
    print(f"{'='*60}\n")
    
    print(f"Secuencia 1 original: {seq1}")
    print(f"Secuencia 2 original: {seq2}")
    print()
    
    # Crear línea de matching
    match_line = ""
    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i]:
            match_line += "|"  # Match
        elif aligned1[i] == '-' or aligned2[i] == '-':
            match_line += " "  # Gap
        else:
            match_line += "."  # Mismatch
    
    print(f"Alineamiento 1: {aligned1}")
    print(f"                {match_line}")
    print(f"Alineamiento 2: {aligned2}")
    print()
    print(f"Puntuación total: {score}")
    
    # Estadísticas
    matches = count_matches(aligned1, aligned2)
    mismatches = count_mismatches(aligned1, aligned2)
    gaps = count_gaps(aligned1) + count_gaps(aligned2)
    identity = calculate_identity(aligned1, aligned2)
    
    print(f"\nEstadísticas:")
    print(f"  - Matches:    {matches}")
    print(f"  - Mismatches: {mismatches}")
    print(f"  - Gaps:       {gaps}")
    print(f"  - Identidad:  {identity:.1f}%")
    print(f"{'='*60}\n")


def print_alignment_colored(seq1, seq2, aligned1, aligned2, score, algorithm=""):
    """
    Imprime el alineamiento con colores (requiere colorama)
    """
    try:
        from colorama import Fore, Style, init
        init(autoreset=True)
        
        print(f"\n{'='*60}")
        print(f"{algorithm:^60}")
        print(f"{'='*60}\n")
        
        print(f"Secuencia 1 original: {seq1}")
        print(f"Secuencia 2 original: {seq2}")
        print()
        
        # Colorear alineamiento
        colored_align1 = ""
        colored_align2 = ""
        match_line = ""
        
        for i in range(len(aligned1)):
            char1 = aligned1[i]
            char2 = aligned2[i]
            
            if char1 == char2:
                # Match - Verde
                colored_align1 += Fore.GREEN + char1
                colored_align2 += Fore.GREEN + char2
                match_line += Fore.GREEN + "|"
            elif char1 == '-' or char2 == '-':
                # Gap - Rojo
                colored_align1 += Fore.RED + char1
                colored_align2 += Fore.RED + char2
                match_line += " "
            else:
                # Mismatch - Amarillo
                colored_align1 += Fore.YELLOW + char1
                colored_align2 += Fore.YELLOW + char2
                match_line += Fore.YELLOW + "."
        
        print(f"Alineamiento 1: {colored_align1}")
        print(f"                {match_line}")
        print(f"Alineamiento 2: {colored_align2}")
        print()
        print(f"Puntuación total: {score}")
        
        # Estadísticas USANDO utils.py
        matches = count_matches(aligned1, aligned2)
        mismatches = count_mismatches(aligned1, aligned2)
        gaps = count_gaps(aligned1) + count_gaps(aligned2)
        identity = calculate_identity(aligned1, aligned2)
        
        print(f"\nEstadísticas:")
        print(f"  - Matches:    {Fore.GREEN}{matches}")
        print(f"  - Mismatches: {Fore.YELLOW}{mismatches}")
        print(f"  - Gaps:       {Fore.RED}{gaps}")
        print(f"  - Identidad:  {identity:.1f}%")
        print(f"{'='*60}\n")
        
    except ImportError:
        print_alignment(seq1, seq2, aligned1, aligned2, score, algorithm)


def compare_algorithms(seq1, seq2, nw_result, sw_result):
    """
    Compara los resultados de NW y SW lado a lado
    """
    print(f"\n{'='*80}")
    print(f"{'COMPARACIÓN: NEEDLEMAN-WUNSCH vs SMITH-WATERMAN':^80}")
    print(f"{'='*80}\n")
    
    print(f"Secuencia 1: {seq1}")
    print(f"Secuencia 2: {seq2}")
    print()
    
    print(f"{'NEEDLEMAN-WUNSCH (Global)':^40} | {'SMITH-WATERMAN (Local)':^40}")
    print("-" * 80)
    
    nw_alin1, nw_alin2, _, nw_score = nw_result
    sw_alin1, sw_alin2, _, sw_score = sw_result
    
    print(f"{nw_alin1:^40} | {sw_alin1:^40}")
    print(f"{nw_alin2:^40} | {sw_alin2:^40}")
    print()
    print(f"{'Puntuación: ' + str(nw_score):^40} | {'Puntuación: ' + str(sw_score):^40}")
    print(f"{'Longitud: ' + str(len(nw_alin1)):^40} | {'Longitud: ' + str(len(sw_alin1)):^40}")
    
    # Comparar identidades
    nw_identity = calculate_identity(nw_alin1, nw_alin2)
    sw_identity = calculate_identity(sw_alin1, sw_alin2)
    print(f"{'Identidad: ' + f'{nw_identity:.1f}%':^40} | {'Identidad: ' + f'{sw_identity:.1f}%':^40}")
    
    print(f"{'='*80}\n")