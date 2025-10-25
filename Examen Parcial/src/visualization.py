from .utils import calculate_identity, count_gaps, count_matches, count_mismatches
import os


def save_alignment_to_file(seq1, seq2, aligned1, aligned2, score, algorithm="", filename="output.txt", line_width=80):
    """
    Guarda el alineamiento completo en un archivo de texto con formato legible.
    Divide las líneas largas en múltiples líneas.
    """
    # Crear directorio output si no existe
    os.makedirs('output', exist_ok=True)
    filepath = os.path.join('output', filename)
    
    with open(filepath, 'w', encoding='utf-8') as f:
        # Encabezado
        f.write("="*80 + "\n")
        f.write(f"{algorithm:^80}\n")
        f.write("="*80 + "\n\n")
        
        # Información de secuencias originales
        f.write(f"Secuencia 1 original ({len(seq1)} caracteres):\n")
        f.write(f"{seq1}\n\n")
        
        f.write(f"Secuencia 2 original ({len(seq2)} caracteres):\n")
        f.write(f"{seq2}\n\n")
        
        f.write("="*80 + "\n")
        f.write("ALINEAMIENTO COMPLETO\n")
        f.write("="*80 + "\n\n")
        
        # Dividir el alineamiento en bloques de line_width caracteres
        num_blocks = (len(aligned1) + line_width - 1) // line_width
        
        for block in range(num_blocks):
            start = block * line_width
            end = min(start + line_width, len(aligned1))
            
            # Extraer secciones
            section1 = aligned1[start:end]
            section2 = aligned2[start:end]
            
            # Crear línea de matching
            match_line = ""
            for i in range(len(section1)):
                if section1[i] == section2[i] and section1[i] != '-':
                    match_line += "|"
                elif section1[i] == '-' or section2[i] == '-':
                    match_line += " "
                else:
                    match_line += "."
            
            # Escribir bloque
            f.write(f"Posición {start+1}-{end}:\n")
            f.write(f"Seq1: {section1}\n")
            f.write(f"      {match_line}\n")
            f.write(f"Seq2: {section2}\n")
            f.write("\n")
        
        # Estadísticas
        f.write("="*80 + "\n")
        f.write("ESTADÍSTICAS\n")
        f.write("="*80 + "\n\n")
        
        matches = count_matches(aligned1, aligned2)
        mismatches = count_mismatches(aligned1, aligned2)
        gaps = count_gaps(aligned1) + count_gaps(aligned2)
        identity = calculate_identity(aligned1, aligned2)
        
        f.write(f"Puntuación total:        {score}\n")
        f.write(f"Longitud alineamiento:   {len(aligned1)} caracteres\n")
        f.write(f"Matches:                 {matches}\n")
        f.write(f"Mismatches:              {mismatches}\n")
        f.write(f"Gaps totales:            {gaps}\n")
        f.write(f"Identidad:               {identity:.2f}%\n")

    return filepath


def print_matrix(matrix, seq1, seq2, title="Matriz de Scoring", max_size=15):
    """
    Imprime la matriz de forma legible con encabezados.
    Si las secuencias son muy largas (>max_size), no imprime la matriz.
    """
    if len(seq1) > max_size or len(seq2) > max_size:
        print(f"\n{'='*60}")
        print(f"{title:^60}")
        print(f"{'='*60}")
        print(f"Matriz muy grande ({len(seq1)+1}x{len(seq2)+1}) - omitida para claridad")
        print(f"{'='*60}\n")
        return
    
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


def print_alignment(seq1, seq2, aligned1, aligned2, score, algorithm="", max_display=80):
    """
    Imprime el alineamiento de forma visual con matches marcados.
    Para secuencias largas, muestra solo inicio y fin.
    """
    print(f"\n{'='*60}")
    print(f"{algorithm:^60}")
    print(f"{'='*60}\n")
    
    # Mostrar solo preview de secuencias originales si son muy largas
    if len(seq1) > max_display:
        print(f"Secuencia 1 original: {seq1[:40]}...{seq1[-40:]} ({len(seq1)} caracteres)")
    else:
        print(f"Secuencia 1 original: {seq1}")
    
    if len(seq2) > max_display:
        print(f"Secuencia 2 original: {seq2[:40]}...{seq2[-40:]} ({len(seq2)} caracteres)")
    else:
        print(f"Secuencia 2 original: {seq2}")
    
    print()
    
    # Para alineamientos largos, mostrar inicio y fin
    if len(aligned1) > max_display:
        # Inicio
        match_line_start = ""
        for i in range(40):
            if aligned1[i] == aligned2[i]:
                match_line_start += "|"
            elif aligned1[i] == '-' or aligned2[i] == '-':
                match_line_start += " "
            else:
                match_line_start += "."
        
        # Fin
        match_line_end = ""
        for i in range(len(aligned1) - 40, len(aligned1)):
            if aligned1[i] == aligned2[i]:
                match_line_end += "|"
            elif aligned1[i] == '-' or aligned2[i] == '-':
                match_line_end += " "
            else:
                match_line_end += "."
        
        print(f"Alineamiento 1: {aligned1[:40]}...{aligned1[-40:]}")
        print(f"                {match_line_start}   {match_line_end}")
        print(f"Alineamiento 2: {aligned2[:40]}...{aligned2[-40:]}")
        print(f"                (Total: {len(aligned1)} caracteres)")
    else:
        # Crear línea de matching completa
        match_line = ""
        for i in range(len(aligned1)):
            if aligned1[i] == aligned2[i]:
                match_line += "|"
            elif aligned1[i] == '-' or aligned2[i] == '-':
                match_line += " "
            else:
                match_line += "."
        
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


def print_alignment_colored(seq1, seq2, aligned1, aligned2, score, algorithm="", max_display=80, save_to_file=None):
    """
    Imprime el alineamiento con colores (requiere colorama).
    Para secuencias largas, guarda en archivo y muestra solo estadísticas.
    """
    # Si se especifica archivo o la secuencia es muy larga, guardar en archivo
    if save_to_file or len(aligned1) > max_display:
        if save_to_file is None:
            # Generar nombre automático basado en el algoritmo
            alg_name = algorithm.lower().replace(" ", "_").replace(":", "")
            save_to_file = f"{alg_name}.txt"
        
        filepath = save_alignment_to_file(seq1, seq2, aligned1, aligned2, score, algorithm, save_to_file)
        
        print(f"\n{'='*60}")
        print(f"{algorithm:^60}")
        print(f"{'='*60}")
        print(f"\nAlineamiento completo guardado en: {filepath}")
        print(f"Longitud: {len(aligned1)} caracteres")
        
        # Mostrar solo estadísticas en consola
        matches = count_matches(aligned1, aligned2)
        mismatches = count_mismatches(aligned1, aligned2)
        gaps = count_gaps(aligned1) + count_gaps(aligned2)
        identity = calculate_identity(aligned1, aligned2)
        
        print(f"\nEstadísticas:")
        print(f"  - Puntuación:  {score}")
        print(f"  - Matches:     {matches}")
        print(f"  - Mismatches:  {mismatches}")
        print(f"  - Gaps:        {gaps}")
        print(f"  - Identidad:   {identity:.1f}%")
        print(f"{'='*60}\n")
        return
    
    # Para secuencias cortas, mostrar con colores en consola
    try:
        from colorama import Fore, Style, init
        init(autoreset=True)
        
        print(f"\n{'='*60}")
        print(f"{algorithm:^60}")
        print(f"{'='*60}\n")
        
        print(f"Secuencia 1 original: {seq1}")
        print(f"Secuencia 2 original: {seq2}")
        print()
        
        # Colorear alineamiento completo
        colored_align1 = ""
        colored_align2 = ""
        match_line = ""
        
        for i in range(len(aligned1)):
            char1 = aligned1[i]
            char2 = aligned2[i]
            
            if char1 == char2:
                colored_align1 += Fore.GREEN + char1
                colored_align2 += Fore.GREEN + char2
                match_line += Fore.GREEN + "|"
            elif char1 == '-' or char2 == '-':
                colored_align1 += Fore.RED + char1
                colored_align2 += Fore.RED + char2
                match_line += " "
            else:
                colored_align1 += Fore.YELLOW + char1
                colored_align2 += Fore.YELLOW + char2
                match_line += Fore.YELLOW + "."
        
        print(f"Alineamiento 1: {colored_align1}")
        print(f"                {match_line}")
        print(f"Alineamiento 2: {colored_align2}")
        
        print()
        print(f"Puntuación total: {score}")
        
        # Estadísticas
        matches = count_matches(aligned1, aligned2)
        mismatches = count_mismatches(aligned1, aligned2)
        gaps = count_gaps(aligned1) + count_gaps(aligned2)
        identity = calculate_identity(aligned1, aligned2)
        
        print(f"\nEstadísticas:")
        print(f"  - Matches:    {Fore.GREEN}{matches}")
        print(f"  - Mismatches: {Fore.YELLOW}{mismatches}")
        print(f"  - Gaps:       {Fore.RED}{gaps}")
        print(f"  - Identidad:  {identity:.1f}%")
        print(f"  - Longitud alineamiento: {len(aligned1)} caracteres")
        print(f"{'='*60}\n")
        
    except ImportError:
        print_alignment(seq1, seq2, aligned1, aligned2, score, algorithm)


def compare_algorithms(seq1, seq2, nw_result, sw_result):
    """
    Compara los resultados de NW y SW lado a lado.
    Solo muestra estadísticas para secuencias largas.
    """
    print(f"\n{'='*80}")
    print(f"{'COMPARACIÓN: NEEDLEMAN-WUNSCH vs SMITH-WATERMAN':^80}")
    print(f"{'='*80}\n")
    
    nw_alin1, nw_alin2, nw_score, _ = nw_result
    sw_alin1, sw_alin2, sw_score, _ = sw_result
    
    # Calcular identidades
    nw_identity = calculate_identity(nw_alin1, nw_alin2)
    sw_identity = calculate_identity(sw_alin1, sw_alin2)
    
    print(f"{'Métrica':<25} | {'Needleman-Wunsch (Global)':^25} | {'Smith-Waterman (Local)':^25}")
    print("-" * 80)
    print(f"{'Longitud original':<25} | {len(seq1):^25} | {len(seq2):^25}")
    print(f"{'Longitud alineamiento':<25} | {len(nw_alin1):^25} | {len(sw_alin1):^25}")
    print(f"{'Puntuación':<25} | {nw_score:^25} | {sw_score:^25}")
    print(f"{'Identidad':<25} | {nw_identity:^24.1f}% | {sw_identity:^24.1f}%")
    print(f"{'Matches':<25} | {count_matches(nw_alin1, nw_alin2):^25} | {count_matches(sw_alin1, sw_alin2):^25}")
    print(f"{'Gaps':<25} | {count_gaps(nw_alin1) + count_gaps(nw_alin2):^25} | {count_gaps(sw_alin1) + count_gaps(sw_alin2):^25}")
    
    print(f"{'='*80}\n")