from src.algoritmos_alineamiento import NeedlemanWunsch, SmithWaterman
from src.visualization import print_matrix, print_alignment_colored, compare_algorithms

# Crear instancias
nw = NeedlemanWunsch(match=1, mismatch=-1, gap=-2)
sw = SmithWaterman(match=1, mismatch=-1, gap=-2)

# Secuencias de prueba
seq1 = "AGRCTAR"
seq2 = "AACTA"

# Probar NW
alin1_nw, alin2_nw, score_nw, matriz_nw = nw.alignment(seq1, seq2)
print_matrix(matriz_nw, seq1, seq2, "NEEDLEMAN-WUNSCH - Matriz de Scoring")
print_alignment_colored(seq1, seq2, alin1_nw, alin2_nw, score_nw, "NEEDLEMAN-WUNSCH (Alineamiento Global)")


# Probar SW
alin1_sw, alin2_sw, score_sw, matriz_sw = sw.alignment(seq1, seq2)
print_matrix(matriz_sw, seq1, seq2, "SMITH-WATERMAN - Matriz de Scoring")
print_alignment_colored(seq1, seq2, alin1_sw, alin2_sw, score_sw, "SMITH-WATERMAN (Alineamiento Local)")

# Comparación
compare_algorithms(seq1, seq2, 
                  (alin1_nw, alin2_nw, matriz_nw, score_nw),
                  (alin1_sw, alin2_sw, matriz_sw, score_sw))
