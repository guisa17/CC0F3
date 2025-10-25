from src.algoritmos_alineamiento import NeedlemanWunsch, SmithWaterman

# Crear instancias
nw = NeedlemanWunsch(match=1, mismatch=-1, gap=-2)
sw = SmithWaterman(match=1, mismatch=-1, gap=-2)

# Secuencias de prueba
seq1 = "AGRCTAR"
seq2 = "AACTA"

# Probar NW
alin1_nw, alin2_nw, score_nw, matriz_nw = nw.alignment(seq1, seq2)
print("=== NEEDLEMAN-WUNSCH ===")
print(f"Secuencia 1: {alin1_nw}")
print(f"Secuencia 2: {alin2_nw}")
print(f"Puntuación: {score_nw}")
print()

# Probar SW
alin1_sw, alin2_sw, score_sw, matriz_sw = sw.alignment(seq1, seq2)
print("=== SMITH-WATERMAN ===")
print(f"Secuencia 1: {alin1_sw}")
print(f"Secuencia 2: {alin2_sw}")
print(f"Puntuación: {score_sw}")
