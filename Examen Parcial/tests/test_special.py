from src.algoritmos_alineamiento import NeedlemanWunsch, SmithWaterman
from src.visualization import print_alignment_colored
from src.utils import calculate_identity

print("="*80)
print(" TEST: CASOS ESPECIALES ".center(80, " "))
print("="*80)

nw = NeedlemanWunsch(match=1, mismatch=-1, gap=-2)
sw = SmithWaterman(match=1, mismatch=-1, gap=-2)

# CASO 1: Secuencias Idénticas
seq_identical = "ATCGATCG"
alin1, alin2, score, _ = nw.alignment(seq_identical, seq_identical)
print_alignment_colored(seq_identical, seq_identical, alin1, alin2, score,
                       "NW: Secuencias Identicas")
print(f"Score esperado: {len(seq_identical)} | Obtenido: {score} | Test: {'PASS' if score == len(seq_identical) else 'FAIL'}")

# CASO 2: Sin similitud
seq1, seq2 = "AAAA", "TTTT"
alin1, alin2, score, _ = nw.alignment(seq1, seq2)
print_alignment_colored(seq1, seq2, alin1, alin2, score, "NW: Sin Similitud")
identity = calculate_identity(alin1, alin2)
print(f"Identidad: {identity:.1f}% | Test: {'PASS' if identity == 0.0 else 'FAIL'}")

# CASO 3: Subsecuencia
seq_long, seq_short = "ACGTACGTACGT", "TACGT"
alin1_nw, alin2_nw, score_nw, _ = nw.alignment(seq_long, seq_short)
alin1_sw, alin2_sw, score_sw, _ = sw.alignment(seq_long, seq_short)

print_alignment_colored(seq_long, seq_short, alin1_nw, alin2_nw, score_nw, "NW: Subsecuencia")
print_alignment_colored(seq_long, seq_short, alin1_sw, alin2_sw, score_sw, "SW: Subsecuencia")
print(f"\nSW encontró la subsecuencia perfecta: {len(alin1_sw)} caracteres")
