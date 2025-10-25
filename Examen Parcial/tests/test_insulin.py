from src.algoritmos_alineamiento import NeedlemanWunsch, SmithWaterman
from src.visualization import print_alignment_colored, compare_algorithms
from src.utils import load_sequence_from_file

print("="*80)
print(" TEST: INSULINA - Homo Sapiens vs Gorila ".center(80, " "))
print("="*80)

# Crear instancias
nw = NeedlemanWunsch(match=1, mismatch=-1, gap=-2)
sw = SmithWaterman(match=1, mismatch=-1, gap=-2)

# Cargar secuencias
insulin_human = load_sequence_from_file('examples/insulin-homo-sapiens.txt')
insulin_gorilla = load_sequence_from_file('examples/insulin-gorilla.txt')

print(f"\nSecuencias cargadas:")
print(f"  Homo Sapiens: {len(insulin_human)} aminoacidos")
print(f"  Gorila:       {len(insulin_gorilla)} aminoacidos")

# NEEDLEMAN-WUNSCH
alin1_nw, alin2_nw, score_nw, matriz_nw = nw.alignment(insulin_human, insulin_gorilla)
print_alignment_colored(insulin_human, insulin_gorilla, alin1_nw, alin2_nw, score_nw,
                       "NW: Insulina")

# SMITH-WATERMAN
alin1_sw, alin2_sw, score_sw, matriz_sw = sw.alignment(insulin_human, insulin_gorilla)
print_alignment_colored(insulin_human, insulin_gorilla, alin1_sw, alin2_sw, score_sw,
                       "SW: Insulina")

# COMPARACIÓN
compare_algorithms(insulin_human, insulin_gorilla,
                  (alin1_nw, alin2_nw, score_nw, matriz_nw),
                  (alin1_sw, alin2_sw, score_sw, matriz_sw))
