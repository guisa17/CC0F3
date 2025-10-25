from src.algoritmos_alineamiento import NeedlemanWunsch, SmithWaterman
from src.visualization import print_alignment_colored, compare_algorithms
from src.utils import load_sequence_from_file

print("="*80)
print(" TEST: HEMOGLOBINA - Homo Sapiens vs Conejo ".center(80, " "))
print("="*80)

# Crear instancias
nw = NeedlemanWunsch(match=1, mismatch=-1, gap=-2)
sw = SmithWaterman(match=1, mismatch=-1, gap=-2)

# Cargar secuencias
hemo_human = load_sequence_from_file('examples/hemoglobin-homo-sapiens.txt')
hemo_rabbit = load_sequence_from_file('examples/hemoglobin-rabbit.txt')

print(f"\nSecuencias cargadas:")
print(f"  Homo Sapiens: {len(hemo_human)} aminoacidos")
print(f"  Conejo:       {len(hemo_rabbit)} aminoacidos")

# NEEDLEMAN-WUNSCH
alin1_nw, alin2_nw, score_nw, matriz_nw = nw.alignment(hemo_human, hemo_rabbit)
print_alignment_colored(hemo_human, hemo_rabbit, alin1_nw, alin2_nw, score_nw,
                       "NW: Hemoglobina")

# SMITH-WATERMAN
alin1_sw, alin2_sw, score_sw, matriz_sw = sw.alignment(hemo_human, hemo_rabbit)
print_alignment_colored(hemo_human, hemo_rabbit, alin1_sw, alin2_sw, score_sw,
                       "SW: Hemoglobina")

# COMPARACIÓN
compare_algorithms(hemo_human, hemo_rabbit,
                  (alin1_nw, alin2_nw, score_nw, matriz_nw),
                  (alin1_sw, alin2_sw, score_sw, matriz_sw))
