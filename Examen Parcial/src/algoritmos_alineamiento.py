from .utils import validate_dna, validate_protein

class AlgoritmoAlineamiento:
    def __init__(self, match=1, mismatch=-1, gap=-2):
        """
        - match: puntuación por coincidencia
        - mismatch: penalización por no coincidencia
        - gap: penalización por un espacio
        """
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
    
    def calculate_score(self, a, b):
        if a == b:
            return self.match
        else:
            return self.mismatch
    
    def validate_seqs(self, seq1, seq2, seq_type=None):
        if not seq1 or not seq2:
            raise ValueError("Las secuencias no pueden estar vacías.")

        if seq_type == "dna":
            if not validate_dna(seq1):
                raise ValueError("La secuencia 1 no es una secuencia de ADN válida.")   
            if not validate_dna(seq2):
                raise ValueError("La secuencia 2 no es una secuencia de ADN válida.")

        elif seq_type == "protein":
            if not validate_protein(seq1):
                raise ValueError("La secuencia 1 no es una secuencia de proteínas válida.")
            if not validate_protein(seq2):
                raise ValueError("La secuencia 2 no es una secuencia de proteínas válida.")

        return True


# ===============================================================================
class NeedlemanWunsch(AlgoritmoAlineamiento):
    """
    Algoritmo de Needleman-Wunsch para alineamiento global
    """

    def initialize_matrix(self, seq1, seq2):
        """
        Se llena la primera fila y columna de la matriz con penalizaciones por gaps.
        [0, -2, -4, -6, ...]
        """
        rows = len(seq1) + 1
        cols = len(seq2) + 1

        # Matriz de ceros
        matrix = [[0 for _ in range(cols)] for _ in range(rows)]

        # Inicializar primera fila y columna
        for i in range(rows):
            matrix[i][0] = self.gap * i
        for j in range(cols):
            matrix[0][j] = self.gap * j
        return matrix


    def fill_matrix(self, matrix, seq1, seq2):
        """
        Calcular el máximo de
        - gap vertical
        - gap horizontal
        - match/mismatch diagonal
        """
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                match_mismatch = matrix[i-1][j-1] + self.calculate_score(seq1[i-1], seq2[j-1])
                gap_vertical = matrix[i-1][j] + self.gap
                gap_horizontal = matrix[i][j-1] + self.gap

                matrix[i][j] = max(match_mismatch, gap_vertical, gap_horizontal)
        return matrix
    

    def traceback(self, seq1, seq2, matrix):
        """
        Reconstruir el alineamiento óptimo siguiendo el camino de la matriz
        - empezar desde la esquina inferior derecha -> (0, 0)
        """
        aligned_seq1 = []
        aligned_seq2 = []
        i, j = len(seq1), len(seq2)

        while i > 0 or j > 0:
            # Primera columan, solo hacia arriba ^
            if i > 0 and j == 0:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
            
            # Primera fila, solo hacia izquierda <
            elif j > 0 and i == 0:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1

            else:
                # Calcular de donde vino el valor actual
                score_current = matrix[i][j]
                score_diagonal = matrix[i-1][j-1] + self.calculate_score(seq1[i-1], seq2[j-1])
                score_up = matrix[i-1][j] + self.gap
                score_left = matrix[i][j-1] + self.gap

                # Seguir el camino óptimo
                if score_current == score_diagonal:
                    aligned_seq1.append(seq1[i-1])
                    aligned_seq2.append(seq2[j-1])
                    i -= 1
                    j -= 1
                elif score_current == score_up:
                    aligned_seq1.append(seq1[i-1])
                    aligned_seq2.append('-')
                    i -= 1
                else:  # score_current == score_left
                    aligned_seq1.append('-')
                    aligned_seq2.append(seq2[j-1])
                    j -= 1

        # Las secuencias se construyen al revés
        return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))
    

    def alignment(self, seq1, seq2):
        """
        Ejecuta el alineamiento completo
        """
        self.validate_seqs(seq1, seq2)

        matrix = self.initialize_matrix(seq1, seq2)
        matrix = self.fill_matrix(matrix, seq1, seq2)
        aligned_seq1, aligned_seq2 = self.traceback(seq1, seq2, matrix)
        score = matrix[len(seq1)][len(seq2)]

        return aligned_seq1, aligned_seq2, score, matrix


# ===============================================================================
class SmithWaterman(AlgoritmoAlineamiento):
    """
    Algoritmo de Smith-Waterman para alineamiento local
    """
    
    def initialize_matrix(self, seq1, seq2):
        """
        Se llena la primera fila y columna de la matriz con 0s
        """
        rows = len(seq1) + 1
        cols = len(seq2) + 1
        matrix = [[0 for _ in range(cols)] for _ in range(rows)]
        return matrix
    

    def fill_matrix(self, matrix, seq1, seq2):
        """
        Si el valor es negativo, se coloca 0
        """
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                match_mismatch = matrix[i-1][j-1] + self.calculate_score(seq1[i-1], seq2[j-1])
                gap_vertical = matrix[i-1][j] + self.gap
                gap_horizontal = matrix[i][j-1] + self.gap
                
                # máximo con 0
                matrix[i][j] = max(match_mismatch, gap_vertical, gap_horizontal, 0)
        
        return matrix
    

    def find_max(self, matrix):
        """
        Encuentra la posición del valor máximo en la matriz
        """
        max_value = 0
        max_pos = (0, 0)
        
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if matrix[i][j] > max_value:
                    max_value = matrix[i][j]
                    max_pos = (i, j)
        
        return max_pos, max_value
    

    def traceback(self, seq1, seq2, matrix):
        """
        Se empieza desde el máximo y termina en el primer 0
        """
        max_pos, _ = self.find_max(matrix)
        i, j = max_pos
        
        aligned_seq1 = []
        aligned_seq2 = []
        
        # Seguir hasta encontrar un 0
        while i > 0 and j > 0 and matrix[i][j] > 0:
            score_current = matrix[i][j]
            score_diagonal = matrix[i-1][j-1] + self.calculate_score(seq1[i-1], seq2[j-1])
            score_up = matrix[i-1][j] + self.gap
            score_left = matrix[i][j-1] + self.gap
            
            if score_current == score_diagonal:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif score_current == score_up:
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j-1])
                j -= 1
        
        return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


    def alignment(self, seq1, seq2):
        """
        Ejecuta el alineamiento local completo
        """
        self.validate_seqs(seq1, seq2)
        
        matrix = self.initialize_matrix(seq1, seq2)
        matrix = self.fill_matrix(matrix, seq1, seq2)
        aligned_seq1, aligned_seq2 = self.traceback(seq1, seq2, matrix)
        
        _, score = self.find_max(matrix)
        
        return aligned_seq1, aligned_seq2, score, matrix
