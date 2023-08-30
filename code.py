def align_pairwise_sequences(seq1, seq2, gap_penalty=-2, match_score=1, mismatch_penalty=-1):
    # Initialize the scoring matrix and traceback matrix
    len_seq1 = len(seq1) + 1
    len_seq2 = len(seq2) + 1
    dp_matrix = [[0] * len_seq2 for _ in range(len_seq1)]

    # Initialize the first row and column of the matrix
    for i in range(len_seq1):
        dp_matrix[i][0] = gap_penalty * i
    for j in range(len_seq2):
        dp_matrix[0][j] = gap_penalty * j

    # Fill in the matrix using dynamic programming
    for i in range(1, len_seq1):
        for j in range(1, len_seq2):
            match = dp_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = dp_matrix[i - 1][j] + gap_penalty
            insert = dp_matrix[i][j - 1] + gap_penalty
            dp_matrix[i][j] = max(match, delete, insert)

    # Traceback to find aligned sequences
    aligned_seq1, aligned_seq2 = "", ""
    i, j = len_seq1 - 1, len_seq2 - 1

    while i > 0 and j > 0:
        if dp_matrix[i][j] == dp_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif dp_matrix[i][j] == dp_matrix[i - 1][j] + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    # Complete traceback for remaining characters
    while i > 0:
        aligned_seq1 = seq1[i - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        i -= 1
    while j > 0:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[j - 1] + aligned_seq2
        j -= 1

    return aligned_seq1, aligned_seq2

def progressive_multiple_alignment(sequences, gap_penalty=-2, match_score=1, mismatch_penalty=-1):
    num_sequences = len(sequences)
   
    pairwise_aligned_sequences = [sequences[0]]
    for i in range(1, num_sequences):
        aligned_seq, _ = align_pairwise_sequences(pairwise_aligned_sequences[i - 1], sequences[i], gap_penalty, match_score, mismatch_penalty)
        pairwise_aligned_sequences.append(aligned_seq)
   
    profile_alignment = pairwise_aligned_sequences[0]
    for i in range(1, num_sequences):
        profile_alignment, _ = align_pairwise_sequences(profile_alignment, pairwise_aligned_sequences[i], gap_penalty, match_score, mismatch_penalty)
   
    return profile_alignment

# Collect user input for protein sequences
num_sequences = int(input("Enter the number of protein sequences: "))
sequences = []

for i in range(num_sequences):
    sequence = input(f"Enter protein sequence {i + 1}: ")
    sequences.append(sequence)

# Perform progressive multiple sequence alignment
aligned_sequence = progressive_multiple_alignment(sequences)

# Print the aligned sequence
print("Progressive Multiple Sequence Alignment:")
print(aligned_sequence)
