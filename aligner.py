#!/usr/bin/env python3
"""
Pairwise Sequence Aligner
Implements the Needleman-Wunsch global alignment algorithm.
"""

import argparse
import sys
import os
from typing import Tuple, List, Dict

def needleman_wunsch(
    seq1: str,
    seq2: str,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2
) -> Tuple[int, str, str, List[List[int]]]:
    """
    Computes the global alignment of seq1 and seq2 using the Needleman-Wunsch algorithm.
    
    Args:
        seq1: First sequence string
        seq2: Second sequence string
        match: Score for matching characters
        mismatch: Penalty for mismatching characters
        gap: Penalty for inserting a gap
        
    Returns:
        A tuple of (alignment_score, aligned_seq1, aligned_seq2, score_matrix)
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize the scoring matrix with zeros
    # Dimensions: (m + 1) rows x (n + 1) columns
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize the traceback matrix
    # Stores directions: 'D' (Diagonal), 'U' (Up), 'L' (Left), 'S' (Stop/Start)
    traceback_matrix = [['S'] * (n + 1) for _ in range(m + 1)]
    
    # Fill in the base cases for the first row and first column
    for i in range(1, m + 1):
        score_matrix[i][0] = i * gap
        traceback_matrix[i][0] = 'U'  # Can only come from above (gap in seq2)
        
    for j in range(1, n + 1):
        score_matrix[0][j] = j * gap
        traceback_matrix[0][j] = 'L'  # Can only come from left (gap in seq1)
        
    # Fill in the rest of the dynamic programming matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate scores from the three possible directions
            char1 = seq1[i - 1]
            char2 = seq2[j - 1]
            
            diag_score = score_matrix[i - 1][j - 1] + (match if char1 == char2 else mismatch)
            up_score = score_matrix[i - 1][j] + gap
            left_score = score_matrix[i][j - 1] + gap
            
            # Select the maximum score (with standard tie-breaking: Diagonal > Up > Left)
            max_score = max(diag_score, up_score, left_score)
            score_matrix[i][j] = max_score
            
            # Record the traceback direction
            if max_score == diag_score:
                traceback_matrix[i][j] = 'D'
            elif max_score == up_score:
                traceback_matrix[i][j] = 'U'
            else:
                traceback_matrix[i][j] = 'L'
                
    # Reconstruct the alignment by tracing back from bottom-right (m, n) to top-left (0, 0)
    aligned1 = []
    aligned2 = []
    i, j = m, n
    
    while i > 0 or j > 0:
        direction = traceback_matrix[i][j]
        
        if direction == 'D':
            # Diagonal: align characters from both sequences
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == 'U':
            # Up: gap in seq2, character from seq1
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            i -= 1
        elif direction == 'L':
            # Left: gap in seq1, character from seq2
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1
        else:
            # We reached the top-left cell
            break
            
    # The traceback builds the alignment in reverse order, so we must flip them
    aligned_seq1 = "".join(reversed(aligned1))
    aligned_seq2 = "".join(reversed(aligned2))
    
    return score_matrix[m][n], aligned_seq1, aligned_seq2, score_matrix

def smith_waterman(
    seq1: str,
    seq2: str,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -2
) -> Tuple[int, str, str, List[List[int]]]:
    """
    Computes the local alignment of seq1 and seq2 using the Smith-Waterman algorithm.
    
    Args:
        seq1: First sequence string
        seq2: Second sequence string
        match: Score for matching characters
        mismatch: Penalty for mismatching characters
        gap: Penalty for inserting a gap
        
    Returns:
        A tuple of (alignment_score, aligned_seq1, aligned_seq2, score_matrix)
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize score and traceback matrices with zeros
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    traceback_matrix = [['S'] * (n + 1) for _ in range(m + 1)]
    
    max_score = 0
    max_pos = (0, 0)
    
    # Fill in the dynamic programming matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            char1 = seq1[i - 1]
            char2 = seq2[j - 1]
            
            diag_score = score_matrix[i - 1][j - 1] + (match if char1 == char2 else mismatch)
            up_score = score_matrix[i - 1][j] + gap
            left_score = score_matrix[i][j - 1] + gap
            
            # Local alignment score cannot be negative (floored at 0)
            cell_score = max(0, diag_score, up_score, left_score)
            score_matrix[i][j] = cell_score
            
            if cell_score == 0:
                traceback_matrix[i][j] = 'S'
            elif cell_score == diag_score:
                traceback_matrix[i][j] = 'D'
            elif cell_score == up_score:
                traceback_matrix[i][j] = 'U'
            else:
                traceback_matrix[i][j] = 'L'
                
            # Track the maximum score position in the matrix
            if cell_score > max_score:
                max_score = cell_score
                max_pos = (i, j)
                
    # Reconstruct the alignment by tracing back from max_pos
    aligned1 = []
    aligned2 = []
    i, j = max_pos
    
    while i > 0 and j > 0:
        if score_matrix[i][j] == 0:
            break
            
        direction = traceback_matrix[i][j]
        if direction == 'S':
            break
        elif direction == 'D':
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == 'U':
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            i -= 1
        elif direction == 'L':
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1
            
    aligned_seq1 = "".join(reversed(aligned1))
    aligned_seq2 = "".join(reversed(aligned2))
    
    return max_score, aligned_seq1, aligned_seq2, score_matrix

def print_matrix(seq1: str, seq2: str, score_matrix: List[List[int]]):
    """Prints the scoring matrix in a readable grid format."""
    # Header row with seq2 characters
    header = [" ", "-"] + [char for char in seq2]
    print("\nScoring Matrix:")
    print("      " + "   ".join(f"{h:>4}" for h in header))
    
    # First row (representing the initial gap state)
    row_str = f"   -  " + "   ".join(f"{val:>4}" for val in score_matrix[0])
    print(row_str)
    
    # Subsequent rows with seq1 characters
    for i in range(1, len(score_matrix)):
        row_char = seq1[i - 1]
        row_str = f"   {row_char}  " + "   ".join(f"{val:>4}" for val in score_matrix[i])
        print(row_str)
    print()

def export_matrix_to_csv(seq1: str, seq2: str, score_matrix: List[List[int]], file_path: str):
    """Exports the scoring matrix to a CSV file."""
    try:
        with open(file_path, 'w') as f:
            # Write header row
            f.write("," + ",".join(["-"] + [char for char in seq2]) + "\n")
            
            # Write first row
            f.write("-," + ",".join(str(val) for val in score_matrix[0]) + "\n")
            
            # Write subsequent rows
            for i in range(1, len(score_matrix)):
                row_char = seq1[i - 1]
                f.write(f"{row_char}," + ",".join(str(val) for val in score_matrix[i]) + "\n")
        print(f"Scoring matrix exported to: {file_path}")
    except IOError as e:
        print(f"Error exporting matrix to CSV: {e}", file=sys.stderr)

def read_fasta(file_path: str) -> Dict[str, str]:
    """Reads a FASTA file and returns a dictionary of headers to sequences."""
    sequences = {}
    current_header = None
    current_seq = []
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = "".join(current_seq)
            
    return sequences

def format_alignment(aligned1: str, aligned2: str) -> str:
    """Formats the alignment with match/mismatch indicators between sequences."""
    match_line = []
    for char1, char2 in zip(aligned1, aligned2):
        if char1 == '-' or char2 == '-':
            match_line.append(' ')
        elif char1 == char2:
            match_line.append('|')
        else:
            match_line.append('.')
            
    matches = "".join(match_line)
    
    # Calculate identity percentage
    identities = matches.count('|')
    total_len = len(aligned1)
    identity_pct = (identities / total_len * 100) if total_len > 0 else 0
    
    output = []
    output.append(f"Seq1: {aligned1}")
    output.append(f"      {matches}")
    output.append(f"Seq2: {aligned2}")
    output.append(f"\nIdentity: {identities}/{total_len} ({identity_pct:.1f}%)")
    return "\n".join(output)

def main():
    parser = argparse.ArgumentParser(
        description="Global pairwise sequence aligner using Needleman-Wunsch dynamic programming."
    )
    
    # Sequence Inputs (direct strings or files)
    group1 = parser.add_mutually_exclusive_group(required=False)
    group1.add_argument("-s1", "--seq1", type=str, help="First sequence string")
    group1.add_argument("-f1", "--file1", type=str, help="Path to first FASTA file")
    
    group2 = parser.add_mutually_exclusive_group(required=False)
    group2.add_argument("-s2", "--seq2", type=str, help="Second sequence string")
    group2.add_argument("-f2", "--file2", type=str, help="Path to second FASTA file")
    
    # Scoring Parameters
    parser.add_argument("-m", "--match", type=int, default=1, help="Match score (default: 1)")
    parser.add_argument("-ms", "--mismatch", type=int, default=-1, help="Mismatch penalty (default: -1)")
    parser.add_argument("-g", "--gap", type=int, default=-2, help="Gap penalty (default: -2)")
    
    # Output Options
    parser.add_argument("--mode", choices=["global", "local"], default="global", help="Alignment mode: 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman) (default: global)")
    parser.add_argument("--show-matrix", action="store_true", help="Print the scoring matrix")
    parser.add_argument("--export-csv", type=str, help="Export scoring matrix to CSV file")
    
    args = parser.parse_args()
    
    # Resolve sequence 1
    seq1 = None
    if args.file1:
        try:
            seqs = read_fasta(args.file1)
            if not seqs:
                print(f"Error: No sequences found in {args.file1}", file=sys.stderr)
                sys.exit(1)
            # Use the first sequence in the file
            header = list(seqs.keys())[0]
            seq1 = seqs[header]
            print(f"Loaded Sequence 1 from FASTA: {header} ({len(seq1)} bp/aa)")
        except Exception as e:
            print(f"Error reading sequence 1: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.seq1:
        seq1 = args.seq1.upper()
    
    # Resolve sequence 2
    seq2 = None
    if args.file2:
        try:
            seqs = read_fasta(args.file2)
            if not seqs:
                print(f"Error: No sequences found in {args.file2}", file=sys.stderr)
                sys.exit(1)
            header = list(seqs.keys())[0]
            seq2 = seqs[header]
            print(f"Loaded Sequence 2 from FASTA: {header} ({len(seq2)} bp/aa)")
        except Exception as e:
            print(f"Error reading sequence 2: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.seq2:
        seq2 = args.seq2.upper()
        
    # If parameters not provided, use default sequences for demo
    if not seq1 or not seq2:
        print("No input sequences specified. Running with default demo sequences HELLO and CELLO.")
        seq1 = "HELLO"
        seq2 = "CELLO"
        
    # Run the alignment
    if args.mode == "local":
        score, aligned1, aligned2, matrix = smith_waterman(
            seq1, seq2, match=args.match, mismatch=args.mismatch, gap=args.gap
        )
        mode_label = "LOCAL"
    else:
        score, aligned1, aligned2, matrix = needleman_wunsch(
            seq1, seq2, match=args.match, mismatch=args.mismatch, gap=args.gap
        )
        mode_label = "GLOBAL"
    
    print("\n" + "=" * 40)
    print(f"{mode_label} SEQUENCE ALIGNMENT REPORT")
    print("=" * 40)
    print(f"Sequence 1 Length: {len(seq1)}")
    print(f"Sequence 2 Length: {len(seq2)}")
    print(f"Parameters: Match={args.match}, Mismatch={args.mismatch}, Gap={args.gap}")
    print(f"Alignment Score: {score}")
    print("-" * 40)
    print(format_alignment(aligned1, aligned2))
    print("=" * 40)
    
    # Optional outputs
    if args.show_matrix:
        print_matrix(seq1, seq2, matrix)
        
    if args.export_csv:
        export_matrix_to_csv(seq1, seq2, matrix, args.export_csv)

if __name__ == "__main__":
    main()
