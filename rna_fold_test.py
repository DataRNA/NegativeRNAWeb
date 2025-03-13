from RNA import fold
import sys
import os

def test_rnafold(fasta_file):
    """
    Test RNA secondary structure prediction using ViennaRNA
    
    Args:
        fasta_file (str): Path to FASTA file containing RNA sequences
    """
    try:
        # Check if file exists
        if not os.path.exists(fasta_file):
            print(f"Error: File {fasta_file} not found")
            return
        
        # Read sequences from FASTA
        with open(fasta_file, 'r') as f:
            sequences = []
            current_seq = ""
            current_name = ""
            
            for line in f:
                if line.startswith('>'):
                    if current_seq:
                        sequences.append((current_name, current_seq))
                    current_name = line.strip()[1:]
                    current_seq = ""
                else:
                    current_seq += line.strip()
            
            if current_seq:
                sequences.append((current_name, current_seq))
        
        # Predict for each sequence
        print("\nResults:")
        print("-" * 50)
        
        for name, seq in sequences:
            # Predict secondary structure using ViennaRNA
            structure, mfe = fold(seq)
            
            print(f"\nSequence: {name}")
            print(f"RNA: {seq}")
            print(f"Structure: {structure}")
            print(f"Minimum Free Energy: {mfe} kcal/mol")
            print("-" * 30)
                
    except Exception as e:
        print(f"Error occurred: {str(e)}")

if __name__ == "__main__":
    # Example usage
    if len(sys.argv) > 1:
        fasta_file = sys.argv[1]
    else:
        # Example FASTA file
        fasta_file = "example.fasta"
        # Create example FASTA if it doesn't exist
        if not os.path.exists(fasta_file):
            with open(fasta_file, "w") as f:
                f.write(">example_sequence\n")
                f.write("GGGCUAUUAGCUCAGUUGGUUAGAGCGCACCCCUGAUAAGGGUGAGGUCGCUGAUUCGAAUUCAGCAUAGCCCA\n")
    
    test_rnafold(fasta_file) 