import mxfold2
import sys
import os

def test_mxfold2(fasta_file):
    """
    Test MXfold2 RNA secondary structure prediction
    
    Args:
        fasta_file (str): Path to FASTA file containing RNA sequences
    """
    try:
        # Check if file exists
        if not os.path.exists(fasta_file):
            print(f"Error: File {fasta_file} not found")
            return
        
        # Create MXfold2 model
        model = mxfold2.MXfold2()
        
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
            # Predict secondary structure
            structure, score = model.predict_structure(seq)
            
            print(f"\nSequence: {name}")
            print(f"RNA: {seq}")
            print(f"Structure: {structure}")
            print(f"Score: {score}")
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
    
    test_mxfold2(fasta_file) 