import pandas as pd
import subprocess
import re  # For regular expressions

# Load the list of organisms from the CSV file with the new headers
file_path = "eukaryotic_organisms_with_genomes.csv"
df = pd.read_csv(file_path)

def check_genome_availability(organism_name):
    """Check if an organism has a reference genome available using ncbi-genome-download."""
    try:
        result = subprocess.run([
            "ncbi-genome-download", "--dry-run", "--formats", "genbank",
            "all", "-s", "refseq", "-g", organism_name
        ], capture_output=True, text=True)

        print(f"Results for {organism_name}:")
        print(result.stdout)

        # Look for any GCF reference genome ID in the output using regex
        genome_match = re.search(r'GCF_\d+\.\d+', result.stdout)
        if genome_match:
            genome_id = genome_match.group(0)  # Get the first matching genome ID
            return (True, genome_id)
        else:
            return (False, None)

    except Exception as e:
        print(f"Error checking {organism_name}: {e}")
        return (False, None)

# Apply the check and add the results as a new column without modifying the original columns
df["Genome Check Result"] = df["Organism Name"].apply(check_genome_availability)

# Save the updated DataFrame to a new CSV file
output_file = "eukaryotic_organisms_with_genomes_with_availability.csv"
df.to_csv(output_file, index=False)

print(f"Genome availability check complete. Results saved to {output_file}")
