#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=giulia.pontali@fht.org
#SBATCH --job-name bonferroni_significant
#SBATCH --output bonferroni_significant.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 16G
#SBATCH --time 8:00:00

source ~/.bashrc
module -s load singularity/3.8.5

# Set the source directory
source_dir="/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/qced_sumstats_digits/output"

# Set the destination directory
destination_dir="/scratch/giulia.pontali/pQTL_bonferroni_significant"

# Iterate over subdirectories starting with "seq"
for subdir in "$source_dir"/seq*; do
    # Check if the subdirectory exists
    if [ -d "$subdir" ]; then
        # Get the sequence number from the subdirectory name
        seq_number=$(basename "$subdir" | sed 's/^seq//')

        # Create a combined filename for the destination file
        dest_filename="${destination_dir}/seq${seq_number}.txt"

        # Find files ending with ".gwaslab.tsv.gz" and process them
        for file in "$subdir"/*.gwaslab.tsv.gz; do
            # Check if files exist
            if [ -e "$file" ]; then
                # Extract rows with MLOG10P > (-log10(5e-8/3978))
                zcat "$file" | awk -v threshold=$(awk 'BEGIN{print 10.90069}') 'BEGIN{OFS="\t"} $13 > threshold {print}' > temp_file.txt

                # Check if temp_file_with_seqid.txt has more than one line
                if [ $(wc -l < temp_file.txt) -gt 1 ]; then
                    # Move temp_file_with_seqid.txt to dest_filename
                    mv temp_file.txt "$dest_filename"
                else
                    # Remove temp_file_with_seqid.txt if it has only one line
                    rm temp_file.txt
                fi
            fi
        done
    fi
done
