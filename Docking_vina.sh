#!/bin/bash

# Create required directories
mkdir -p Ligands/minimized
mkdir -p Ligands/pdbqt
mkdir -p Results

# Minimize all .sdf files in Ligands/
cd Ligands/
for file in *.sdf; do
  base=$(basename "$file" .sdf)
  min_file="minimized/min_${base}.sdf"
  obabel "$file" -O "$min_file" --minimize --ff UFF --cg --steps 200 --update 1 --criteria 0.1
  echo "Minimized: $file → $min_file"
done

# Convert minimized SDFs to PDBQT
for min_file in minimized/*.sdf; do
  base=$(basename "$min_file" .sdf)
  pdbqt_file="pdbqt/${base}.pdbqt"
  obabel "$min_file" -O "$pdbqt_file" --gen3d --partialcharge gasteiger
  echo "Converted: $min_file → $pdbqt_file"
done

cd pdbqt/


# Run AutoDock Vina
for f in *.pdbqt; do
  echo "Processing $f"
  /Users/students-ad/Desktop/After_MD_docking/Vina/vina --config ../../conf.txt --ligand "$f" --out ../../Results/"$f" --log ../../Results/"$f.txt"
done

cd ../../Results/

# Rename output files
for f in *.pdbqt; do
  mv "$f" "${f%.pdbqt}.pdbqt"
done


# Output file
output="../results.csv"
> "$output" # Clear the file if it exists
echo "filename,affinity" >> "$output"

# Loop over all .pdbqt.txt files
for file in *.pdbqt.txt; do
    if [[ -f "$file" ]]; then
        # Extract only the first line starting with mode 1 (first field is 1)
        affinity=$(awk '$1 == 1 {print $2; exit}' "$file")
        
        # Only add if an affinity was found
        if [[ -n "$affinity" ]]; then
            echo "$file,$affinity" >> "$output"
        fi
    fi
done

echo "Extraction complete. Results saved to $output"
