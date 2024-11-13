#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 20:00:43 2024

@author: mlskat
"""

import numpy as np
import matplotlib.pyplot as plt
import requests
import csv
import json
import os


##### Step 1: Set parameters ##### 
## 1-1 Working directory
working_dir = "/Users/nori/Dropbox/BBSRC_EvoFPs/FP_phylogenetics/FPs_normal_chrmophore"
os.chdir(working_dir)

## 1-2 Filenames to be used 
# Table of protein names and the species names made by fetch_protein_species.py
protein_species_table = "protein_species_table.csv" 

# List of FP names to be excluded. ex. lanRFP (red with unnusual chromophore)
list_proteins_to_exclude = "exceptions.txt" 

# List of species names to be included
list_to_include = "species_normal_chromophore.txt"

# Downloaded data on FPs. if not found, a new one will be made
all_proteins_file = "all_proteins.json" 

## 1-3 Output file
cyan_fasta_output = "normal_chromophore_cyan.fasta" # FASTA file to save the sequences
green_fasta_output = "normal_chromophore_green.fasta" # FASTA file to save the sequences
red_fasta_output = "normal_chromophore_red.fasta" # FASTA file to save the sequences
all_fasta_output = "normal_chromophore_all.fasta"

## 1-4 Parameters
# Minimum length of the FP sequence (shorter sequences will be skipped)
min_seq_length = 200
max_seq_length = 250

#### Step 2: Check if the FP data has been downloaded and if not, download it. ####
# Check if we already have the file saved locally
if os.path.exists(all_proteins_file):
    # Load all proteins data from the local file
    with open(all_proteins_file, "r") as file:
        all_proteins = json.load(file)
    print("Loaded protein data from local file.")
else:
    # If the file doesn't exist, download the data from the FPbase API
    url = "https://www.fpbase.org/api/proteins/?format=json"
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        all_proteins = response.json()
        
        # Save the downloaded data to a file for future use
        with open(all_proteins_file, "w") as file:
            json.dump(all_proteins, file)
        print("Downloaded and saved protein data to local file.")
    else:
        print("Failed to download data from FPbase API.")
        all_proteins = []


##### Step 3: Load lists of species to include and proteins to exclude ####
proteins_to_exclude = set()
with open(list_proteins_to_exclude, "r") as exclude_file:
    for line in exclude_file:
        proteins_to_exclude.add(line.strip())  # Add each line to the set, stripping whitespace
      
        
def species_to_species_code(species_name):
    parts = species_name.split()
    if len(parts) == 2:
        genus, specific = species_name.split()
        if specific=='sp.':
            species_code = genus
        else:
            species_code = f"{genus[0].upper()}{specific[:3].lower()}"
    else:
        species_code=parts[0]
    return species_code
    
      
#species_to_include = set()
species_codes_to_include = set()
with open(list_to_include, "r") as include_file:
    for line in include_file:
        species_name = line.strip()
        species_code = species_to_species_code(species_name)
        species_codes_to_include.add(species_code)

#### Step 4: Load species name mapping from protein_species_table.csv and create 4-letter codes ####
species_codes = {}
with open(protein_species_table, "r") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        species_name = row["species"]
        species_code = species_to_species_code(species_name)
        species_codes[row["name"]] = species_code  # Map protein name to species code


#### Step 5: Parse data and write to FASTA ####
with open(all_fasta_output, "w") as all_fasta_file, \
     open(cyan_fasta_output, "w") as cyan_fasta_file, \
     open(green_fasta_output, "w") as green_fasta_file, \
     open(red_fasta_output, "w") as red_fasta_file:
    for protein in all_proteins:
        name = protein.get("name")
        print(name)
        # Skip if the protein is in the exclusion list
        if name in proteins_to_exclude:
            continue

        # Get species code from species_codes dictionary and skip if it is not in the list to include
        species_code = species_codes.get(name, "Unknown")
        if species_code not in species_codes_to_include:
            continue

        # Extract the sequence, skip if there is no sequence
        sequence = protein.get("seq")
        if not sequence:
            continue
        
        # Skip if the sequence is too short
        if len(sequence) < min_seq_length or len(sequence) > max_seq_length :
            continue 

        # Check if protein is a "basic" switch type and get ex_max and em_max if available
        ex_max, em_max = "-", "-"
        switch_type=protein.get("switch_type")
        if switch_type == "b":
            if "states" in protein and protein["states"]:
                state = protein["states"][0]  # Assume the first state is the "default"
                ex_max = state.get("ex_max", "-")
                em_max = state.get("em_max", "-")
        if em_max=="-":
            continue
        if em_max is None:
            continue
        
        # Construct the FASTA header
        #safer_name = name.replace(".","")
        #safer_name = safer_name.replace(":","")
        #safer_name = safer_name.replace("|","")
       # safer_name = safer_name.replace(" ","_")
        header = f">'{species_code}|{name}|{ex_max}/{em_max}'"
        
        # Write to a FASTA file corresponding to the max emission wavelength
        all_fasta_file.write(f"{header}\n{sequence}\n")
        if em_max <485:
            cyan_fasta_file.write(f"{header}\n{sequence}\n")
        elif 485 <= em_max < 550:
            green_fasta_file.write(f"{header}\n{sequence}\n")
        elif em_max >= 550:
            red_fasta_file.write(f"{header}\n{sequence}\n")
    
    #cmFP512
    header=">'Cmem|cmFP512|503/512' tr|Q5ZQQ5|Q5ZQQ5_CERMM Green fluorescent protein FP512 OS=Cerianthus membranaceus OX=208460 PE=1 SV=1"
    sequence="MSQLDNNLSVSVYMKGNVNNHEFEYDGIGGGDPNSGQFSLKTKLRGGKPLPFSYDIITMGFQYGFRAFTKYPEGIADYFKGSFPEAFQWNRRIEFEDGGVINMSSDITYKDKVLHGDVWALGVNFPPNGPVMKNEIVMEEPAEETLTAKNGVLVGFCPKAYLLKDGSYYYGHMTTFYRSKKSGQPLPGFHFIKHRLVKTKVEPGFKMVEQAEYATAHVCDLPKPN"
    green_fasta_file.write(f"{header}\n{sequence}\n")

    #ccalRFP2
    header=">'Ccal|ccalRFP2|558/516' tr|Q1ALD5|Q1ALD5_CORYC GFP-type red fluorescent protein OS=Corynactis californica OX=44298 PE=2 SV=1"
    sequence="MSLSKQVLPHDVRMRYHMDGCVNGHSFTIEGEGAGKPYEGKKTLKLRVTKGGPLPFAFDILSATFTYGNRCFCEYPEDMPDYFKQSLPEGYSWERTMMYEDGGCGTSSAHIRLEKNCFVHQSTFLGVNFPANGPVMQKKALNWEPSSELITPCDGILKGDVTMFLMLEGGHRQKCQFTTSYKASKAVKMPPNHIIEHVVWKGEDSDGFQIKEHAVAKHFTVDVKET"
    red_fasta_file.write(f"{header}\n{sequence}\n")
    
    #GFPxm
    header=">'Amac|GFPxm/476/496' tr|Q8WP95|Q8WP95_9CNID Green fluorescent protein OS=Aequorea macrodactyla OX=147615 GN=GFPxm PE=2 SV=1"
    sequence="MSKGEELFTGIVPVLIELDGDVHGHKFSVRGEGEGDADYGKLEIKFICTTGKLPVPWPTLVTTFSYGIQCFARYPEHMKMNDFFKSAMPEGYIQERTIFFQDDGKYKTRGEVKFEGDTLVNRIELKGMDFKEDGNILGHKLEYNFNSHNVYIMPDKANNGLKVNFKIRHNIEGGGVQLADHYQTNVPLGDGPVLIPINHYLSTQTAISKDRNETRDHMVFLEFFSACGHTHGMDELYK"
    green_fasta_file.write(f"{header}\n{sequence}\n")    
    
print(f"Sequences saved to '{cyan_fasta_output}', '{green_fasta_output}' , and '{red_fasta_output}'")
