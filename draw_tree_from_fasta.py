#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 23:39:52 2024

@author: nori
"""
import os
import subprocess
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from ete3 import Tree, TreeStyle, TextFace, NodeStyle
#import matplotlib.pyplot as plt
#from reportlab.lib.pagesizes import A3
#from reportlab.pdfgen import canvas
import re
import shutil


#### Step 1: Preparation ####
working_directory = "/Users/nori/Dropbox/BBSRC_EvoFPs/FP_phylogenetics/FPs_normal_chrmophore"
os.chdir(working_directory)

input_fasta_file = "normal_chromophore_all.fasta"

# offset of the Em_max wavelenth for adjusting the peak wavelenth to the apparent color
wl_offset = 15

#### Step 2: Check the FPs are already aligned and if not do alignment
## get the file name of the aligned FASTA file
# Split the filename into name and extension
name, ext = os.path.splitext(input_fasta_file)
# Insert "_aligned" and reassemble
aligned_fasta_file = f"{name}_aligned{ext}"

if os.path.exists(aligned_fasta_file):
    alignment = AlignIO.read(aligned_fasta_file, "fasta")
    print(f"Loaded a pre-existing alignment file {aligned_fasta_file}")
else:
    print(f"{aligned_fasta_file} not found. Performing MAFFT.")
    os.environ["PATH"] += os.pathsep + "/opt/homebrew/bin"
    mafft_path = shutil.which("mafft")
    bash_command = f"{mafft_path} --auto {input_fasta_file} > {aligned_fasta_file}"
    print(bash_command)
    subprocess.run(bash_command, shell=True)
    alignment = AlignIO.read(aligned_fasta_file, "fasta")
    print(f"Loaded the newly made alignment {aligned_fasta_file}")


#### Step 3: Compute distance matrix ####
#calculator = DistanceCalculator('identity')
calculator = DistanceCalculator('blosum62')
distance_matrix = calculator.get_distance(alignment)


#### Step 4: Construct phylogenetic tree using Neighbor Joining ####
## make a phylogenetic tree
constructor = DistanceTreeConstructor()
phylo_tree = constructor.nj(distance_matrix)

## save the tree in the Newick format 
# def quote_node_labels(tree):
#     # Regex to detect labels that need quoting (contains special characters)
#     special_chars = re.compile(r"[^a-zA-Z0-9_]")
#     for clade in tree.find_clades():
#         if clade.name and special_chars.search(clade.name):
#             clade.name = f'"{clade.name}"' # Add single quotes around the name


# quote_node_labels(phylo_tree)

output_tree_file = f"{name}_tree.nw"
Phylo.write(phylo_tree, output_tree_file, "newick")


#### Step 5: Define functions for customizing the tree visuals with ETE Toolkit 

## convert the wavelength of light to RGB
def wavelength_to_rgb(wavelength):
    if 380 <= wavelength < 440:
        R = -(wavelength - 440) / (440 - 380)
        G = 0.0
        B = 1.0
    elif 440 <= wavelength < 490:
        R = 0.0
        G = (wavelength - 440) / (490 - 440)
        B = 1.0
    elif 490 <= wavelength < 510:
        R = 0.0
        G = 1.0
        B = -(wavelength - 510) / (510 - 490)
    elif 510 <= wavelength < 580:
        R = (wavelength - 510) / (580 - 510)
        G = 1.0
        B = 0.0
    elif 580 <= wavelength < 645:
        R = 1.0
        G = -(wavelength - 645) / (645 - 580)
        B = 0.0
    elif 645 <= wavelength <= 780:
        R = 1.0
        G = 0.0
        B = 0.0
    else:
        R = G = B = 0.0  # Wavelength is outside the visible range

    if 380 <= wavelength < 420:
        factor = 0.3 + 0.7 * (wavelength - 380) / (420 - 380)
    elif 645 <= wavelength <= 780:
        factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 645)
    else:
        factor = 1.0

    R = int((R * factor) * 255)
    G = int((G * factor) * 255)
    B = int((B * factor) * 255)
    
    return (R, G, B)

## customize the tree visuals
def customize_tree(tree, font_size=6, line_width=1):
    for node in tree.traverse():
        # Customize node style to adjust line thickness and remove circles
        nstyle = NodeStyle()
        nstyle["size"] = 0  # Hide circles at nodes
        nstyle["vt_line_width"] = line_width  # Set vertical line width
        nstyle["hz_line_width"] = line_width  # Set horizontal line width
        node.set_style(nstyle)

        # Check the leaf name for Ex_max and Em_max and if successful, set the text color
        if node.is_leaf():
            match = re.search(r'\|\d+/\d+$', node.name)
            if match:
                em = int(match.group(0).split('/')[-1])  # Extract the [Em] value
                em += wl_offset
                color = wavelength_to_rgb(em)  # Get RGB color based on [Em]
                hex_color = '#%02x%02x%02x' % color  # Convert to hex format
            else:
                hex_color = "#000000" 
            node.add_face(TextFace(node.name, fsize=font_size, fgcolor=hex_color), column=0)


#### Step 6: Visualize the tree ####

## 6-1 Read the tree
tree = Tree(output_tree_file , format=1, quoted_node_names=True)


## 6-2 Set the root
# Here we use the node between LanFP1 (lancelet) and CpYGFP (crustaceans)
ancester = tree.get_common_ancestor("Bflo|LanFP1|500/510","Cpop|CpYGFP|508/518")
tree.set_outgroup(ancester)

## 6-3 Apply small font size to the tree leaves
customize_tree(tree, font_size=6, line_width=1)  # Adjust font size as desired

## 6-4 Set the TreeStyle and draw into a PDF file
ts = TreeStyle()
ts.show_leaf_name = False  # Disable default leaf names to use custom TextFaces
ts.scale = 600  # Adjust scale to fit the tree better on A3; higher values shrink the tree
output_tree_pdf = f"{name}_tree.pdf"
tree.render(output_tree_pdf, w=1000, units="px", tree_style=ts)

