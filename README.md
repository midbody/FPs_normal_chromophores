# FPs_normal_chromophores
This is to make FASTA sequence files for 'all', 'cyan', 'green', 'red' FPs of GFP-like FPs available at FPbase https://www.fpbase.org/
Since FPbase covers other luminescent proteins such as phytochrome and engineered constructs (tandem etc.), we need to select appropriate sequences for proper alignment. The criteria are
1. Sequence length is > 200 and < 250
2. The species name is found in "species_normal_chromophore.txt". 
This is to exclude bacteria- or plant-derived FPs not homologous to avGFP. Some animal species were not included. For example, KO (https://www.fpbase.org/protein/ko/) from Verrillofungia concinna (a stony coral) forms a chromophore with 4 residues. This species was excluded.
3. The names of the FP to be excluded are listed in "exceptions.txt". 
This is to exclude protein-by-protein. laRFP has unnusual extension of pi-system by conjugation with Y at n+2 position. zFP538, ZsYellow1, mPapaya, and mPapaya07 has KYG triplet, in which Lys forms an additional ring structure. cEGFP and cp-mKate are circularly permutated.

'normal_chromophore_xxx.fasta' files contain sequences before alignment

a. max_em < 485: 'cyan'

b. 485 <= max_em < 550: 'green'

c. max_em >= 550: 'red'

'normal_chromophore_all_aligned.fasta' contain sequences after alignment with MAFFT

'normal_chromophore_all_tree.nw' is the phylogenetic tree of the selected FP sequences (neighbour-joining with BLOSUM62)

'normal_chromophore_all_tree.pdf' is a graphical representation of the above tree. Out-group was set at the branching point between the lencets and crustaceans.

The headers in the FASTA files and the texts on tree leaves are in the format [Species code]|[FPbase entry name]|[peak Ex]/[peak Em]

The species code is "Hsap" for "Homo sapiens" and "Discosoma" for "Discosoma sp."
