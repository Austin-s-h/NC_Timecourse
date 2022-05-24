# Running FIMO for TFAP2A on the ESC vs NC Specific dual-motif locations
# First get FASTA for motif locations.
bedtools getfasta -fi /data/Austin/workdir/genome/Gallus_gallus.GRCg6a.dna.toplevel_withchr.fa -bed ESC_Specific_motifs.bed > ESC_Specific_motifs.fa
bedtools getfasta -fi /data/Austin/workdir/genome/Gallus_gallus.GRCg6a.dna.toplevel_withchr.fa -bed NC_Specific_motifs.bed > NC_Specific_motifs.fa

# Now FIMO.
fimo --thresh 1e-3 -o TFAP2A_ESC_Specific ../Imports/MA0003.4.meme ESC_Specific_motifs.fa
fimo --thresh 1e-3 -o TFAP2A_NC_Specific ../Imports/MA0003.4.meme NC_Specific_motifs.fa

fimo --thresh 1e-3 -o PAX7_ESC_Specific ../Imports/MA0780.1.meme ESC_Specific_motifs.fa
fimo --thresh 1e-3 -o PAX7_NC_Specific ../Imports/MA0780.1.meme NC_Specific_motifs.fa

fimo --thresh 1e-3 -o NANOG_ESC_Specific ../Imports/UN0383.1.meme ESC_Specific_motifs.fa
fimo --thresh 1e-3 -o NANOG_NC_Specific ../Imports/UN0383.1.meme NC_Specific_motifs.fa

fimo --thresh 1e-3 -o KLF4_ESC_Specific ../Imports/MA0039.4.meme ESC_Specific_motifs.fa
fimo --thresh 1e-3 -o KLF4_NC_Specific ../Imports/MA0039.4.meme NC_Specific_motifs.fa

fimo --thresh 1e-3 -o ZIC1_ESC_Specific ../Imports/MA0696.1.meme ESC_Specific_motifs.fa
fimo --thresh 1e-3 -o ZIC1_NC_Specific ../Imports/MA0696.1.meme NC_Specific_motifs.fa

fimo --thresh 1e-3 -o FOXD3_ESC_Specific ../Imports/MA0041.1.meme ESC_Specific_motifs.fa
fimo --thresh 1e-3 -o FOXD3_NC_Specific ../Imports/MA0041.1.meme NC_Specific_motifs.fa

