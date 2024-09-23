Filter vcf file:
	bcftools filter -O z -o RPE-1.variants.filtered.vcf.gz -i 'QUAL > 100 && INFO/DP > 10 && TYPE="SNP"' RPE-1.variants.vcf.gz

vcf file was broken down into the different chromosomes (i.e. chr1.vcf, chr2.vcf, etc.) using: split_vcf_by_chr.pl

filtered SNPs on non-decoy chromosomes in hs38d5.fa were replaced with N using fasta_N_sub.py
    dependencies are pandas, math, subprocess, pathlib, and os
number of lines in hs38d5_Nsub.fa is 51471479

to extract only decoy:

    chr sed -e '1,51471479d' hs38d5.fa > decoy_chr.fa

Combine the output of fasta_N_sub.py with the non-substituted decoy chromosomes:

    cat hs38d5_Nsub.fa decoy_chr.fa > ~/RPE1-SNPs/hs38d5_Nsub.with-decoy.fa

    vcf number of SNPs is 3832729
    vcf number of SNPs excluding decoy chr is 3799777
    hs38d5.fa 'N' count is 165046383
    expected 'N' count for hs38d5_Nsub.with-decoy.fa is therefore 165046383 + 3799777 = 168846160
    hs38d5_Nsub.with-decoy.fa 'N' count is 168846160

Build the T->C conversion genome needed for hisat-3n (see https://daehwankimlab.github.io/hisat2/hisat-3n/):
    hisat-3n-build --base-change T,C ~/RPE1-SNPs/hs38d5_RPE1_Nsub.with-decoy.fa ~/RPE1-SNPs/HISAT-3N_SLAMseq/hs38d5_RPE1_Nsub.with-decoy.TC

Then can run SLAM_hisat-3n.sh to perform the alignment.
    You will need to modify various things within this script before running.

Then can run SLAM_HTSeq_FR.py
    You will need to modify various things within this script before running.
    Dependencies are pandas, os, collections, and HTSeq

Then can run SLAM-DESeq.r
    You will need to modify various things within this script before running.
    You will also need to create a SLAM_summary.txt type file (see example).
    Dependencies are DESeq2, dplyr, data.table, stringr, and comprehenr
