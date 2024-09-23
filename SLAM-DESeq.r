#!/usr/bin/env Rscript

# i was having some issues with this i think due to an updated DESEQ2 version, but I am hoping that it is now operational as of april 9, 2024

hisat_out <- "/hisat-3n-FR"

BED <- "WG"
BED_FILE <- "/gencode.v25/gencode.v25.annotation_genes.bed" # this is derived from the gtf file used in HTSeq

SLAM_summary_path <- "/SLAM_summary.txt" # this helps identify which files are replicates. See SLAM_summary.txt for example.

total_filter <- 20
nascent_filter <- 2 # 15min 4sU labeling
nascent_filter <- 20 # 120min 4sU labeling
nascent_filter <- 10 # 45min 4sU labeling

# Shouldn't need to modify anything below this point when running
library(DESeq2)
library(data.table)
library(dplyr)
library(stringr)
library(comprehenr)

summary_df <- read.delim(SLAM_summary_path, sep = "\t", header = TRUE)

COUNT_FILE_LIST <- unlist(list.dirs(hisat_out,full.names = FALSE, recursive = FALSE))
SUMMARY_LIST <- summary_df$Library
COUNT_FILES <- to_vec(for(FILE in COUNT_FILE_LIST) if (FILE %in% SUMMARY_LIST) FILE)

out_dir <- paste0(hisat_out,"/testDESEQ-",BED)
dir.create(out_dir)

#imports the information from the bed file used for calling counts by HTSeq
#this is used as the backbone for counts
BED_TABLE <- read.delim(BED_FILE, sep = "\t")
colnames(BED_TABLE) <- c("chr", "start", "end",
                        "gene_id", "idk", "Strand", "gene_name")
counts_all <- BED_TABLE[c("chr", "start", "end",
                        "gene_id", "Strand", "gene_name")]

counts_all["id"] <- paste0(counts_all$chr,":",counts_all$start,"-",counts_all$end,"_",counts_all$gene_id)

#loop to extract counts from the files
#also extracts the sum of total_read_counts column to use as library size
data_cols <- c()
treatments <- c()

coldata <- data.frame()
coldata_tc <- data.frame()

for (lib in COUNT_FILES){
    sample_treatment <- paste0(filter(summary_df, Library == lib)$Treatment[1],"_",filter(summary_df, Library == lib)$Treatment.length[1])
    replicate <- paste0("R",filter(summary_df, Library == lib)$Replicate[1])

    data_cols <- append(data_cols, paste0(lib, ":",sample_treatment,"-",replicate))
    treatments <- append(treatments,sample_treatment)

    # import HTseq count file
    file <- paste0(hisat_out, "/", lib, "/counts/HTSeq_counts_gene.tsv")
    table <- read.delim(file, header = TRUE, sep = "\t")

    # subset the table based on the important columns, including the nascent or total reads column
    counts <- table[c("chr", "start", "end", "gene_id", "TC_read_counts","total_read_counts")]
    colnames(counts) <- c("chr", "start", "end", "gene_id", paste0(lib, ":",sample_treatment,"-",replicate,"_tc"),paste0(lib, ":",sample_treatment,"-",replicate,"_total"))
    counts$start <- counts$start + 1

    # create a merged table containing the additional information from the bed file used for SLAMdunk
    # this will also iteratively add count information for each file
    counts_all <- merge(counts_all, counts,
                        by = c("chr", "start", "end", "gene_id"))

    # add sum of total_read_counts column to reads_df
    coldata[paste0(lib, ":",sample_treatment,"-",replicate,'_total'),"lib_sizes"] <- sum(table["total_read_counts"])
    coldata[paste0(lib, ":",sample_treatment,"-",replicate,'_total'),"condition"] <- sample_treatment
    coldata_tc[paste0(lib, ":",sample_treatment,"-",replicate,'_tc'),"lib_sizes"] <- sum(table["total_read_counts"])
    coldata_tc[paste0(lib, ":",sample_treatment,"-",replicate,'_tc'),"condition"] <- sample_treatment
}

treatments <- unique(treatments)

# for some reason rearranges the differential analysis order if you provide strings rather than numbers
rep <- 1
for (i in treatments){
    coldata_tc$condition <- str_replace(coldata_tc$condition,i,toString(rep))
    coldata$condition <- str_replace(coldata$condition,i,toString(rep))
    rep <- rep + 1
}
coldata_tc <- transform(coldata_tc, condition = as.numeric(condition))
coldata <- transform(coldata, condition = as.numeric(condition))
coldata_tc <- transform(coldata_tc, condition = LETTERS[condition])
coldata <- transform(coldata, condition = LETTERS[condition])

# creates a list of pairwise comparisons to make based on the different treatment conditions
combos <- combn(unique(treatments), 2)

write.table(counts_all,
    (paste0(out_dir, "/", BED, "_merged_counts.txt")),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE)

dds <- DESeqDataSetFromMatrix(
    countData = counts_all[c("id", paste(data_cols, "_tc", sep=""))],
    colData = coldata_tc[paste(data_cols, "_tc", sep=""), ],
    design =  ~ condition,
    tidy = TRUE
    )

dds.total <- DESeqDataSetFromMatrix(
    countData = counts_all[c("id", paste(data_cols, "_total", sep=""))],
    colData = coldata[paste(data_cols, "_total", sep=""), ],
    design =  ~ condition,
    tidy = TRUE
    )

featureData <- data.frame(basepairs=(counts_all$end - counts_all$start))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds.total) <- DataFrame(mcols(dds.total), featureData)

#  Run deseq main command.
dds.total <- DESeq(dds.total,)
sizeFactors(dds) <- sizeFactors(dds.total) # apply size factors to tc counts

size_factors <- data.frame(sizeFactors(dds))
colnames(size_factors) <- c('sizeFactor')
size_factors$reciprocal <- 1/size_factors$sizeFactor

counts_fpkm <- as.data.frame(fpkm(dds, robust = TRUE))
counts_fpkm$id <- rownames(counts_fpkm)
counts_fpkm <- merge(counts_all[c("id", "chr", "start", "end",
                    "gene_id", "Strand", "gene_name")],counts_fpkm, by="id")
counts_fpkm$id <- NULL

write.table(counts_fpkm,
    (paste0(out_dir, "/", BED, "_nascentRNA_FPKM.txt")),
    sep = "\t",
    col.names = colnames(counts_fpkm),
    row.names = FALSE,
    quote = FALSE)

counts_fpkm <- as.data.frame(fpkm(dds.total, robust = TRUE))
counts_fpkm$id <- rownames(counts_fpkm)
counts_fpkm <- merge(counts_all[c("id", "chr", "start", "end",
                    "gene_id", "Strand", "gene_name")],counts_fpkm, by="id")
counts_fpkm$id <- NULL

write.table(counts_fpkm,
    (paste0(out_dir, "/", BED, "_totalRNA_FPKM.txt")),
    sep = "\t",
    col.names = colnames(counts_fpkm),
    row.names = FALSE,
    quote = FALSE)

for (pairs in seq_len(ncol(combos))){
    # this determines which columns in the counts_all_filtered correspond
    # to the data of interest for differential analysis
    # this also ensures that they are in the correct order

    use_cols <- c(data_cols[grepl(paste0(":",combos[1, pairs],"-R"), data_cols)],data_cols[grepl(paste0(":",combos[2, pairs],"-R"), data_cols)])
    min.rep <- min(length(c(data_cols[grepl(combos[1, pairs], data_cols)])),length(c(data_cols[grepl(combos[2, pairs], data_cols)])))

    dds <- DESeqDataSetFromMatrix(
        countData = counts_all[c("id", paste(use_cols, "_tc", sep=""))],
        colData = coldata_tc[paste(use_cols, "_tc", sep=""), ],
        design =  ~ condition,
        tidy = TRUE
        )

    dds.total <- DESeqDataSetFromMatrix(
        countData = counts_all[c("id", paste(use_cols, "_total", sep=""))],
        colData = coldata[paste(use_cols, "_total", sep=""), ],
        design =  ~ condition,
        tidy = TRUE
        )

    featureData <- data.frame(basepairs = (counts_all$end - counts_all$start),gene=counts_all$gene_name)
    mcols(dds) <- DataFrame(mcols(dds), featureData)
    mcols(dds.total) <- DataFrame(mcols(dds.total), featureData)

    #  Run deseq main command.
    dds.total <- DESeq(dds.total)
    sizeFactors(dds) <- sizeFactors(dds.total) # apply size factors to tc counts
    dds <- DESeq(dds)

    # remove lowly expressed genes based on total rather than nascent
    keep <- rowSums(counts(dds.total) >= total_filter) >= min.rep
    dds.total <- dds.total[keep]
    keep <- rowSums(counts(dds) >= nascent_filter) >= min.rep
    dds <- dds[keep]


    for (i in c("nascent", "total")){
        # just for ease, I am going to separate out the different analyses
        # based on whether they include DMSO as one of the pairs
        if (grepl("DM",combos[1, pairs], fixed=TRUE) == TRUE || grepl("DM",combos[2, pairs], fixed=TRUE) == TRUE) {
            out_loc <- paste0(out_dir, "/DMSO_deseq-analysis")
        } else {
            out_loc <- paste0(out_dir, "/nonDMSO_deseq-analysis")
        }
        dir.create(out_loc)

        # get the normalized counts for GSEA
        if (i == "nascent") {norm_counts <- as.data.frame(counts(dds, normalized = TRUE))}
        if (i == "total") {norm_counts <- as.data.frame(counts(dds.total, normalized = TRUE))}
        norm_counts$id <- rownames(norm_counts)
        rownames(norm_counts) <- NULL
        norm_counts <- merge(counts_all[c("gene_name","id")],norm_counts, by="id",all.x=FALSE,all.y=FALSE)
        norm_counts$id <- "NA"

        if (i == "nascent") {norm_counts <- norm_counts[,c("gene_name","id",paste(use_cols, "_tc", sep=""))]}
        if (i == "total") {norm_counts <- norm_counts[,c("gene_name","id",paste(use_cols, "_total", sep=""))]}

        # rename the replicate columns as the same thing
        col_list <- c()
        for (treat in c(combos[1, pairs],combos[2, pairs])){
            treat_cols <- c(use_cols[grepl(treat, use_cols)])
            col_list <- append(col_list,rep(treat,length(treat_cols)))
        }

        setnames(norm_counts,old=c("gene_name","id"),new=c("NAME","DESCRIPTION"))

        dir.create(paste0(out_dir, "/GSEA"))
        dir.create(paste0(out_dir, "/GSEA/",i,"RNA_expression_GSEA"))
        out <- paste0(out_dir, "/GSEA/",i,"RNA_expression_GSEA/",combos[1, pairs], "_vs_", combos[2, pairs])
        dir.create(out)
        write.table(norm_counts,
            (paste0(out, "/", i, "RNA_norm.txt")),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

        # create a class file
        cls <- data.frame(matrix(ncol = length((col_list)), nrow = 0))
        cls[1,1:3] <- list(length(col_list), 2, 1)
        cls[2,"X1"] <- "#"
        cls[2,2:3] <- c(combos[1, pairs],combos[2, pairs])
        cls[3,1:length(col_list)] <- col_list

        write.table(cls,
            (paste0(out, "/GSEA.cls")),
            sep = "\t",
            na = "",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

        if (i == "nascent") {res <- results(dds)}
        if (i == "total") {res <- results(dds.total)}

        res$id <- rownames(res)
        results_df <- as.data.frame(res)

        if (i == "nascent") {col_list <- paste(use_cols, "_tc", sep = "")}
        if (i == "total") {col_list <- paste(use_cols, "_total", sep = "")}

        results_df <- merge(counts_all[c("id", "chr", "start", "end",
                            "gene_id", "Strand", "gene_name", col_list)],results_df, by="id",all=FALSE)
        results_df$id <- NULL

        # now we save the toptags output
        out <- paste0(combos[1, pairs], "_vs_", combos[2, pairs], "_", BED,"_",i)
        write.table(results_df[order(results_df$padj), ],
            (paste0(out_loc, "/", out, "RNA_differential_analysis-DEseq.txt")),
            sep = "\t",
            col.names = colnames(results_df),
            row.names = FALSE,
            quote = FALSE)

    }

}