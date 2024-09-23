#!/usr/bin/env python

master_dir = '~/hisat-3n-FR' # the output directory from SLAM_hisat-3n
min_MQ = 10 # minimum MAPQ
gtf_path = '~/gencode.v25/gencode.v25.annotation.gtf' # path to gtf gile

import pandas as pd
import os
import collections
import HTSeq

def main(master_dir):
    file_list = os.listdir(master_dir)

    for i in file_list:

        bam_dir = master_dir + '/' + i
        name_sorted =  bam_dir + '/map/hisat2_aligned.nodup.name_sorted.bam'

        counts(name_sorted=name_sorted,bam_dir=bam_dir)

def counts(bam_dir,name_sorted):
    out_dir = bam_dir + '/counts'
    os.makedirs(out_dir,exist_ok=True)

    features = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    almnt_file = HTSeq.BAM_Reader(name_sorted)

    gtf_file = HTSeq.GFF_Reader(gtf_path)
    
    rep = 0
    for use_feature in ['gene']:
        feature_df = pd.DataFrame([])
        feature_rep = 0
        for feature in gtf_file:
            feature_rep += 1
            if feature_rep == 1000: break
            if feature.type == use_feature:
                # getting all of the information I want for each of the features
                features[feature.iv] += feature.attr["gene_id"]
                chr = str(feature.iv).split(':')[0]
                start = str(feature.iv).split(':')[1].split('/')[0].split(',')[0].replace('[','')
                end = str(feature.iv).split(':')[1].split('/')[0].split(',')[1].replace(')','')
                strand = str(feature.iv).split(':')[1].split('/')[1]
                data = pd.DataFrame([chr,start,end,strand,feature.attr["gene_id"],feature.attr["gene_name"]])
                feature_df = pd.concat([feature_df,data.transpose()])
        feature_df.columns = ['chr','start','end','strand',"gene_id","gene_name"]
        counts = collections.Counter()
        counts_TC = collections.Counter()
        converted_T = collections.Counter()
        unconverted_T = collections.Counter()
        total_T = collections.Counter()
        misc_counts = collections.Counter()  

        for bundle in HTSeq.pair_SAM_alignments(almnt_file, bundle=True ):
            if rep == 1000: break

            if len(bundle) != 1:
                continue  # Skip multiple alignments

            first_almnt, second_almnt = bundle[0]  # extract pair

            if (first_almnt is None) or (second_almnt is None):
                if (first_almnt is None): almnt = second_almnt
                if (second_almnt is None): almnt = first_almnt

                if (not almnt.aligned):
                    misc_counts["_unmapped"] += 1
                    continue
                # remove supplemental reads
                if first_almnt.supplementary or second_almnt.supplementary:
                    misc_counts["_supplementary"] += 1
                    continue

                try: Yf_almnt = almnt.optional_field('Yf')
                except: Yf_almnt = None

                # removes read pairs if one of them is missing the Yf field
                if Yf_almnt is None:
                    misc_counts["_noYf"] += 1
                    continue

                # NH tag reflects multiple alignments I think. Remove read pairs if one of the reads aligns multiple times.
                if (almnt.optional_field("NH") > 1) or (almnt.optional_field("NH") > 1):
                    misc_counts['multi_align'] += 1
                    continue

                # remove read pair if one of the reads doesn't meet the minimum read quality set by the user
                if (almnt.aQual < min_MQ) or (almnt.aQual < min_MQ):
                    misc_counts['low_mq'] += 1
                    continue

            else:
                almnt = None

                # removes read pair if both of the reads do not align
                # if only one read aligns, it switches to only using that read
                if (not first_almnt.aligned) and (not second_almnt.aligned):
                    misc_counts["_unmapped"] += 1
                    continue
                if (not first_almnt.aligned): almnt = second_almnt
                if (not second_almnt.aligned): almnt = first_almnt

                # removes read pair if both of the reads are supplementary
                # if only one read is supplementary, it switches to only using its read pair
                if almnt is None:
                    if (first_almnt.supplementary) and (second_almnt.supplementary):
                        misc_counts["_supplementary"] += 1
                        continue
                    if (first_almnt.supplementary): almnt = second_almnt
                    if (second_almnt.supplementary): almnt = first_almnt
                elif almnt.supplementary:
                    misc_counts["_supplementary"] += 1
                    continue
                
                # first have to check if the reads have the Yf field
                # removes read pair if neither of the reads contain the Yf field
                # if only one read has the Yf field, it switches to only using that read
                try: Yf_first = first_almnt.optional_field('Yf')
                except: Yf_first = None

                try: Yf_second = second_almnt.optional_field('Yf')
                except: Yf_second = None

                try: Yf_almnt = almnt.optional_field('Yf')
                except: Yf_almnt = None

                if almnt is None:
                    if (Yf_first is None) and (Yf_second is None):
                        misc_counts["_noYf"] += 1
                        continue
                    if (Yf_first is None): almnt = second_almnt
                    if (Yf_second is None): almnt = first_almnt
                elif (Yf_almnt is None):
                    misc_counts["_noYf"] += 1
                    continue

                # NH tag reflects multiple alignments I think. Remove read pair if both of the reads align multiple times.
                if almnt is None:
                    if (first_almnt.optional_field("NH") > 1) and (second_almnt.optional_field("NH") > 1):
                        misc_counts['multi_align'] += 1
                        continue
                elif (almnt.optional_field("NH") > 1):
                        misc_counts['multi_align'] += 1
                        continue         

                # removes read pair if both of the reads have a mapping quality below the user-defined value
                # if only one read exceeds the required mapping qualuty, it switches to only using that read
                if almnt is None:
                    if (first_almnt.aQual < min_MQ) and (second_almnt.aQual < min_MQ):
                        misc_counts['low_mq'] += 1
                        continue
                    if (first_almnt.aQual < min_MQ): almnt = second_almnt
                    if (second_almnt.aQual < min_MQ): almnt = first_almnt
                elif (almnt.aQual < min_MQ):
                    misc_counts['low_mq'] += 1
                    continue

            rep += 1
            gene_ids = set()
            if (almnt is None) or ((almnt is not None) and (almnt.pe_which == "first")):
                if almnt is None: use_almnt = first_almnt
                else: use_almnt = almnt
                for cigop in use_almnt.cigar:
                    if cigop.type != "M":
                        continue
                    for iv, val in features[cigop.ref_iv].steps():
                        gene_ids |= val

            if (almnt is None) or ((almnt is not None) and (almnt.pe_which == "second")):
                if almnt is None: use_almnt = second_almnt
                else: use_almnt = almnt
                for cigop in use_almnt.cigar:
                    if cigop.type != "M":
                        continue
                    for iv, val in features[invert_strand(cigop.ref_iv)].steps():
                        gene_ids |= val

                        # this is how the second read was dealt with in the count.py file. Not super sure of the application of itertools to it.
                        # the first strand was done how I did it above or at least I think
                        # iv_seq = itertools.chain(
                        #         iv_seq,
                        #         (invert_strand(co.ref_iv) for co in r[1].cigar
                        #         if co.type in com and co.size > 0))


            if len(gene_ids) == 1:
                gene_id = list(gene_ids)[0]
                counts[gene_id] += 1

                if almnt is None:
                    if (first_almnt.optional_field("Yf") > 0) or (second_almnt.optional_field("Yf") > 0):
                        counts_TC[gene_id] += 1
                    Yf_sum = first_almnt.optional_field("Yf") + second_almnt.optional_field("Yf")
                    Zf_sum = first_almnt.optional_field("Zf") + second_almnt.optional_field("Zf")
                else:
                    if (almnt.optional_field("Yf") > 0): counts_TC[gene_id] += 1
                    Yf_sum = almnt.optional_field("Yf")
                    Zf_sum = almnt.optional_field("Zf")

                converted_T[gene_id] += Yf_sum
                unconverted_T[gene_id] += Zf_sum
                total_T[gene_id] += Yf_sum + Zf_sum
            elif len(gene_ids) == 0:
                misc_counts["_no_feature"] += 1
            else:
                misc_counts["_ambiguous"] += 1

        data_df = pd.DataFrame([])
        for gene_id in counts:
            if total_T[gene_id] > 0: proportion = converted_T[gene_id]/total_T[gene_id]
            else: proportion = 0

            data = pd.DataFrame({'gene_id':gene_id,'TC_read_counts':counts_TC[gene_id],'total_read_counts':counts[gene_id],'total_converted_T':converted_T[gene_id],'total_unconverted_T':unconverted_T[gene_id],'total_T':total_T[gene_id],'conversion_rate':proportion},index=[0])
            data_df = pd.concat([data_df,data])
            
        output_df = feature_df.merge(data_df,on='gene_id',how='right')
        output_df.to_csv(out_dir + '/HTSeq_counts_' + use_feature + '.tsv',sep='\t',header=True,index=False)

        misc_df = pd.DataFrame([])
        for notation in misc_counts:
            data = pd.DataFrame({notation:misc_counts[notation]},index=[0])
            misc_df = pd.concat([misc_df,data])
        misc_df.to_csv(out_dir + '/HTSeq_excluded_counts_' + use_feature + '.tsv',sep='\t',header=True,index=False)

def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2

if __name__ == "__main__":
    main(master_dir)