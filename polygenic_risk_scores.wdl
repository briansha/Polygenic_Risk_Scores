version 1.0

## Version 07-02-2021
##
## This workflow calculates polygenic risk scores for different summary statistics files.
##
## Cromwell version support - Successfully tested on v65
##
## Distributed under terms of the MIT License
## Copyright (c) 2021 Brian Sharber
## Contact <brian.sharber@vumc.org>

workflow polygenic_risk_scores {

    input {
        String plink2_docker = "briansha/plink2:polygenic"
        File sum_stats_files
    }
	Array[Array[String]] sum_stats_file_lines = read_tsv(sum_stats_files)

    scatter(items in sum_stats_file_lines) {
        call Calculate_Scores {
          input:
            docker = plink2_docker,
            formatted_sum_stats_file_name = items[0] + ".txt",
            concat_file_name = items[0] + "_concat.txt",
            output_prefix = items[0],
            mismatch_file_prefix = items[0],
            scored_file_prefix = items[0],
            sum_stats_file = items[1]
        }
    }

    output {
        Array[File] output_formatted_sum_stats_file = Calculate_Scores.output_formatted_sum_stats_file
        Array[File] output_concat_file = Calculate_Scores.output_concat_file
        Array[File] output_mismatch_file = Calculate_Scores.output_mismatch_file
        Array[File] output_scored_file = Calculate_Scores.output_scored_file
    }

    meta {
    	author : "Brian Sharber"
        email : "brian.sharber@vumc.org"
        description : "This workflow calculates polygenic risk scores for different summary statistics files."
    }
}

task Calculate_Scores {
    input {
        File vcf_col_split  # Python file for splitting columns in a vcf.
        File sum_stats_file # .vcf.gz file
        File bim_file_23
        File bed
        File bim
        File fam
        Int nGwas
        String formatted_sum_stats_file_name    # name for the generated formatted sum stats file.
        String mismatch_file_prefix             # prefix for the mismatched file.
        String scored_file_prefix               # prefix for the scored file.
        String concat_file_name                 # name for the generated concatenated file.
        String output_prefix                    # prefix for the output files for the --output_dir parameter in PRScs.
        File ref_dir_zipped                     # zipped directory containing the reference panels.

        String docker # Docker image containing PLINK2.
        Int disk = 200
        Float memory = 16.0
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0
    }
    String ref_dir = basename(ref_dir_zipped, ".zip")
    String vcf_prefix = basename(sum_stats_file, ".vcf.gz")
    String vcf_generated_es_file = vcf_prefix + "_ES.txt"
    String vcf_generated_lp_file = vcf_prefix + "_LP.txt"
    String vcf_generated_meta_file = vcf_prefix + "_Meta.txt"
    String sum_stats_file_in_current_directory = basename(sum_stats_file)

    # Move sum_stats_file to the current directory and use sum_stats_file_in_current_directory to indicate it.
    # - else, the ES, LP, and Meta files will have gs://... tacked onto them and not be in the current directory.
    command <<<
        unzip ~{ref_dir_zipped}
        git clone https://github.com/getian107/PRScs.git
        mv ~{sum_stats_file} .
        chmod 700 ~{vcf_col_split}
        python3 ~{vcf_col_split} ~{sum_stats_file_in_current_directory} ES LP

        python3 - ~{formatted_sum_stats_file_name} ~{vcf_generated_es_file} ~{vcf_generated_lp_file} ~{vcf_generated_meta_file} << 'END_SCRIPT'
        import pandas as pd
        import sys, os
        es = pd.read_csv(sys.argv[2], header=0, sep='\t', quotechar='"', error_bad_lines=False)
        lp = pd.read_csv(sys.argv[3], header=0, sep='\t', quotechar='"', error_bad_lines=False)
        meta = pd.read_csv(sys.argv[4], header=0, sep='\t', quotechar='"', error_bad_lines=False)
        meta_subset=meta[["ID", "REF", "ALT"]]
        sum_stats=pd.concat([meta_subset, es, lp], axis=1, ignore_index=True)
        sum_stats.columns =['SNP', 'A1', 'A2', 'BETA', 'P']
        sum_stats['P']=sum_stats['P']*(-1)
        sum_stats['P'] = 10 ** sum_stats['P']
        sum_stats =sum_stats[sum_stats["P"]<=0.05]
        sum_stats.to_csv(sys.argv[1],sep='\t', index=False, header=True)
        END_SCRIPT

        python3 PRScs/PRScs.py \
        --ref_dir=~{ref_dir} \
        --bim_prefix=~{sub(bim_file_23,'\\.bim$','')} \
        --sst_file=~{formatted_sum_stats_file_name} \
        --n_gwas=~{nGwas} \
        --out_dir=~{output_prefix}

        cat ~{output_prefix}*.txt > ~{concat_file_name}
        plink2 -bfile ~{sub(bed,'\\.bed$','')} -score ~{concat_file_name} 2 4 6 --rm-dup -out ~{mismatch_file_prefix}
        plink2 -bfile ~{sub(bed,'\\.bed$','')} -score ~{concat_file_name} 2 4 6 --exclude ~{mismatch_file_prefix}.rmdup.mismatch  -out ~{scored_file_prefix}
    >>>

    output {
        File output_formatted_sum_stats_file = "${formatted_sum_stats_file_name}"
        File output_concat_file = "${concat_file_name}"
        File output_mismatch_file = "${mismatch_file_prefix}.rmdup.mismatch"
        File output_scored_file = "${scored_file_prefix}.sscore"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
	disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}
