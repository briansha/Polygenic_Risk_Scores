version 1.0

## Version 02-14-2022
##
## Snapshot 21 - uses a single sum stats file.
## Snapshot 28 - similar to Snapshot 29 - but uses a python.py file and feeds arguments into the file instead of
##               using a HereDoc.
## Snapshot 29 - uses a tsv file containing GCP paths for multiple sum stats files.
##
## This workflow calculates polygenic risk scores for different summary statistics files.
##
## Cromwell version support - Successfully tested on v76
##
## Distributed under terms of the MIT License
## Copyright (c) 2021 Brian Sharber
## Contact <brian.sharber@vumc.org>

workflow polygenic_risk_scores {

    input {
        String plink2_docker = "briansha/plink2:polygenic"
        File sum_stats_files
        Boolean formatting_flag = true                     # Whether you want the sum_stats_file formatted or not.

        File bim_file_23
        File bed
        File bim
        File fam
        Int nGwas
        File ref_dir_zipped

    }
	Array[Array[String]] sum_stats_file_lines = read_tsv(sum_stats_files)

    if (formatting_flag) {
    	scatter(items in sum_stats_file_lines) {
            call Calculate_Scores {
              input:
                docker = plink2_docker,
                formatted_sum_stats_file_name = items[0] + ".txt",
                concat_file_name = items[0] + "_concat.txt",
                output_prefix = items[0],
                mismatch_file_prefix = items[0],
                scored_file_prefix = items[0],
                sum_stats_file = items[1],
                bim_file_23 = bim_file_23,
                bed = bed,
                bim = bim,
                fam = fam,
                nGwas = nGwas,
                ref_dir_zipped = ref_dir_zipped
            }
          }
      }

      if (!formatting_flag) {
          scatter(items in sum_stats_file_lines) {
            call Calculate_Scores_Already_Formatted {
              input:
                docker = plink2_docker,
                formatted_sum_stats_file_name = items[0] + ".txt",
                concat_file_name = items[0] + "_concat.txt",
                output_prefix = items[0],
                scored_file_prefix = items[0],
                sum_stats_file = items[1],
                bim_file_23 = bim_file_23,
                bed = bed,
                bim = bim,
                fam = fam,
                nGwas = nGwas,
                ref_dir_zipped = ref_dir_zipped
            }
          }
      }

    output {
          Array[File]? output_formatted_sum_stats_file = Calculate_Scores.output_formatted_sum_stats_file
          Array[File]? output_concat_file = Calculate_Scores.output_concat_file
          Array[File]? output_mismatch_file = Calculate_Scores.output_mismatch_file
          Array[File]? output_scored_file = Calculate_Scores.output_scored_file
          Array[File]? output_concat_file_already_formatted = Calculate_Scores_Already_Formatted.output_concat_file
          Array[File]? output_scored_file_already_formatted = Calculate_Scores_Already_Formatted.output_scored_file
    }

    meta {
    	author : "Brian Sharber"
        email : "brian.sharber@vumc.org"
        description : "This workflow calculates polygenic risk scores for different summary statistics files."
    }
}

task Calculate_Scores {
    input {
        File sum_stats_file
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

    command <<<
        unzip ~{ref_dir_zipped}
        git clone https://github.com/getian107/PRScs.git
        python3 - ~{sum_stats_file} ~{formatted_sum_stats_file_name} << 'END_SCRIPT'
        #!/usr/bin/python3
        # -*- coding: utf-8 -*-

        #setting up the environment
        import pandas as pd
        import os, sys

        #format the sum stats file
        sum_stats = pd.read_csv(sys.argv[1], compression='gzip', header=0, sep='\t', quotechar='"', error_bad_lines=False)
        sig_subset =sum_stats[sum_stats["p_value"]<=0.05]
        sig_subset2= sig_subset[["variant_id", "A1", "A2", "BETA", "p_value"]]
        sig_subset2=sig_subset2.rename(columns={'variant_id': 'SNP', 'p_value': 'P'})
        sig_subset2.to_csv(sys.argv[2], sep='\t', index=False, header=True)
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

task Calculate_Scores_Already_Formatted {
    input {
        File sum_stats_file
        File bim_file_23
        File bed
        File bim
        File fam
        Int nGwas
        String formatted_sum_stats_file_name    # name for the generated formatted sum stats file.
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

    command <<<
        unzip ~{ref_dir_zipped}
        git clone https://github.com/getian107/PRScs.git

        python3 PRScs/PRScs.py \
        --ref_dir=~{ref_dir} \
        --bim_prefix=~{sub(bim_file_23,'\\.bim$','')} \
        --sst_file=~{sum_stats_file} \
        --n_gwas=~{nGwas} \
        --out_dir=~{output_prefix}

        cat ~{output_prefix}*.txt > ~{concat_file_name}
        plink2 -bfile ~{sub(bed,'\\.bed$','')} -score ~{concat_file_name} 2 4 6 -out ~{scored_file_prefix}
    >>>

    output {
        File output_concat_file = "${concat_file_name}"
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
