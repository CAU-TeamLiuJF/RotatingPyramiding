set -e

name="mao_edit_5m5f_first"
repeat=100
output="${name}_100_0816"

module load /apps/modulefile/program/r/4.3.3

Rscript run.R -g 5 -s $1 -r ${repeat} -p G5.txt -n ${name} -o ${output}

Rscript get_result.R -g 5 -s $1 -r ${repeat} -o ${output}

Rscript get_target_gene.R -g 5 -s $1 -r ${repeat} -o ${output}