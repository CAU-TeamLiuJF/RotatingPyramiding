set -e

strategy_name=$1
repeat=5
output="20250202_100"

module load /apps/modulefile/program/r/4.3.3

Rscript run.R -r ${repeat} -s ${strategy_name} -o ${output}

Rscript get_result.R -r ${repeat} -s ${strategy_name} -o ${output}

Rscript get_target_gene.R -r ${repeat} -s ${strategy_name} -o ${output}