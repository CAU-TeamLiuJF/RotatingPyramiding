set -e
module load /apps/modulefile/program/r/4.3.3
repeat=100
output="20250202_100"
strategies=("Baseline" "Cascading" "NonCascading" "Rotating")

strategy_name=$1
for strategy_name in "${strategies[@]}"
do
#Rscript run.R -r ${repeat} -s ${strategy_name} -o ${output}
Rscript get_result.R -r ${repeat} -s ${strategy_name} -o ${output}
#Rscript get_target_gene.R -r ${repeat} -s ${strategy_name} -o ${output}
done
