set -e

repeat=100
output="edit_all_100"

module load /apps/modulefile/program/r/4.3.3

for i in 1 2 3 4;
do
Rscript run.R -g 5 -s ${i} -r ${repeat} -p G5.txt -o ${output}
done
#Rscript run.R -g 5 -s 5 -r ${repeat} -p G5.txt -o ${output}
for i in 1 2 3 4;
do
Rscript get_result.R -g 5 -s ${i} -r ${repeat} -o ${output}
done
