echo "# Started" > list_of_commands.txt
exp_unique_num=$(awk 'NR==1 {next} {print $1}' ../../Merged_Data_5.csv | cut -d "," -f 4 | uniq | wc -l) #stores unique count of experiments 
traits=1 #num traits param   #64
seed=5 #set the seed 
K=10 # batch size/ number of folds
obj_best=0.99 # objective threshold
Iter=${1:-10000000} #number of iterations


#Iter=10000000000000    #10M^2
#Iter=10000000000000000000000000000   #9x3 0s + 0 in 10M^4 
#Iter=100000000000000000000000000000000000000000000000000000000  #10M^8
#Iter=10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 #10M^16
Iter=100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 #10M^32

#Itername is just to name the directory the output goes into
Itername="10M^32"




echo "Iter = $Iter">&2
for ((trait=1; trait<=traits; trait=trait*2)); do #every iteration it doubles traits variable
	for ((i=1; i<=$exp_unique_num; i=i+10)); do 
		if [ $(($exp_unique_num - $i)) -lt 10 ]; then 
			echo "Almost done"; 
			last_j=$exp_unique_num
			unset exp_unique_index
		else 
			last_j=$(($i+$K-1))
		fi
		declare -a exp_unique_index
		for ((j=i; j<=last_j; j++)); do
			exp_index=$(awk 'BEGIN {FS=","} NR==1 {next} $4 == "'$j'" {print NR-1}' ../../Merged_Data_5.csv)
			arr=($exp_index)
			start=${arr[0]}
			end=${arr[-1]}
			current_exp=$(($(($j-$i))*3))
			exp_unique_index[$current_exp]=$j
			exp_unique_index[$(($current_exp+1))]=$start
			exp_unique_index[$(($current_exp+2))]=$end
		done
		echo "cd ./seed_${seed}traits_${traits}_iter${Itername}; Rscript Run_shell_C.R ../../../Merged_Data_5.csv $seed $trait $obj_best $Iter ${exp_unique_index[@]}" >> list_of_commands.txt
	done
done
