# params: Iter obj traits seed
if [ $# -lt 4 ]; then
    echo "Usage: $0 [Iter] [obj] [traits] [seed]"
    echo ""
    exit 0
fi

echo "# Started" >& 2
exp_unique_num=$(awk 'NR==1 {next} {print $1}' ../Merged_Data_5.csv | cut -d "," -f 4 | uniq | wc -l) #stores unique count of experiments 
traits=${3:-1} #num traits param   #64
seed=${4:-5} #set the seed 
K=10 # batch size/ number of folds


obj_best=0.99 # objective threshold
Iter=${1:-10000000} #number of iterations
obj=${2:-1} #the objective function

case "$obj" in
    1) obj_best=0.99   ;;
    2) obj_best=0.03   ;;
    3) obj_best=0.0005 ;;
    *) echo "Error: obj must be 1, 2, or 3" >&2; exit 1;;
esac


#Itername is just to name the directory the output goes into
Itername=$((Iter/1000000))


echo "Iter = $Iter">&2
mkdir -p "./seed_${seed}traits_${traits}_iter${Itername}M_Obj${obj}"
cp ./Run_shell_C.R "./seed_${seed}traits_${traits}_iter${Itername}M_Obj${obj}/"


for ((trait=1; trait<=traits; trait=trait*2)); do #every iteration it doubles traits variable
	for ((i=1; i<=$exp_unique_num; i+=10)); do 
		if [ $(($exp_unique_num - $i)) -lt 10 ]; then 
			echo "Almost done" >&2 
			last_fold=$exp_unique_num
			unset exp_unique_index
		else 
			last_fold=$(($i+$K-1))
		fi
		declare -a exp_unique_index
		for ((fold=i; fold<=last_fold; fold++)); do
			exp_index=$(awk 'BEGIN {FS=","} NR==1 {next} $4 == "'$fold'" {print NR-1}' ../Merged_Data_5.csv)
			arr=($exp_index)
			start=${arr[0]}
			end=${arr[-1]}
			current_exp=$(($(($fold-$i))*3))
			exp_unique_index[$current_exp]=$fold
			exp_unique_index[$(($current_exp+1))]=$start
			exp_unique_index[$(($current_exp+2))]=$end
		done
		echo "cd ./seed_${seed}traits_${traits}_iter${Itername}M_Obj${obj}; Rscript Run_shell_C.R ../../Merged_Data_5.csv $seed $trait $obj_best $Iter ${exp_unique_index[@]}"
	done
done
