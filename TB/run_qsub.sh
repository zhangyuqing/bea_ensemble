for curr_n in 20 40
do
	for curr_m in 0 3 5
	do
		for curr_v in {1..8}
		do
			curr_name=N${curr_n}m${curr_m}v${curr_v}
			qsub -P combat -j y -o logs/${curr_name}.txt -N ${curr_name} -cwd -b y -pe omp 8 qsub_simpipe_crossmod_syseval.qsub $curr_n $curr_m $curr_v
		done
	done
done
