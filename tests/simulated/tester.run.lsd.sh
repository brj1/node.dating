for i in $(seq 1 50); do
	~/Downloads/lsd-0.2/bin2/src/lsd -i HIV_${i}_rooted.tre -o HIV_${i}_rooted_lsd -d HIV_${i}_dates.txt -c -n 1
done
