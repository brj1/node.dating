for i in $(seq 1 50); do
	java -classpath ~/random/tempest dr.app.tempest.Hacker HIV_${i}_rooted.tre HIV_${i}_rooted.nex
done
