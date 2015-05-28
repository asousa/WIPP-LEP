#!/bin/bash
# -------------------------------------------
# jh.sh
# 
# Checks the status of jobs in the queue.
# Usage: ./jh <prefix>
# 		- Finds 
#prefix=$1
#echo $prefix

me=`whoami`

if [ $# -eq 1 ]; then
	# stats=`qstat | grep $me | grep $1 | awk -F' ' '{print $5}'`
	stats=`qstat -u $me | grep $1 | awk -F' ' '{print $10}'`
else
	stats=`qstat -u $me | awk -F' ' '{print $10}'`
fi

len=`echo ${stats[@]} | wc -w`

ccount=$(grep -o "C" <<< "$stats" | wc -l)
rcount=$(grep -o "R" <<< "$stats" | wc -l)
qcount=$(grep -o "Q" <<< "$stats" | wc -l)

echo "Completed: " $ccount
echo "Running: " $rcount
echo "Queued: " $qcount

if [ $rcount -eq 0 ] && [ $qcount -eq 0 ]; then
	res=0
else
	res=1
fi

exit $res