(
cd $1
for thing in $(ls -1 | grep div-ps-step)
do
	echo $(python -c 'import sys; print float(sys.argv[-1][11:-3])*0.1' \
	    $thing) \
	    $(python ../vf.py $thing)
done
) | sort -n

