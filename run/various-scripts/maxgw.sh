(
cd $1
for thing in $(ls -1 | grep power)
do
	echo $(python -c 'import sys; print float(sys.argv[-1][10:-4])*0.1' \
	    $thing) \
	    $(python ../max.py $thing)
done
) | sort -n

