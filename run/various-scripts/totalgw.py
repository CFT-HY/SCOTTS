# Total energy in gravitational waves
# note, does not include (32 pi) factor. G is 1.
import sys, string, math

total = 0.0
notcounted = 0.0
count = 0
L = int(sys.argv[2])
vol = float(L*L*L)

for line in open(sys.argv[1],'r').readlines():
	bits = string.split(string.strip(line))
	count = count + int(bits[2])

	if float(bits[0]) > math.pi:
		notcounted = notcounted + float(bits[1])
		continue
	total = total + float(bits[1])

# Norm to volume
print total/vol
sys.stderr.write('%s: count: %d, total: %lf; not counted: %f\n' % (sys.argv[1],count,total/vol,notcounted/vol))
