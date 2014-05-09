# Peak of gravitational waves. Note, does not include 8G factor (or 4 pi?)
# G is 1
import sys, string, math

max = None

for line in open(sys.argv[-1],'r').readlines():
	bits = string.split(string.strip(line))
	if not len(bits) == 3:
		continue
	k = float(bits[0])
	pow = float(bits[1])

	if max == None:
		max = pow
	elif max < pow:
		max = pow

print max


