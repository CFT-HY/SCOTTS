import math, sys, string

last = 0
last2 = 0

i = 0

for line in sys.stdin.readlines():
    bits = string.split(string.strip(line))
    this = int(bits[8])
    print i,(this - last2)/2.0
    last2 = last
    last = this
    i = i +1
