import sys, string, math

prev_k = 0.000001

total = 0.0

for line in open(sys.argv[1]).readlines():
    bits = string.split(string.strip(line))
    
    this_k = float(bits[0])
    up = float(bits[1])

    dlogk = math.log(this_k) - math.log(prev_k)

    total = total + dlogk*up*this_k*this_k

    prev_k = this_k

print sys.argv[1][11:-4], total
