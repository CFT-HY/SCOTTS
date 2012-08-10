import sys, string, math

nbins = 50

mink = None;
maxk = None;

vals = []

for line in sys.stdin.readlines():

    bits = string.split(string.strip(line))

    if len(bits) < 2:
        continue

    k = float(bits[0])
    power = float(bits[-1])

    if k < mink or mink == None:
        mink = k

    if k > maxk or maxk == None:
        maxk = k

    vals.append((k, power))

sys.stderr.write('Mink: %f, Maxk: %f\n' % (mink, maxk))

dk = (maxk - mink)/float(nbins)

binned = []

for i in range(nbins+1):
    binned.append(0.0)

for item in vals:
#     print item
    k = item[0]
    bin_number = int(math.floor((k - mink)/dk))
#     print 'bin', `bin_number`

    binned[bin_number] = binned[bin_number] + item[1]

thisk = mink
for thisbin in binned:
    print thisk, thisbin
    thisk = thisk + dk
