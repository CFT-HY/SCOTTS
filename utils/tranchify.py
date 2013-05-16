# tranchify.py
#
# Slightly jury-rigged routine to compute (eg) peak
# of power spectrum as function of time.
import sys, string

rowent = int(sys.argv[1])
sys.stderr.write('Row: %d\n' % rowent)

times_re = {}
times_im = {}

for i in range(0,5000,100):
    times_re[i] = []
    times_im[i] = []

k_example = None
dt = 0.1

for item in sys.argv[2:]:
    therow = string.split(string.strip((open(item,'r').readlines())[rowent]))
    t = int(therow[4])
    k = float(therow[0])

    if k_example == None:
        k_example = k

    re = float(therow[1])
    im = float(therow[2])

    act_i = t % 5000
    
    times_re[act_i].append(re)
    times_im[act_i].append(im)

t0 = sum(times_re[0])/len(times_re[0])

for i in range(0,5000,100):
    print k_example*dt*i, (1.0/t0)*sum(times_re[i])/len(times_re[i]), sum(times_re[i])/len(times_re[i]), sum(times_im[i])/len(times_im[i]), i, len(times_re[i]), len(times_im[i])
