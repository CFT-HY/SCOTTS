# nucproc.py
#
# Outputs a pre-computed list of timesteps at which nucleation
# should be attempted. The number of timesteps (duration), the rate (beta),
# a normalisation factor (p0) and the timestep size (dt) are all adjustable.
import math, random, time, sys, string

# Initial time
t = 0.0

# Timestep size
dt = 0.1

# Nucleation rate
beta = 0.0125

# Number of timesteps
duration = 20000

# Normalisation
p0 = 0.01

# End timestep
end = dt*float(duration)

# Random seed
random.seed(time.time())

# Volume (may be considered a fudge factor)
vol = 1024*1024*1024





bublist = []

shift = 0

for i in range(duration):
	sys.stderr.write('%d %d\n' % (i, len(bublist)))


	prob = p0*math.exp(beta*(t-end))
	t = t + dt

	p_nobub = math.pow(1.0-prob,vol)


	chance = random.uniform(0,1)

	if chance > p_nobub:
		p_onebub = vol*math.pow(1.0-prob,vol-1)*prob

#		sys.stderr.write('i %d p_nobub = %lf\n' % (i, p_nobub))
		bublist.append(`i`)
		if len(bublist) == 10000:
			break

#		sys.stderr.write('chance %g prob %g p_nobub %g p_onebub %g twoormore %g\n' % (chance,prob,p_nobub,p_onebub,1.0-(p_nobub + p_onebub)))

		if chance > (p_nobub + p_onebub):
			p_twobub = vol*math.pow(1.0-prob,vol-2)*prob*prob

			if len(bublist) == 0:
				shift = i
			bublist.append(`i-shift`)
			if len(bublist) == 10000:
				break


			if chance > (p_nobub + p_onebub + p_twobub):
				bublist.append(`i-shift`)
				if len(bublist) == 10000:
					break



#			p_twobub = math.pow(1.0-prob,vol-2)*prob*prob

#		else:
#			print 'one bubble'	
		
	else:
		pass



for bub in bublist:
	print bub

# print string.join(bublist,',')
