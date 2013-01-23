# Prototype nucleation process
import math, random, time, sys, string

t = 0.0

dt = 0.1

beta = 0.0125

duration = 20000

p0 = 0.01

end = dt*float(duration)

random.seed(time.time())

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
		if len(bublist) == 500:
			break

#		sys.stderr.write('chance %g prob %g p_nobub %g p_onebub %g twoormore %g\n' % (chance,prob,p_nobub,p_onebub,1.0-(p_nobub + p_onebub)))

		if chance > (p_nobub + p_onebub):
			if len(bublist) == 0:
				shift = i
			bublist.append(`i-shift`)
			if len(bublist) == 500:
				break
#			p_twobub = math.pow(1.0-prob,vol-2)*prob*prob

#		else:
#			print 'one bubble'	
		
	else:
		pass



#for bub in bublist:
#	print bub

print string.join(bublist,',')
