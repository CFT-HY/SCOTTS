# Prototype nucleation process
import math, random, time

t = 0.0

dt = 0.1

beta = 0.008

duration = 8000

end = dt*float(duration)

random.seed(time.time())

for i in range(duration):

	prob = math.exp(beta*(t-end))
	t = t + dt
	if random.uniform(0,1) < prob:
		print i # t  #  , prob, 1.0
	else:
		pass

#		print t, prob, 0.0
