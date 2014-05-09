import sys, string, math

totalv = 0.0
count = 0

for line in open(sys.argv[1],'r').readlines():
    bits = string.split(line)
    count = count + int(bits[-1])
    totalv = totalv + float(bits[-2])

for line in open(string.replace(sys.argv[1],'div','rot'),'r').readlines():
    bits = string.split(line)
    count = count + int(bits[-1])
    totalv = totalv + float(bits[-2])


print totalv/(float(count)*(1024.0*1024.0*1024.0))
    
