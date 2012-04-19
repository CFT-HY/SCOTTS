#!/usr/bin/env python
#
# animate.py
#
# David Weir 2012
#
# quick and dirty animated plots using Tkinter and gnuplot
# not very pretty, and not the 'right' way to do it, but it works well
#
import sys, os, time, string, array

from Tkinter import *
import Image, ImageTk


def main():

	similarity = False
	animation = False

	if not len(sys.argv) in [3,4,5]:
		sys.stderr.write('Usage: %s <output> <step> [similarity] [animation]\n\n' % 
				 sys.argv[0])
		sys.stderr.write('Try step of 10 to start.\n\n')
		sys.stderr.write('The x-axis is in the frame of the grid, unless\n' +
				 '\'similarity\' is selected. This also changes\n' +
				 'the y-axis to be appropriate to shocks.\n\n')
		sys.stderr.write('Using \'animation\' will generate all the images\n' +
				 'then attempt to stitch into a movie with ffmpeg.\n' +
				 'The resulting movie will be movie-$PID.avi.\n\n')
		sys.stderr.write('Use the EXPERIMENTAL animation flag at your own risk!!!\n\n')
		sys.exit(100)

	if len(sys.argv) > 3:
		if not string.strip(string.lower(sys.argv[3])) in ['similarity','animation']:
			sys.stderr.write('Optional arguments are \'similarity\', \'animation\'\n')
			sys.exit(100)
		elif string.strip(string.lower(sys.argv[3])) == 'similarity':
			sys.stderr.write('-- similarity solution options enabled\n')
			similarity = True
		else:
			sys.stderr.write('-- generating animation\n')
			animation = True

	if len(sys.argv) > 4:
		if not string.strip(string.lower(sys.argv[4])) in ['similarity','animation']:
			sys.stderr.write('Optional arguments are \'similarity\', \'animation\'\n')
			sys.exit(100)
		elif string.strip(string.lower(sys.argv[4])) == 'similarity':
			sys.stderr.write('-- similarity solution options enabled\n')
			similarity = True
		else:
			sys.stderr.write('-- generating animation\n')
			animation = True

	scaling = False

        filename = sys.argv[1]
	rate = 1
	dt = 0.01
	increment = int(sys.argv[2])

	root = Tk()

	thelabel = Label(root)
	thelabel.pack()

	fh = open(filename, mode='rb')

	i = 0

	lvalue = array.array('i')
	lvalue.read(fh,1)
	L = lvalue.tolist()[0]
	sys.stderr.write('L=%d\n' % L)

	integer_offset = fh.tell()

	# How big is a record?
	tvalue = array.array('d')
	xvalues = array.array('d')
	y1values = array.array('d')
	y2values = array.array('d')
	y3values = array.array('d')
	y4values = array.array('d')


	try:
		tvalue.read(fh, 1)
		xvalues.read(fh, L)
		y1values.read(fh, L)
		y2values.read(fh, L)
		y3values.read(fh, L)
		y4values.read(fh, L)
	except EOFError:
		sys.stderr.write('Can\'t even read a single record!\n')
		sys.exit(100)

	record_offset = fh.tell()-integer_offset
	fh.seek(integer_offset)

	sys.stderr.write('Record is %d bytes\n' % record_offset)

	if animation:
		os.mkdir('/tmp/movietmp-%d/' % os.getpid())
		sys.stderr.write('Made /tmp/movietmp-%d/ for animation stills\n' % os.getpid())

	try:

		while True:

			if not (i % increment == 0):
				fh.seek(record_offset,1)
				i = i + 1
				continue

			tvalue = array.array('d')
			xvalues = array.array('d')
			y1values = array.array('d')
			y2values = array.array('d')
			y3values = array.array('d')
			y4values = array.array('d')


			try:
				tvalue.read(fh, 1)
				xvalues.read(fh, L)
				y1values.read(fh, L)
				y2values.read(fh, L)
				y3values.read(fh, L)
				y4values.read(fh, L)
			except EOFError:
				sys.stderr.write('Run out of file between ' +
						 'records %d and %d\n'
							 % (i-increment,i))
				break

			t = tvalue.tolist()[0]

			xv = xvalues.tolist()


			y1v = y1values.tolist()
			y2v = y2values.tolist()
			y3v = y3values.tolist()
			y4v = y4values.tolist()
			
					
			if not(i % increment == 0):
				continue

			
			ofile1 = open('/tmp/animate.%d.temp1' % os.getpid(),'w')
			ofile2 = open('/tmp/animate.%d.temp2' % os.getpid(),'w')
			ofile3 = open('/tmp/animate.%d.temp3' % os.getpid(),'w')
			ofile4 = open('/tmp/animate.%d.temp4' % os.getpid(),'w')

			for l in range(len(xv)):
				ofile1.write('%lf %lf\n' % (xv[l], y1v[l]))
				ofile2.write('%lf %lf\n' % (xv[l], y2v[l]))
				ofile3.write('%lf %lf\n' % (xv[l], y3v[l]))
				ofile4.write('%lf %lf\n' % (xv[l], y4v[l]))

			xmax = max(xv)

			ofile1.close()
			ofile2.close()
			ofile3.close()
			ofile4.close()

				
			gnuplot_fh = os.popen('gnuplot -p','w')
			gnuplot_fh.write('set term png\n')
			gnuplot_fh.write('set output ' +
					 '\'/tmp/animate.%d.png\'\n' 
					 % os.getpid())

			gnuplot_fh.write('set title \'t=%lf\'\n' % t)


				
			# suitable for velocity profile

			if similarity:
				gnuplot_fh.write('set xlabel \'x/t\'\n')
				gnuplot_fh.write('set xrange [0:1]\n')
				gnuplot_fh.write('set yrange [0:0.25]\n')
			else:
				gnuplot_fh.write('set xlabel \'x\'\n')
				gnuplot_fh.write('set xrange [0:%f]\n' % xmax)
				gnuplot_fh.write('set yrange [0:1]\n')
				
			# also suitable for velocity profile
			if similarity:
				gnuplot_fh.write(('plot ' +
						  '\'/tmp/animate.%d.temp1\' ' +
						  'u (($1-0.5)/(1+%lf)):($2/100) ' +
						  'w lines title \'E/100\', ' + 
						  '\'/tmp/animate.%d.temp2\' ' +
						  'u (($1-0.5)/(1+%lf)):($2/10) ' +
						  'w lines title \'P/10\', ' +
						  '\'/tmp/animate.%d.temp3\' ' +
						  'u (($1-0.5)/(1+%lf)):2 ' +
						  'w lines title \'v\', ' + 
						  '\'/tmp/animate.%d.temp4\' ' +
						  'u (($1-0.5)/(1+%lf)):($2/25) ' +
						  'w lines title \'phi/25\'\n') %
						 (os.getpid(),t,os.getpid(),t,
						  os.getpid(),t,os.getpid(),t))

			else:
				
				gnuplot_fh.write(('plot ' +
						  '\'/tmp/animate.%d.temp1\' ' +
						  'u 1:2 ' +
						  'w lines title \'E\', ' + 
						  '\'/tmp/animate.%d.temp2\' ' +
						  'u 1:2 ' +
						  'w lines title \'P\', ' +
						  '\'/tmp/animate.%d.temp3\' ' +
						  'u 1:2 ' +
						  'w lines title \'v\', ' +
						  '\'/tmp/animate.%d.temp4\' ' +
						  'u 1:($2/25.0) ' +
						  'w lines title \'phi/25\'\n') % 
						 (os.getpid(),os.getpid(),os.getpid(),os.getpid()))

			gnuplot_fh.close()

			time.sleep(dt)

#			im = Image.open('/tmp/animate.%d.png' 
#					% os.getpid())
#			new_im = im.convert('L')

#			tkimage = ImageTk.PhotoImage(new_im,
#						     palette=256)

			tkimage = ImageTk.PhotoImage(file='/tmp/animate.%d.png'
						     % os.getpid())

			thelabel.config(image=tkimage)
			
			root.update_idletasks()

			root.update()

			if animation:
				os.system(('cp /tmp/animate.%d.png ' +
					  '/tmp/movietmp-%d/animate.%05d.png') %
					  (os.getpid(), os.getpid(), i))
				sys.stderr.write('Written to ' +
						 '/tmp/movietmp-%d/animate.%05d.png\n' %
						 (os.getpid(), i))


			i = i + 1

		if animation:
			sys.stderr.write('Making a movie...\n')
#			os.system(('ffmpeg -r 10 -b 1000 -i /tmp/movietmp-%d/animate.%%05d.png '
#				   + 'movie-%d.mp4\n') %
#				  (os.getpid(), os.getpid()))

			os.system(('mencoder \"mf:///tmp/movietmp-%d/animate*png\" -mf ' +
				   'fps=10 -o movie-%d.avi -ovc lavc -lavcopts ' +
				   'vcodec=msmpeg4v2:vbitrate=800\n') %
				  (os.getpid(), os.getpid()))

			sys.stderr.write('Your movie should be in movie-%d.avi\n' %
					 os.getpid())
			
	except KeyboardInterrupt:
		sys.stderr.write('Giving up gracefully...\n')

	finally:
		try:
			os.remove('/tmp/animate.%d.temp1' % os.getpid())
			os.remove('/tmp/animate.%d.temp2' % os.getpid())
			os.remove('/tmp/animate.%d.temp3' % os.getpid())
			os.remove('/tmp/animate.%d.temp4' % os.getpid())
			os.remove('/tmp/animate.%d.png' % os.getpid())

			sys.stderr.write('Successfully removed temp files\n')

			if animation:
				os.system('rm -rvf /tmp/movietmp-%d/' % os.getpid())

		except OSError, e:
			sys.stderr.write('OS error: %s\n' % `e`)

if __name__ == '__main__':
    main()
