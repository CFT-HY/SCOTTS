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
	if not len(sys.argv) == 3:
		sys.stderr.write('Usage: %s <output> <step>\n' % sys.argv[0])
		sys.stderr.write('Try step of 10 to start.\n')
		sys.exit(100)

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
	yvalues = array.array('d')


	try:
		tvalue.read(fh, 1)
		xvalues.read(fh, L)
		yvalues.read(fh, L)
	except EOFError:
		sys.stderr.write('Can\'t even read a single record!\n')
		sys.exit(100)

	record_offset = fh.tell()-integer_offset
	fh.seek(integer_offset)

	sys.stderr.write('Record is %d bytes\n' % record_offset)

	try:

		while True:

			if not (i % increment == 0):
				fh.seek(record_offset,1)
				i = i + 1
				continue

			tvalue = array.array('d')
			xvalues = array.array('d')
			yvalues = array.array('d')


			try:
				tvalue.read(fh, 1)
				xvalues.read(fh, L)
				yvalues.read(fh, L)
			except EOFError:
				sys.stderr.write('Run out of file between ' +
						 'records %d and %d\n'
						 % (i-increment,i))
				break

			t = tvalue.tolist()[0]

			xv = xvalues.tolist()
			yv = yvalues.tolist()
			
			
			if not(i % increment == 0):
				continue

			
			ofile = open('/tmp/animate.%d.temp' % os.getpid(),'w')

			for l in range(len(xv)):
				ofile.write('%lf %lf\n' % (xv[l], yv[l]))

			ofile.close()

				
			gnuplot_fh = os.popen('gnuplot -p','w')
			gnuplot_fh.write('set term png\n')
			gnuplot_fh.write('set output ' +
					 '\'/tmp/animate.%d.png\'\n' 
					 % os.getpid())
			gnuplot_fh.write('unset key\n')
			gnuplot_fh.write('set title \'t=%lf\'\n' % t)

			gnuplot_fh.write('set xlabel \'r/t\'\n')
			gnuplot_fh.write('set ylabel \'v\'\n')
				
			# suitable for velocity profile
			gnuplot_fh.write('set xrange [0:0.25]\n')
			gnuplot_fh.write('set yrange [0:0.0015]\n')
				
			# also suitable for velocity profile
			gnuplot_fh.write(('plot ' +
					 '\'/tmp/animate.%d.temp\' ' +
					 'u (($1-0.5)/(1+%lf)):2 ' +
					 'w lines\n') % (os.getpid(),t))

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

			i = i + 1
			
	except KeyboardInterrupt:
		sys.stderr.write('Giving up gracefully...\n')

	finally:
		try:
			os.remove('/tmp/animate.%d.temp' % os.getpid())
			os.remove('/tmp/animate.%d.png' % os.getpid())
			sys.stderr.write('Successfully removed temp files\n')
		except OSError:
			pass

if __name__ == '__main__':
    main()
