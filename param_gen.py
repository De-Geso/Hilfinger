import numpy as np

nfiles = 1
fpre = "kn"

# Initialize column values
xmin = 40
xmax = 40
ymin = 0.5
ymax = 10

incx = 1
incy = 1/16

nx = int(abs(xmax-xmin)/incx+1)
ny = int(abs(ymax-ymin)/incy+1)

xx = np.linspace(xmin, xmax, nx)
yy = np.linspace(ymin, ymax, ny)

counter = 0

fnum = 0
fname = f"{fpre}_{fnum}.txt"
f = open(fname, "w")

for i in range(nx):
	for j in range(ny):
		f.write(f"{xx[i]} {yy[j]}\n")
		counter += 1
		if (counter % np.ceil(nx*ny/nfiles) == 0 and nfiles > 1):
			fnum += 1
			fname = f"{fpre}_{fnum}.txt"
			f = open(fname, "w")
