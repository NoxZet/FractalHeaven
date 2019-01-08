from PIL import Image
import numpy
from time import time


# Version running calculations on native Python using numpy just for ranges
def mandelbrot(x1, x2, y1, y2, xsteps, ysteps, cyclesteps):
    complete = []
    for y in numpy.linspace(y1, y2, ysteps, False):
        row = []
        for x in numpy.linspace(x1, x2, xsteps, False):
            c = x + y*1j
            z = 0
            for n in range(cyclesteps):
                z = z**2 + c
                if abs(z) > 2:
                    break
            row.append(n)
        complete.append(row)
    return complete


# Version using numpy, noticeably faster
def mandelbrotNumpy(x1, x2, y1, y2, xsteps, ysteps, cyclesteps, cap=2.0, seed=0.0):
    xs = numpy.linspace(x1, x2, xsteps, False, dtype=numpy.complex128)
    ys = numpy.linspace(y1, y2, ysteps, False, dtype=numpy.complex128)
    constants = xs + (ys[:, None])*1j
    ns = numpy.full(constants.shape, numpy.complex128(seed), dtype=numpy.complex128)
    res = numpy.zeros(constants.shape, dtype=int)
    for n in range(cyclesteps):
        res += numpy.less(2.0, abs(ns))
        ns = ns**2 + constants
        ns[numpy.isnan(ns)] = 32
    return res


# Render array into a PIL image and save
def renderMandelbrot(mbArray, targetFile, cycle):
    width = len(mbArray[0])
    height = len(mbArray)
    outImg = Image.new("L", (width, height))
    for y in range(height):
        for x in range(width):
            outImg.putpixel((x, y), int(mbArray[y][x]/cycle*256))
    outImg.save(targetFile)


# Takes about 174 seconds on Ryzen 5 1600X @ 3.87 GHz
time1 = time()
# renderMandelbrot(mandelbrotNumpy(-0.72, -0.62, 0.41, 0.51, 3000, 3000, 512), "mbNumpy.png", 512)
time2 = time()
print("Numpy took: ", time2-time1, "seconds")


#renderMandelbrot(mandelbrot(-2.25, +0.75, -1.25, 1.25, 600, 480, 60), "mbOverall.png", 60)
#renderMandelbrot(mandelbrot(-1, -0.5, 0, 0.5, 400, 400, 80), "mb.png", 80)
#renderMandelbrot(mandelbrot(-0.72, -0.62, 0.41, 0.51, 500, 500, 200), "mbClose2.png", 200)
#renderMandelbrot(mandelbrot(-0.72, -0.62, 0.41, 0.51, 1500, 1500, 300), "mbClose2.png", 300)