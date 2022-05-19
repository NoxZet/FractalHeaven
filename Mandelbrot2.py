from PIL import Image
import numpy
from time import time
from sys import argv, exit


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
    xs = numpy.linspace(x1, x2, xsteps, False, dtype=numpy.clongfloat)
    ys = numpy.linspace(y1, y2, ysteps, False, dtype=numpy.clongfloat)
    constants = xs + (ys[:, None])*1j
    ns = numpy.full(constants.shape, numpy.clongfloat(seed), dtype=numpy.clongfloat)
    res = numpy.zeros(constants.shape, dtype=int)
    for n in range(cyclesteps):
        res += numpy.less(cap, abs(ns))
        ns = ns**2 + constants
        ns[numpy.isnan(ns)] = 32
    return res

# Render array into a PIL image and save
def renderMandelbrot(mbArray, targetFile, pallete):
    width = len(mbArray[0])
    height = len(mbArray)
    outImg = Image.new("RGB", (width, height))
    for y in range(height):
        for x in range(width):
            outImg.putpixel((x, y), pallete[int(mbArray[y][x])])
    outImg.save(targetFile)


# Takes about 174 seconds on Ryzen 5 1600X @ 3.87 GHz
# time1 = time()
# renderMandelbrot(mandelbrotNumpy(-0.72, -0.62, 0.41, 0.51, 3000, 3000, 512), "mbNumpy.png", 512)
# time2 = time()
# print("Numpy took: ", time2-time1, "seconds")

scale = 1#0.00001698708 #1
xOff = -0.8148014326160998
yOff = 0.19100385298787578
cycles = 2048
frameSize = 800
for j in range(84):
    if j >= 0 and j % int(argv[1]) == int(argv[2]):
        print(j, 3 * scale / frameSize)
        mapped = []
        if j < 57:
            for i in range(256 * 8):
                if i < 256 * 7:
                    tint = i // 6
                    mapped.append((tint, 0, 0))
                elif i < 256 * 7.3:
                    tint = int((i - 256 * 7) * 3.333333)
                    mapped.append((255 - tint, 0, tint))
                elif i < 256 * 7.8:
                    tint = int((i - 256 * 7.3) * 2)
                    mapped.append((0, tint, 255 - tint))
                else:
                    tint = int((i - 256 * 7.8) * 2.5)
                    mapped.append((tint, 255 - tint, 0))
        else:
            reduction = (j - 56) / 200
            for i in range(256 * 8):
                if i < 256 * 7:
                    tint = i // 6
                    mapped.append((tint, 0, 0))
                elif i < 256 * 7.15:
                    tint = int((i - 256 * 7) * 3.333333)
                    mapped.append((255 - tint, 0, tint))
                elif i < 256 * (7.3 - reduction):
                    tint = 128 + int((i - 256 * 7.15) * (0.5 / (0.15 - reduction)))
                    mapped.append((255 - tint, 0, tint))
                elif i < 256 * (7.8 - reduction * 3):
                    tint = int((i - 256 * (7.3 - reduction)) / (0.5 - reduction * 2))
                    mapped.append((0, tint, 255 - tint))
                else:
                    tint = int((i - 256 * (7.8 - reduction * 3)) * 2.5)
                    mapped.append((tint, 255 - tint, 0))
        # print(mapped)
        renderMandelbrot(mandelbrotNumpy(-1.5 * scale + xOff, 1.5 * scale + xOff, -1.5 * scale + yOff, 1.5 * scale + yOff, frameSize, frameSize, int(cycles)), "mandelBrot5/mb" + str(j) + ".png", mapped)
    scale *= 0.6855
#renderMandelbrot(mandelbrot(-2.25, +0.75, -1.25, 1.25, 600, 480, 60), "mbOverall.png", 60)
#renderMandelbrot(mandelbrot(-1, -0.5, 0, 0.5, 400, 400, 80), "mb.png", 80)
#renderMandelbrot(mandelbrot(-0.72, -0.62, 0.41, 0.51, 500, 500, 200), "mbClose2.png", 200)
#renderMandelbrot(mandelbrot(-0.72, -0.62, 0.41, 0.51, 1500, 1500, 300), "mbClose2.png", 300)