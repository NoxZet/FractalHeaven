from PIL import Image
import numpy
from time import time
import colorsys


# Gets magnitude of mandelbulb point on given coordinates with given parameters
def mandelbulbMag(x, y, z, cyclesteps, power, cap, seed=(0, 0, 0)):
    xv, yv, zv = seed
    radius = numpy.sqrt(xv ** 2 + yv ** 2 + zv ** 2)
    for n in range(cyclesteps):
        radiusPower = numpy.power(radius, power)

        phi = power
        if xv == 0:
            phi *= numpy.pi / 2
        else:
            phi *= numpy.arctan(yv / xv)

        theta = power
        if radius == 0:
            theta *= numpy.pi / 2
        else:
            theta *= numpy.arccos(zv / radius)

        xv = radiusPower * numpy.sin(theta) * numpy.cos(phi) + x
        yv = radiusPower * numpy.sin(theta) * numpy.sin(phi) + y
        zv = radiusPower * numpy.cos(theta) + z
        radius = numpy.sqrt(xv ** 2 + yv ** 2 + zv ** 2)
        if numpy.isnan(radius) or radius > cap:
            break
    return n


# Returns 3D array of points' magnitude. Mind the xsteps*ysteps*zsteps*cyclesteps complexity
def mandelbulb(x1, x2, y1, y2, z1, z2, xsteps, ysteps, zsteps, cyclesteps, power=8, cap=2.0, seed=(0, 0, 0)):
    complete = []
    for z in numpy.linspace(z1, z2, zsteps, False):
        print(z)
        plane = []
        complete.append(plane)
        for y in numpy.linspace(y1, y2, ysteps, False):
            row = []
            plane.append(row)
            for x in numpy.linspace(x1, x2, xsteps, False):
                row.append(mandelbulbMag(x, y, z, cyclesteps, power, cap, seed))
    return complete


# Returns magnitude of 3D vector
def vectorMagnitude(v1):
    return numpy.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)


# Casts ray from rayPoint, going in direction given by rayVector
# the rays go in straight line and only bounce in -rayVector direction
# starts at rayPoint+rayVector*rayStart until rayPoint+rayVector*rayEnd
# second step is rayPoint+rayVector*(rayStart+rayStep)
# bulbParams: tuple(cyclesteps, power, cap, [seed])
# threshold: magnitude of point needed to be taken into account
# intensity: sum of magnitudes needed for ray to "bounce"
def raycast(rayPoint, rayVector, rayStart, rayEnd, rayStep, bulbParams, threshold, intensity):
    rem = intensity
    collisions = {}
    for n in range(rayStart, rayEnd, rayStep):
        rayOffset = (rayVector[0]*n, rayVector[1]*n, rayVector[2]*n)
        collisionPoint = (rayPoint[0]+rayOffset[0], rayPoint[1]+rayOffset[1], rayPoint[2]+rayOffset[2])
        mag = mandelbulbMag(*collisionPoint, *bulbParams)
        if mag < threshold:
            continue
        collisions[vectorMagnitude(rayOffset)] = numpy.minimum(mag, rem)
        rem -= mag
        if rem <= 0:
            break
    return collisions


def renderMandelbrot(mbArray, targetFile, cycle):
    width = len(mbArray[0])
    height = len(mbArray)
    outImg = Image.new("L", (width, height))
    for y in range(height):
        for x in range(width):
            outImg.putpixel((x, y), int(mbArray[y][x]/cycle*256))
    outImg.save(targetFile)


# Multiplies
def tupleMpToInt(tupleIn, mp):
    outList = []
    for val in list(tupleIn):
        outList.append(int(val*mp))
    return tuple(outList)


def renderRaycast(rtArray, targetFile, low, hueRange, maxMag):
    width = len(rtArray[0])
    height = len(rtArray)
    outImg = Image.new("RGB", (width, height))
    for y in range(height):
        for x in range(width):
            hsv = ((rtArray[y][x][1] - low)/hueRange, 1, rtArray[y][x][0]/maxMag)
            rgb = colorsys.hsv_to_rgb(*hsv)
            outImg.putpixel((x, y), tupleMpToInt(rgb, 255))
    outImg.save(targetFile)


def raycastMandelbulb(xsteps, ysteps, zsteps, bulbsteps=6, angle=0):
    vectorStep = 200 // zsteps
    complete = []
    startPoint = (numpy.sin(angle)*4, 0, -numpy.cos(angle)*4)
    for y in numpy.linspace(-0.3, 0.3, ysteps, True):
        row = []
        complete.append(row)
        ysin = numpy.sin(y)/50
        for x in numpy.linspace(-0.3, 0.3, xsteps, True):
            res = raycast(startPoint, (numpy.sin(-angle+x)/50, ysin, numpy.cos(-angle+x)/50), 100, 300, vectorStep, (bulbsteps, 8, 2.0), bulbsteps-1, 15)
            if len(res) == 0:
                row.append((0, 0))
                continue
            vals = list(res.keys())
            row.append((len(vals), sum(vals)/float(len(vals))))
    return complete


# Redundant functions that render only the upper half and mirror it
def renderRaycastHalf(rtArray, targetFile, low, hueRange, maxMag):
    width = len(rtArray[0])
    height = len(rtArray)
    heightFull = height*2 - 1
    outImg = Image.new("RGB", (width, heightFull))
    for y in range(height):
        for x in range(width):
            hsv = ((rtArray[y][x][1] - low)/hueRange, 1, rtArray[y][x][0]/maxMag)
            rgb = colorsys.hsv_to_rgb(*hsv)
            rgbAdjusted = tupleMpToInt(rgb, 255)
            outImg.putpixel((x, y), rgbAdjusted)
            outImg.putpixel((x, heightFull-y-1), rgbAdjusted)
    outImg.save(targetFile)


# Needs generalization
def raycastMandelbulbHalf(xsteps, ysteps, zsteps, bulbsteps=6, angle=0):
    vectorStep = 200 // zsteps
    complete = []
    startPoint = (numpy.sin(angle)*4, 0, -numpy.cos(angle)*4)
    ysteps = (ysteps+1) // 2
    for y in numpy.linspace(-0.3, 0, ysteps, True):
        row = []
        complete.append(row)
        ysin = numpy.sin(y)/50
        for x in numpy.linspace(-0.3, 0.3, xsteps, True):
            res = raycast(startPoint, (numpy.sin(-angle+x)/50, ysin, numpy.cos(-angle+x)/50), 100, 300, vectorStep, (bulbsteps, 8, 2.0), bulbsteps-1, 15)
            if len(res) == 0:
                row.append((0, 0))
                continue
            vals = list(res.keys())
            row.append((len(vals), sum(vals)/float(len(vals))))
    return complete


for n in range(0, 20, 1):
    # renderRaycast(raycastMandelbulb(201, 201, 200, 9), "rtx"+str(n)+".png", 2.7, 1.3, 3)
    renderRaycastHalf(raycastMandelbulbHalf(201, 201, 200, 10, n/10), "rtx"+str(n)+".png", 2.7, 1.6, 3)


# renderRaycastHalf(raycastMandelbulbHalf(201, 201, 200, 10), "rtxHALF.png", 2.7, 1.6, 3)

# numpy.set_printoptions(threshold=100000)
# result = mandelbulb(-1, 1, -1, 0, -1, 1, 80, 40, 80, 7, 6, 2.0, (0.5, 0.2, 0.1))
# for n in range(80):
#     renderMandelbrot(result[n], "mb"+str(n)+".png", 6)