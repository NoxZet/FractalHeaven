from PIL import Image
import numpy
import math
from time import time
import colorsys


# Can run on any vector. Needs manual casting.
def tupleMagnitude(tuple):
    return numpy.sqrt(numpy.sum(numpy.square(tuple)))


# Much faster for 3D vectors
def vectorMagnitude(v1):
    return math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)


def normalVector(v1, v2):
    return (v1[1]*v2[2] - v1[2]*v2[1],
            v1[2]*v2[0] - v1[0]*v2[2],
            v1[0]*v2[1] - v1[1]*v2[0])


def multiplyVector(v1, mp):
    return (v1[0]*mp, v1[1]*mp, v1[2]*mp)


def dotProductVectors(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]


def normalizeVector(v1):
    mag = vectorMagnitude(v1)
    return (v1[0]/mag, v1[1]/mag, v1[2]/mag)


def addVectors(v1, v2):
    return (v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])


def addVectorsMultiplied(v1, v2, mp2):
    return (v1[0]+v2[0]*mp2, v1[1]+v2[1]*mp2, v1[2]+v2[2]*mp2)


def subVectors(v1, v2):
    return (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])


def multiplyVectorComponents(v1, v2):
    return (v1[0]*v2[0], v1[1]*v2[1], v1[2]*v2[2])


# This function might seem heavy but it wins against numpy because of conversion overhead
def solveMatrix(depth, left, consts):
    if depth == 1:
        if left[0][0] == 0:
            return None
        return [consts[0]/left[0][0]]
    newLeft = []
    newConsts = []
    res = None
    lowest = None
    lowestVal = None
    for i in range(depth):
        if left[i][0] == 0:
            continue
        lowest = i
        lowestVal = left[i][0]
    for i in range(depth):
        if i == lowest:
            continue
        if left[i][0] == 0:
            newLeft.append(left[i][1:])
            newConsts.append(consts[i])
        else:
            newLeftSingle = []
            for j in range(1, depth):
                newLeftSingle.append(left[i][j]-left[lowest][j]/lowestVal*left[i][0])
            newLeft.append(newLeftSingle)
            newConsts.append(consts[i]-consts[lowest]/lowestVal*left[i][0])
        if len(newLeft) == depth-1:
            res = solveMatrix(depth-1, newLeft, newConsts)
    if res is None:
        return None
    newConst = consts[lowest]
    for j in range(1, depth):
        newConst -= res[j-1]*left[lowest][j]
    return [newConst/lowestVal] + res


# Line given by starting point
class PhysLine:
    def __init__(self, startPoint, v1, mp1=None, mp2=None):
        self.p1 = startPoint
        self.v1 = v1
        self.mp1 = mp1
        self.mp2 = mp2

    def isValidMp(self, mp):
        if self.mp1 is not None and self.mp1 > mp:
            return False
        if self.mp2 is not None and self.mp2 < mp:
            return False
        return True

    def getMpPoint(self, mp):
        return addVectorsMultiplied(self.p1[0], self.v1[0], mp)


class PhysPlane:
    def __init__(self, startPoint, v1, v2):
        self.p1 = startPoint
        self.v1 = v1
        self.v2 = v2

    # Returns None if there is no intersect
    # Returns (line.v1 multiplier, self.v1 multiplier, self.v2 multiplier) if there is collision within bounds
    # For faster bound checking, bounds are [0, 1] for both vectors
    def lineIntersect(self, line):
        algLeft = [[self.v1[0], self.v2[0], -line.v1[0]],
                   [self.v1[1], self.v2[1], -line.v1[1]],
                   [self.v1[2], self.v2[2], -line.v1[2]]]
        algRight = [line.p1[0] - self.p1[0],
                    line.p1[1] - self.p1[1],
                    line.p1[2] - self.p1[2]]
        result = solveMatrix(3, algLeft, algRight)
        if result is None: return None
        if 0 <= result[0] <= 1 and 0 <= result[1] <= 1 and line.isValidMp(result[2]):
            return (result[2], result[0], result[1])
        else: return None

    def getNormalOnPoint(self, collisionPoint):
        return normalizeVector(normalVector(self.v1, self.v2))

    def impactAngle(self, line, collisionPoint):
        n1 = self.getNormalOnPoint(collisionPoint)
        return math.asin(dotProductVectors(line.v1, n1)/(vectorMagnitude(line.v1)*vectorMagnitude(n1)))

    # Doesn't check bounds/collision
    # collisionPoint is unused but included so physics can be duck-typed :^)
    def bounceVector(self, v1, collisionPoint):
        n1 = self.getNormalOnPoint(collisionPoint)
        return addVectorsMultiplied(v1, n1, -2*dotProductVectors(v1, n1))


class PhysCircle(PhysPlane):
    def lineIntersect(self, line):
        result = super().lineIntersect(line)
        if result is None: return None
        if (result[1]*2-1)**2 + (result[2]*2-1)**2 <= 1:
            return result
        return None


class PhysSphere:
    def __init__(self, startPoint, radius):
        self.p1 = startPoint
        self.r1 = radius

    # Returns None if there is no intersect
    # Returns (line.v1 multiplier on collision, line.v1 multiplier on second collision)
    def lineIntersect(self, line):
        # solve to get t, multiplier for line vector, that is (l.p + l.v*t) will be intersection points
        # sqrt((s.p0 - l.p0 - l.v0*t)^2 + (s.p1 - l.p1 - l.v1*t)^2 + (s.p2 - l.p2 - l.v2*t)^2) = s.r^2
        # polynomial:
        # (l.v0^2 + l.v1^2 + l.v2^2) t^2
        # -2((s.p0 - p.p0)*v0 + (s.p1 - p.p1)*v1 + (s.p2 - p.p2)*v2) t
        # (s.p0 - p.p0)^2 + (s.p1 - p.p1)^2 + (s.p1 - p.p1)^2 - r^2
        degree2 = line.v1[0]**2 + line.v1[1]**2 + line.v1[2]**2
        degree1 = -2*((self.p1[0] - line.p1[0])*line.v1[0]
                     +(self.p1[1] - line.p1[1])*line.v1[1]
                     +(self.p1[2] - line.p1[2])*line.v1[2])
        degree0 = ((self.p1[0] - line.p1[0])**2
                  +(self.p1[1] - line.p1[1])**2
                  +(self.p1[2] - line.p1[2])**2 - self.r1**2)
        discrim = degree1**2 - 4*degree2*degree0
        if discrim < 0: return None
        root = math.sqrt(discrim)
        res = ((-degree1-root)/(2*degree2), (-degree1+root)/(2*degree2))
        return res if line.isValidMp(res[0]) else None
        '''try:
            roots = numpy.roots([degree2, degree1, degree0])
            if True in numpy.iscomplex(roots): return None
            res = (numpy.amin(roots), numpy.amax(roots))
            return res if line.isValidMp(res[0]) else None
        except: return None'''

    def getNormalOnPoint(self, collisionPoint):
        return normalizeVector(subVectors(collisionPoint, self.p1))

    def impactAngle(self, line, collisionPoint):
        n1 = self.getNormalOnPoint(collisionPoint)
        return math.pi/2 - abs(math.asin(dotProductVectors(line.v1, n1)/(vectorMagnitude(line.v1)*vectorMagnitude(n1))))

    # Collision point needed
    def bounceVector(self, v1, collisionPoint):
        n1 = self.getNormalOnPoint(collisionPoint)
        return addVectorsMultiplied(v1, n1, -2*dotProductVectors(v1, n1))


# Global lights are lights that are applied only when ray is leaving the scene
# Intensity is typically an rgb tuple 0.0 - 1.0
class LightGlobalCone:
    # Vector must be normalized
    def __init__(self, vector, maxAngle, intensity=(1, 1, 1), bounceMin=0):
        self.v1 = vector
        self.ang = maxAngle
        self.lint = intensity
        self.bounceMin = bounceMin

    def lightAngle(self, ray):
        return math.pi / 2 + math.asin(dotProductVectors(ray.v1, self.v1) / (vectorMagnitude(ray.v1) * vectorMagnitude(self.v1)))

    # Bounds is tuple that specifies what line.v1 multiplier should be rendered
    # Global lights only process if bounds is None
    def processLight(self, ray, bounds):
        if bounds is not None or ray.bounces<self.bounceMin: return None
        angle = self.lightAngle(ray)
        if angle >= self.ang: return None
        targetIntensity = 1 - angle/self.ang
        ray.lbaked = addVectors(ray.lbaked, multiplyVector(multiplyVectorComponents(ray.lint, self.lint), targetIntensity))


class LightGlobalHighlight:
    def __init__(self, intensity=(1, 1, 1), bounceMin=1, tinted=True):
        self.lint = intensity
        self.bounceMin = bounceMin
        self.tinted = tinted

    def processLight(self, ray, bounds):
        if ray.bounces<self.bounceMin: return None
        ray.lbaked = addVectors(ray.lbaked, multiplyVectorComponents(ray.lint, self.lint) if self.tinted else self.lint)


def standardProcessImpact(obj, ray, impact):
    if impact is None:
        impact = obj.lineIntersect(ray)[0]
    impactP = addVectorsMultiplied(ray.p1, ray.v1, impact)
    bounceV = obj.bounceVector(ray.v1, impactP)
    ray.applyColorBounce(obj.lsub)
    ray.p1 = impactP
    ray.v1 = bounceV
    ray.bounces += 1


class ObjectPlane(PhysPlane):
    def __init__(self, startPoint, v1, v2, color):
        self.p1 = startPoint
        self.v1 = v1
        self.v2 = v2
        self.lsub = color

    def processImpact(self, ray, impact=None):
        standardProcessImpact(self, ray, impact)


class ObjectCircle(PhysCircle):
    def __init__(self, startPoint, v1, v2, color):
        self.p1 = startPoint
        self.v1 = v1
        self.v2 = v2
        self.lsub = color

    def processImpact(self, ray, impact=None):
        standardProcessImpact(self, ray, impact)


class ObjectSphere(PhysSphere):
    def __init__(self, startPoint, radius, color):
        self.p1 = startPoint
        self.r1 = radius
        self.lsub = color

    # This is object method and not ray method so that each object can apply anything it wants to the rays
    # Modifies the ray and returns list of any extra rays that were generated (empty list if no rays were generated)
    # Returns None if there was no impact
    def processImpact(self, ray, impact=None):
        standardProcessImpact(self, ray, impact)


class Ray(PhysLine):
    def __init__(self, startPoint, v1, sourcePixel, maxBounces=1000, maxGenerations=1, generation=1):
        super().__init__(startPoint, v1, 0.000001)
        self.sourcePixel = sourcePixel
        self.bouncesMax = maxBounces
        # Used to limit ray splitting
        self.generation = generation
        self.maxGenerations = maxGenerations
        # TODO light processing to allow for limited distance light sources (create bounce list)
        # When ray passes a light sources, lint is multiplied by its intensity and added to baked
        self.lbaked = (0, 0, 0)
        self.lint = (1, 1, 1)
        self.bounces = 0

    def applyColorBounce(self, colorBounce):
        self.lint = (max(0, self.lint[0]*colorBounce[0]),
                     max(0, self.lint[1]*colorBounce[1]),
                     max(0, self.lint[2]*colorBounce[2]))
        return self.lint


class Scene:
    def __init__(self, objectStack=[], lightStack=[]):
        self.objectStack = objectStack
        self.lightStack = lightStack

    def addObjects(self, objects):
        try:
            for obj in objects:
                self.objectStack.append(obj)
        except:
            self.objectStack.append(objects)

    def addLights(self, lights):
        try:
            for obj in lights:
                self.lightStack.append(obj)
        except:
            self.lightStack.append(lights)


# Returns list of three ObjectPlane that form hexahedron on given vectors
def generateHex(beg, v1, v2, v3, color):
    sumV = addVectors(addVectors(addVectors(v1, v2), v3), beg)
    return [
        ObjectPlane(sumV, v1, v2, color),
        ObjectPlane(sumV, v1, v3, color),
        ObjectPlane(sumV, v2, v3, color),
        ObjectPlane(sumV, multiplyVector(v1, -1), multiplyVector(v2, -1), color),
        ObjectPlane(sumV, multiplyVector(v1, -1), multiplyVector(v3, -1), color),
        ObjectPlane(sumV, multiplyVector(v2, -1), multiplyVector(v3, -1), color)
    ]


# Camera angle is 2-tuple (rotation left-right, rotation up-down)
# Default facing is (0, 1, 0)
def raytraceScene(scene, size, cameraPoint, cameraAngle, cameraPixelAngle, maxBounces=5):
    pixels = []
    rayStack = []
    for y in range(size[1]):
        row = []
        pixels.append(row)
        angleVert = size[1]*cameraPixelAngle/2 - cameraPixelAngle*y + cameraAngle[1]
        for x in range(size[0]):
            row.append((0, 0, 0))
            angleHor = size[0]*cameraPixelAngle/2 - cameraPixelAngle*x + cameraAngle[0]
            zcam = math.sin(angleVert)
            ycam = math.cos(angleHor)*math.cos(angleVert)
            xcam = math.sin(angleHor)*math.cos(angleVert)
            rayStack.append(Ray(cameraPoint, (xcam, ycam, zcam), (x, y), maxBounces, 1))
    rays = len(rayStack)
    raysDone = 0
    for ray in rayStack:
        raysDone += 1
        if raysDone % 10000 == 0:
            print(raysDone/rays)
        obj: ObjectSphere
        while True:
            objectCollisions = []
            for obj in scene.objectStack:
                collision = obj.lineIntersect(ray)
                if collision is not None:
                    objectCollisions.append((collision[0], obj))
            if len(objectCollisions) == 0 or ray.bounces >= ray.bouncesMax:
                for light in scene.lightStack:
                    light.processLight(ray, None)
                pixels[ray.sourcePixel[1]][ray.sourcePixel[0]] = ray.lbaked
                break
            objectCollisions.sort(key=lambda pair: pair[0])
            for light in scene.lightStack:
                light.processLight(ray, (0, objectCollisions[0][0]))
            objectCollisions[0][1].processImpact(ray, objectCollisions[0][0])
    return pixels


def floatToIntRGB(rgb):
    return (max(0,min(255,int(rgb[0]*255))),
            max(0,min(255,int(rgb[1]*255))),
            max(0,min(255,int(rgb[2]*255))))


def renderRaytrace(pixels, size, targetFile):
    outImg = Image.new("RGB", (size[0], size[1]))
    for y in range(size[1]):
        for x in range(size[0]):
            outImg.putpixel((x, y), floatToIntRGB(pixels[y][x]))
    outImg.save(targetFile)
