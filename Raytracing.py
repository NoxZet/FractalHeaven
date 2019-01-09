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
    return numpy.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)


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


# Line given by starting point
class PhysLine:
    def __init__(self, startPoint, v1, mp1=None, mp2=None):
        self.p1 = startPoint
        self.v1 = v1
        self.mp1 = mp1
        self.mp2 = mp2

    def isValidMp(self, mp):
        if self.mp1 is None or self.mp2 is None: return True
        return self.mp1 < mp < self.mp2

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
        try:
            result = numpy.linalg.solve(algLeft, algRight)
            if 0 <= result[0] <= 1 and 0 <= result[1] <= 1 and line.isValidMp(result[2]):
                return tuple(result)
            else: return False
        except: return False

    # Doesn't check bounds/collision
    # collisionPoint is unused but included so physics can be duck-typed :^)
    def bounceVector(self, v1, collisionPoint):
        n1 = normalizeVector(normalVector(self.v1, self.v2))
        return addVectorsMultiplied(v1, n1, -2*dotProductVectors(v1, n1))


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
        degree1 = -2((self.p1[0] - line.p1[0])*line.v1[0]
                    +(self.p1[1] - line.p1[1])*line.v1[1]
                    +(self.p1[2] - line.p1[2])*line.v1[2])
        degree0 = ((self.p1[0] - line.p1[0])**2
                  +(self.p1[1] - line.p1[1])**2
                  +(self.p1[2] - line.p1[2])**2 - self.r1**2)
        try:
            roots = numpy.roots([degree2, degree1, degree0])
            return (numpy.amin(roots), numpy.amax(roots))
        except: return None

    # Collision point needed
    def bounceVector(self, v1, collisionPoint):
        n1 = normalizeVector(subVectors(collisionPoint, self.p1))
        return addVectorsMultiplied(v1, n1, -2*dotProductVectors(v1, n1))


testline = PhysLine((1, 0, 0), (0, 1, 0))
testplane = PhysPlane((0, 0, -1), (4, 4, 0), (0, 0, 2))
print(testplane.bounceVector(testline.v1))