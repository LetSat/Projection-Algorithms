########### CameraProjections.py ###########
#
# Release 0.1
#
# First working release
#
# TODO: * Error Handling
#       * Modify function signatures to
#         reduce repetitive calculations
#       * Documentation
#       * Add camera parameter and cartesian
#         helper functions (Chandler)
#       * Verify range of output
# 
############################################


########### Imports and Declarations ###########

from math import *

# Earth's radius in km
EARTH_RADIUS = 6371.0


########### polarDistance ###########
# Discover the distance between the center of the earth in an image and an arbitrary point on the
# earth using polar coordinates

# 11-20-17
# Brian Howell

# Brief: Find the distance between the center of the earth in the image and the given point using
#        polar coordinates
#
# param[in]: alt - altitude of camera from surface of earth (km)
#            clat - lattitude of center of earth in image
#            clong - longitude of center of earth in image
#            plat - lattitude of point on earth to find polar distance to
#            plong - longitude of point on earth to find polar distance to
#
# return: magnitude - magnitude of vector from center of earth to point in the image scaled
#                     from 0 to 1, where 0 is center and 1 is edge of earth
#         angle - angle of polar vector from center of earth to point in the image scaled from
#                 -180 to 180
def polarDistance (alt, clat, clong, plat, plong):
    #print ""
    
    # Find vector from center of earth to center point
    earthToCenter = pointVector(clat, clong)
    #print "earthToCenter", earthToCenter
    
    # Find vector from center of earth to other point
    earthToPoint = pointVector(plat, plong)
    #print "earthToPoint", earthToPoint
    
    # Find the difference of the two vectors
    centerToPoint = diffVectors(earthToPoint, earthToCenter)
    #print "centerToPoint", centerToPoint
    
    # Find vector from center of image to camera FIXED
    #centerToCam = centerCamVector(alt + EARTH_RADIUS, earthToCenter)
    centerToCam = centerCamVector(alt, earthToCenter)
    #print "centerToCam", centerToCam
    
    # Find vector from point to camera
    pointToCam = diffVectors(centerToCam, centerToPoint)
    #print "pointToCam", pointToCam
    
    # Find angle between centerToCam and pointToCam
    pointAngle = abs(vectorAngleSimple(centerToCam, pointToCam))
    #print "pointAngle", pointAngle
    
    # Find angle between camera and edge point of viewing
    view = viewAngle(alt)
    #print "view", view
    
    # Find magnitude of final polar coordinate based on angles
    # See https://en.wikipedia.org/wiki/3D_projection#Perspective_projection
    magnitude = tan(pointAngle) / tan(view)
    
    #Find vector between center of Earth and point ALT: 90 LONG: DC
    earthToTop = [0,0,1]
    
    #Project centerToPoint and earthToTop to the viewing plane
    projCenterToPoint = projectVectorToPlane(earthToCenter, earthToPoint)
    #print "projCenterToPoint", projCenterToPoint
    projEarthToTop = projectVectorToPlane(earthToCenter, earthToTop)
    #print "projEarthToTop", projEarthToTop
    
    #Find the angle between the projected vectors
    pointAngle = getPointAngle(projEarthToTop,projCenterToPoint, earthToCenter)
    #print "pointAngle", pointAngle
    
    #Reference this angle from 90 for final polar angle
    #angle = radians(90) - pointAngle
    angle = degrees(pointAngle)
    
    return [magnitude, angle]

########### latLong ###########
# Discover the Latitude and Longitude of a point given the polar vector from the center of a satellite
# image to an arbitrary point on the image.

# Brief: Find the Latitude and Longitude of a point on an image of Earth given by a polar vector from
#        the center of the image.
# param[in]: alt - altitude of camera from surface of earth (km)
#            clat - lattitude of center of earth in image
#            clong - longitude of center of earth in image
#            pmag - magnitude of the vector to the point (0..1)
#            pang - angle of vector to the point (-180..180)
#
# return: plat - latitude of the given point
#         plong - Longitude of the given point
def latLong(alt, clat, clong, pmag, pang):
    
    pang = radians(pang+180)
    clat = radians(clat)
    clong = radians(clong)
    
    #Reference Vector for rotation
    refVector = [EARTH_RADIUS, 0, 0]
    
    # Find angle between camera and edge of earth in image
    viewAngle = asin(EARTH_RADIUS / (EARTH_RADIUS + alt))
    #print "viewAngle", degrees(viewAngle)

    pointAngle = atan(pmag*tan(viewAngle))
    
    #print "pointAngle", degrees(pointAngle)
    
    pointX = xCoord(alt,pointAngle)
    
    pointY = (pointX - EARTH_RADIUS - alt)*tan(pointAngle)
    
    refPoint1 = [pointX, pointY, 0]
    #print "refPoint1", refPoint1
        
    refPoint2 = rotateVectorX(refPoint1, pang)
    #print "refPoint2", refPoint2
    
    refPoint3 = rotateVectorY(refPoint2, -clat)
    #print "refPoint3", refPoint3
    
    finalPoint = rotateVectorZ(refPoint3, clong)
    #print "finalPoint", finalPoint
    
    return vectorToLatLong(finalPoint)
    

########### polarDistance Dependencies ###########

def addVectors(vectorA, vectorB):
    
    x = vectorA[0] + vectorB[0]
    y = vectorA[1] + vectorB[1]
    z = vectorA[2] + vectorB[2]
    
    return [x, y, z]
    
def diffVectors(vectorA, vectorB):
    
    x = vectorA[0] - vectorB[0]
    y = vectorA[1] - vectorB[1]
    z = vectorA[2] - vectorB[2]
    
    return [x, y, z]

def projectVectorToPlane(centerNorm, point): #! Project point onto plane described by centerNorm
    
    mag = vectorMagnitude(centerNorm)
    cross = crossProduct(point, centerNorm)
    cross = [cross[0]/mag,cross[1]/mag,cross[2]/mag]
    cross = crossProduct(centerNorm, cross)
    cross = [cross[0]/mag,cross[1]/mag,cross[2]/mag]
    return cross

def projectVectorToVector(vectorA, vectorB): # Project a onto b
    
    unitB = unitVector(vectorB)
    a1 = dotProduct(vectorA, unitB)
    proj = [unitB[0]*a1, unitB[1]*a1, unitB[2]*a1]

    return proj

def centerCamVector(alt, earthToCenter):
    
    #Find unit vector in direction from center of earth to center of image
    unitVector = [0.0, 0.0, 0.0]
    unitVector[0] = earthToCenter[0] / EARTH_RADIUS
    unitVector[1] = earthToCenter[1] / EARTH_RADIUS
    unitVector[2] = earthToCenter[2] / EARTH_RADIUS
    
    #Find vector from center of image to camera
    centerToCam = [0.0, 0.0, 0.0]
    centerToCam[0] = unitVector[0] * alt
    centerToCam[1] = unitVector[1] * alt
    centerToCam[2] = unitVector[2] * alt
    
    return centerToCam

def vectorAngleSimple(vectorA, vectorB):
    cross = crossProduct(vectorA, vectorB)
    crossMag = vectorMagnitude(cross)
    dot = dotProduct(vectorA, vectorB)
    angle = atan2(crossMag, dot)
    return angle
    
def getPointAngle(vectorA, vectorB, planeNorm):
    cross = crossProduct(vectorA, vectorB)
    crossMag = vectorMagnitude(cross)
    
    dot = dotProduct(vectorA, vectorB)
    angle = atan2(crossMag, abs(dot))
    
    crossSign = vectorMagnitude(addVectors(cross, unitVector(planeNorm))) - crossMag
    
    if crossSign > 0:
        if dot < 0:
            angle = 3*pi/2 - angle
        else:
            angle = angle - 3*pi/2
    else:
        if dot < 0:
            angle = angle - pi/2
        else:
            angle = pi/2 - angle
            
        
    return angle
    
def viewAngle(alt):
    return asin(EARTH_RADIUS / (EARTH_RADIUS + alt))
    
def crossProduct(a, b):

    return [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]
    
def dotProduct(vectorA, vectorB):

    x = vectorA[0] * vectorB[0]
    y = vectorA[1] * vectorB[1]
    z = vectorA[2] * vectorB[2]
    
    return x + y + z
    
def unitVector(a):

    mag = vectorMagnitude(a)
    
    return [a[0]/mag, a[1]/mag, a[2]/mag]
    
def vectorMagnitude(vector):
    
    return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])


########### latLong dependencies ###########

def vectorToLatLong(vector):
    lat = asin(vector[2]/EARTH_RADIUS)
    long = atan2(vector[1], vector[0])
    
    return [degrees(lat),degrees(long)]

    
def rotateVectorZ(vector, angle):
    sine = sin(angle)
    cosine = cos(angle)
    
    x = vector[0] * cosine - vector[1] * sine
    y = vector[0] * sine   + vector[1] * cosine
    z = vector[2]
    
    return [x,y,z]


def rotateVectorY(vector, angle):
    sine = sin(angle)
    cosine = cos(angle)
    
    x =  vector[0] * cosine + vector[2] * sine
    y =  vector[1]
    z = -vector[0] * sine + vector[2] * cosine
    
    return [x,y,z]

def rotateVectorX(vector, angle):
    sine = sin(angle)
    cosine = cos(angle)
    
    x = vector[0]
    y = vector[1] * cosine - vector[2] * sine
    z = vector[1] * sine   + vector[2] * cosine
    
    return [x,y,z]

def xCoord(alt, pA):
    r = EARTH_RADIUS
    
    tang = tan(pA)
    tansq = tang*tang
    term1 = sqrt(r*r - 2*alt*r*tansq - alt*alt*tansq)
    x1 = (-term1 + tansq*(r+alt))/(tansq+1)
    x2 = (term1 + tansq*(r+alt))/(tansq+1)
    
    return x1 if x1 > x2 else x2


########### Shared dependencies ###########


# Brief: Find the 3D vector between the center of the earth in the image and the given point using
#        x,y,z coordinates
#
# param[in]: lat - lattitude of point on earth to find vector for
#            long - longitude of point on earth to find vector for
#
# param[out]: x - the x component of the resulting vector
#             y - the y component of the resulting vector
#             z - the z component of the resulting vector
def pointVector(lat, long):
    
    lat = radians(lat)
    long = radians(long)
    x = EARTH_RADIUS * cos(lat) * cos(long)
    y = EARTH_RADIUS * cos(lat) * sin(long)
    z = EARTH_RADIUS * sin(lat)
    
    return [x,y,z]
