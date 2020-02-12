from math import *
import numpy as np
from matrix import matrix
from parametricObject import parametricObject
from graphicsWindow import graphicsWindow
from parametricSphere import parametricSphere
from parametricTorus import parametricTorus
from wireMesh import wireMesh

class cameraMatrix:
    
    def __init__(self,UP,E,G,nearPlane=10.0,farPlane=50.0,width=640,height=480,theta=90.0):
        __Mp = self.__setMp(nearPlane,farPlane)
        __T1 = self.__setT1(nearPlane,theta,width/height)
        __S1 = self.__setS1(nearPlane,theta,width/height)
        __T2 = self.__setT2()
        __S2 = self.__setS2(width,height)
        __W2 = self.__setW2(height)

        self.__UP = UP.normalize()
        self.__N = (E - G).removeRow(3).normalize()
        self.__U = self.__UP.removeRow(3).transpose().crossProduct(self.__N.transpose()).normalize().transpose()
        self.__V = self.__N.transpose().crossProduct(self.__U.transpose()).transpose()
        self.__Mv = self.__setMv(self.__U,self.__V,self.__N,E)
        self.__C = __W2*__S2*__T2*__S1*__T1*__Mp
        self.__M = self.__C*self.__Mv

    def __setMv(self,U,V,N,E):
        __Mv = matrix(np.identity(4))
        __Mv.set(0,0,U.get(0,0))
        __Mv.set(1,0,V.get(0,0))
        __Mv.set(2,0,N.get(0,0))
        __Mv.set(0,1,U.get(1,0))
        __Mv.set(1,1,V.get(1,0))
        __Mv.set(2,1,N.get(1,0))
        __Mv.set(0,2,U.get(2,0))
        __Mv.set(1,2,V.get(2,0))
        __Mv.set(2,2,N.get(2,0))
        E = E.removeRow(3).transpose()
        __Mv.set(0,3,(-E*U).get(0,0))
        __Mv.set(1,3,(-E*V).get(0,0))
        __Mv.set(2,3,(-E*N).get(0,0))
        return __Mv

    def __setMp(self,nearPlane,farPlane):
        a = -(farPlane+nearPlane)/(farPlane-nearPlane)
        b = -(2*farPlane*nearPlane)/(farPlane-nearPlane)
        __Mp = matrix(np.zeros((4,4)))
        __Mp.set(0,0,nearPlane)
        __Mp.set(1,1,nearPlane)
        __Mp.set(2,2,a)
        __Mp.set(2,3,b)   
        __Mp.set(3,2,-1)
        return __Mp

    def __setT1(self,nearPlane,theta,aspect):
        t = nearPlane * tan((pi/180)*(theta/2))
        b = -t
        r = aspect * t
        l = -r
        __T1 = matrix(np.identity(4))
        __T1.set(0,3,-(r+l)/2)
        __T1.set(1,3,-(t+b)/2)
        return __T1

    def __setS1(self,nearPlane,theta,aspect):
        t = nearPlane * tan((pi/180)*(theta/2))
        b = -t
        r = aspect * t
        l = -r
        __S1 = matrix(np.identity(4))
        __S1.set(0,0,2/(r-l))
        __S1.set(1,1,2/(t-b))
        return __S1

    def __setT2(self):
        __T2 = matrix(np.identity(4))
        __T2.set(0,3,1.0)
        __T2.set(1,3,1.0)
        return __T2

    def __setS2(self,width,height):
        __S2 = matrix(np.identity(4))
        __S2.set(0,0,width/2)
        __S2.set(1,1,height/2)
        return __S2

    def __setW2(self,height):
        __W2 = matrix(np.identity(4))
        __W2.set(1,1,-1)
        __W2.set(1,3,height)
        return __W2

    def worldToViewingCoordinates(self,P):
        return self.__Mv*P

    def viewingToImageCoordinates(self,P):
        return self.__C*P

    def imageToPixelCoordinates(self,P):
        return P.scalarMultiply(1.0/P.get(3,0))

    def worldToImageCoordinates(self,P):
        return self.__M*P

    def worldToPixelCoordinates(self,P):
        return self.__M*P.scalarMultiply(1.0/(self.__M*P).get(3,0))

    def getUP(self):
        return self.__UP

    def getU(self):
        return self.__U

    def getV(self):
        return self.__V

    def getN(self):
        return self.__N

    def getMv(self):
        return self.__Mv

    def getC(self):
        return self.__C

    def getM(self):
        return self.__M


class parametricCircle(parametricObject):
    
    def __init__(self,T=matrix(np.identity(4)),radius=1.0,color=(0,0,0),reflectance=(0.0,0.0,0.0),uRange=(0.0,0.0),vRange=(0.0,0.0),uvDelta=(0.0,0.0)):
        super().__init__(T,color,reflectance,uRange,vRange,uvDelta)
        self.__radius = radius
    
    #   r*u*cos(v)
    #   r*u*sin(v)
    #       0
    #       1

    def getPoint(self,u,v):
        __P = matrix(np.ones((4,1)))
        __P.set(0,0,self.__radius*u*cos(v))
        __P.set(1,0,self.__radius*u*sin(v))
        __P.set(2,0,0.0)
        return __P

    def setRadius(self,radius):
        self.__radius = radius

    def getRadius(self):
        return self.__radius


class parametricCone(parametricObject):
    
    def __init__(self,T=matrix(np.identity(4)),height=1.0,radius=1.0,color=(0,0,0),reflectance=(0.0,0.0,0.0),uRange=(0.0,0.0),vRange=(0.0,0.0),uvDelta=(0.0,0.0)):
        super().__init__(T,color,reflectance,uRange,vRange,uvDelta)
        self.__height = height
        self.__radius = radius

    #   height*(1-u)/height)*radius*cos(v)
    #   height*(1-u)/height)*radius*sin(v)
    #               height*u
    #                   1

    def getPoint(self,u,v):
        __P = matrix(np.ones((4,1)))
        __P.set(0,0,(self.__height*(1-u)/self.__height)*self.__radius*cos(v))
        __P.set(1,0,(self.__height*(1-u)/self.__height)*self.__radius*sin(v))
        __P.set(2,0,(self.__height*u))
        return __P
    
    def setHeight(self,height):
        self.__height = height

    def getHeight(self):
        return self.__height

    def setRadius(self,radius):
        self.__radius = radius

    def getRadius(self):
        return self.__radius


class parametricCylinder(parametricObject):
    
    def __init__(self,
        T=matrix(np.identity(4)),height=1.0,radius=1.0,color=(0,0,0),reflectance=(0.0,0.0,0.0),uRange=(0.0,0.0),vRange=(0.0,0.0),uvDelta=(0.0,0.0)):
        super().__init__(T,color,reflectance,uRange,vRange,uvDelta)
        self.__height = height
        self.__radius = radius
    
    #   r*cos(v)
    #   r*sin(v)
    #     h*u
    #      1

    def getPoint(self,u,v):
        __P = matrix(np.ones((4,1)))
        __P.set(0,0,self.__radius*cos(v))
        __P.set(1,0,self.__radius*sin(v))
        __P.set(2,0,(self.__height*u))
        return __P

    def setHeight(self,height):
        self.__height = height

    def getHeight(self):
        return self.__height

    def setRadius(self,radius):
        self.__radius = radius

    def getRadius(self):
        return self.__radius


class parametricPlane(parametricObject):
    
    def __init__(self,T=matrix(np.identity(4)),width=1.0,height=1.0,color=(0,0,0),reflectance=(0.0,0.0,0.0),uRange=(0.0,0.0),vRange=(0.0,0.0),uvDelta=(0.0,0.0)):
        super().__init__(T,color,reflectance,uRange,vRange,uvDelta)
        self.__width = width
        self.__height = height
    
    #   u*w
    #   v*h
    #    0
    #    1

    def getPoint(self,u,v):
        __P = matrix(np.ones((4,1)))
        __P.set(0,0,self.__width*u)
        __P.set(1,0,self.__height*v)
        __P.set(2,0,0.0)
        return __P

    def setWidth(self,width):
        self.__width = width

    def getWidth(self):
        return self.__width

    def setHeight(self,height):
        self.__height = height

    def getHeight(self):
        return self.__height


#Set up constants required for the camera and the rendering process
#Near and far planes
NP = 10.0
FP = 50.0

#Image size
WIDTH = 2800
HEIGHT = 1600

#Image aspect
THETA = 45.0
ASPECT = WIDTH/HEIGHT

#Vector in the up direction
Px = 0.0
Py = 0.0
Pz = 1.0

#Position of camera
Ex = 40.0
Ey = 40.0
Ez = 85.0

#Gaze point
Gx = 0.0
Gy = 0.0
Gz = 0.0

P = matrix(np.ones((4,1)))
E = matrix(np.ones((4,1)))
G = matrix(np.ones((4,1)))

#Set up the up vector
P.set(0,0,Px)
P.set(1,0,Py)
P.set(2,0,Pz)

#Set up the viewing point
E.set(0,0,Ex)
E.set(1,0,Ey)
E.set(2,0,Ez)

#Set up the gaze point
G.set(0,0,Gx)
G.set(1,0,Gy)
G.set(2,0,Gz)

#Plane attributes
planeT = matrix(np.identity(4))
planeT.set(0,3,0.0)
planeT.set(1,3,-50.0)
planeCol = (255,255,0)
planeRef = (0.0,0.0,0.0)
planeWidth = 20.0
planeLength = 20.0

#Circle attributes
circleT = matrix(np.identity(4))
circleT.set(0,3,-40.0)
circleT.set(1,3,40.0)
circleCol = (0,255,255)
circleRef = (0.0,0.0,0.0)
circleRadius = 10.0

#Sphere attributes
sphereT = matrix(np.identity(4))
sphereCol = (255,0,0)
sphereRef = (0.0,0.0,0.0)
sphereRadius = 10.0

#Torus attributes
torusT = matrix(np.identity(4))
torusCol = (0,255,0)
torusRef = (0.0,0.0,0.0)
torusInnerRadius = 20.0
torusOuterRadius = 5.0

#Cone attributes
coneT = matrix(np.identity(4))
coneT.set(0,3,-40.0)
coneT.set(1,3,0.0)
coneCol = (0,0,255)
coneRef = (0.0,0.0,0.0)
coneHeight = 20.0
coneRadius = 10.0

#Cylinder attributes
cylinderT = matrix(np.identity(4))
cylinderT.set(0,3,40.0)
cylinderT.set(1,3,0.0)
cylinderCol = (255,0,255)
cylinderRef = (0.0,0.0,0.0)
cylinderHeight = 20.0
cylinderRadius = 10.0

window = graphicsWindow(WIDTH,HEIGHT) #Open a graphics window
camera = cameraMatrix(P,E,G,NP,FP,WIDTH,HEIGHT,THETA) #Set camera viewing system

plane = parametricPlane(planeT,planeWidth,planeLength,planeCol,planeRef,(0.0,1.0),(0.0,1.0),(1.0/10.0,1.0/10.0))
circle = parametricCircle(circleT,circleRadius,circleCol,circleRef,(0.0,1.0),(0.0,2.0*pi),(1.0/10.0,pi/18.0))
sphere = parametricSphere(sphereT,sphereRadius,sphereCol,sphereRef,(0.0,2.0*pi),(0.0,pi),(pi/18.0,pi/18.0))
cone = parametricCone(coneT,coneHeight,coneRadius,coneCol,coneRef,(0.0,1.0),(0.0,2.0*pi),(1.0/10.0,pi/18.0))
cylinder = parametricCylinder(cylinderT,cylinderHeight,cylinderRadius,cylinderCol,cylinderRef,(0.0,1.0),(0.0,2.0*pi),(1.0/10.0,pi/18.0))
torus = parametricTorus(torusT,torusInnerRadius,torusOuterRadius,torusCol,torusRef,(0.0,2.0*pi),(0.0,2.0*pi),(pi/18.0,pi/9.0))

window.drawSegments(wireMesh(plane,camera).getSegList(),planeCol)
window.drawSegments(wireMesh(circle,camera).getSegList(),circleCol)
window.drawSegments(wireMesh(sphere,camera).getSegList(),sphereCol)
window.drawSegments(wireMesh(cone,camera).getSegList(),coneCol)
window.drawSegments(wireMesh(cylinder,camera).getSegList(),cylinderCol)
window.drawSegments(wireMesh(torus,camera).getSegList(),torusCol)

window.saveImage("testImage.png")
window.showImage()