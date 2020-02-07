from math import *
import numpy as np
from matrix import matrix
from parametricObject import parametricObject

class parametricPlane(parametricObject):

    def __init__(self,T=matrix(np.identity(4)),width=1.0,height=1.0,color=(0,0,0),reflectance=(0.0,0.0,0.0),uRange=(0.0,0.0),vRange=(0.0,0.0),uvDelta=(0.0,0.0)):
        super().__init__(T,color,reflectance,uRange,vRange,uvDelta)
        self.__width = width
        self.__height = height
    
    #   u*w
    #   v*h
    #    0
    #    1

    def getPoint(self,u,v,width,height):
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
