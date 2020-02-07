from math import *
import numpy as np
from matrix import matrix

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
        u = U*U
        __Mv = matrix(np.zeros((2,4)))
        __Mv.set(0,0,u)
        __Mv.set(0,1,-E*U)
        __Mv.set(1,0,V)
        __Mv.set(1,1,-E*V)
        __Mv.set(2,0,N)
        __Mv.set(2,1,-E*N)
        __Mv.set(3,1,1.0)
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
        __T1.set(1,1,-(t+b)/2)
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
