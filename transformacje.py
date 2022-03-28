# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:36:58 2022

@author: Paula
"""
from math import sin, cos, sqrt, atan, degrees


class Transformacje:
    def __init__(self, model: str = "wgs84"):
        #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = 2 * self.flattening - self.flattening ** 2
        print(model, self.b)
        
    def xyz2flh(self, X, Y, Z):
        
        import math as m
        a = 6378137.000
        e2 = 0.00669438002290
        l = m.atan2(Y,X)
        p = m.sqrt(X**2+Y**2) #zmienna pomocnicza
        f = m.atan(Z/(p*(1-e2)))
        
        while 1:
            N = Np(f,a, e2)
            h = ( p/m.cos(f) )- N
            f1 = f
            f = m.atan(Z/p *(1-(e2 * N / (N + h)))**(-1))
            if abs(f-f1) < (0.000001/206265):
                    break 
                
        return(f, l,h) #helip
    
    

if __name__=="__main__":
    #utworzenie obiektu
    geo = Transformacje(model = "grs80")
    #dane XYZ goecentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    fi, lam, h = geo.xyz2flh(X, Y, Z)
    print(fi, lam, h)