# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:36:58 2022

@author: Paula
"""
#from math import sin, cos, sqrt, atan, degrees,pi
import math as m
import numpy as np

############################################################################


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
        self.e2 = 2 * self.flattening - self.flattening ** 2
        print(model, self.b)
    
    def Np(self,f):
    
        '''
        Funkcja liczy promień krzywizny w pierwszym wertykale
        
        Argumenty:
            a = 6378137.000
            e2 = 0.00669438002290
        '''
        
        N = self.a/m.sqrt(1-(self.e2)*m.sin(f)**2)
        return(N)
    
    
    def XYZ2flh(self, X, Y, Z):
        
        '''
        Funkcja przelicza współrzędne kartezjańskie geocentryczne XYZ na 
        współrzędne krzywoliniowe fi (szerokosc geodezyjna), lambda (długosc geodezyjna), h (wysokosc elipsoidalna)
        
        Argumenty:
            X,Y,Z - współrzędne kartezjańskie (w metrach)        / typ: float lub int
        
        Wyniki:
            f,l,h - współrzędne krzywoliniowe       / typ: float lub int
            f, l (w radianach)
            h {elipsoidalne} ( w metrach)
        '''
    
        l = m.atan2(Y,X)
        p = m.sqrt(X**2+Y**2) #zmienna pomocnicza
        f = m.atan(Z/(p*(1-self.e2)))
        
        def Np(self,f):
            N = self.a/m.sqrt(1-(self.e2)*m.sin(f)**2)
            return(N)
        
        while 1:
            N = Np(self,f)
            h = ( p/m.cos(f) )- N
            f1 = f
            f = m.atan(Z/p *(1-(self.e2 * N / (N + h)))**(-1))
            if abs(f-f1) < (0.000001/206265):
                    break 

        return(f, l,h) #helip
    
    
    def flh2XYZ(self, f,l,h):
        
        '''
        Funkcja przelicza współrzędne krzywoliniowe f (szerokosc geodezyjna), l (długosc geodezyjna), h (wysokosc elipsoidalna)
        na współrzędne kartezjańskie geocentryczne XYZ.
        
        Argumenty:
            f,l,h - współrzędne krzywoliniowe       / typ: float lub int
            f, l (w radianach)
            h {elipsoidalne} ( w metrach)
           
        Wyniki:
             X,Y,Z - współrzędne kartezjańskie (w metrach)        / typ: float lub int
        '''
        
        def Np(self,f):
            N = self.a/m.sqrt(1-(self.e2)*m.sin(f)**2)
            return(N)
        
        N = Np(self,f)
        X = (N+h)*m.cos(f)*m.cos(l)
        Y = (N+h)*m.cos(f)*m.sin(l)
        Z = (N*(1-self.e2)+h)*m.sin(f)
        
        return (X,Y,Z)
    
    
    def fl2neu(self, f, l):
        
        R = np.array([(-m.sin(f)*m.cos(l), -m.sin(l), m.cos(f) * m.cos(l)),
                      (-m.sin(f)*m.sin(l), m.cos(l), m.cos(f)*m.sin(l)),
                      (m.cos(f), 0, m.sin(f))])
        n = R[:,0]
        e = R[:,1]
        u = R[:,2]
        return(n, e, u)
   
   
    def rad2deg (self,B):
    
        '''
        Funkcja przelicza radiany na stopnie
            Argumenty:
                B - kąt w radianach                    / typ: float lub int
                
            Wyniki:
               deg - kąt przeliczony na stopnie        / typ: float lub int
               
        '''
        deg = B*180/(m.pi)
        return(deg)
    
    
    
    
    def  fl2u2000 (self,f,l,L0):
        '''
        Funkcja przelicza współrzędne krzywoliniowe na xy na płaszczyźnie Gaussa-Krugera 
        i następnie do układu płaskiego 2000
        
        Argumenty:
            f,l - współrzędne krzywoliniowe (w radianach)       
            f(szerokosc geodezyjna)                          / typ: float lub int
            l(długosc geodezyjna)                            / typ: float lub int
            
            L0 - południk osiowy zależy od układu do jakiego bedziemy przeliczac lub z jakiego przeliczamy 
                 (w radianch)
            
                 Dla układu 2000: 15stopni, 18stopni, 21stopni, 24stopnie
             
        Wyniki:
            x2000, y2000 - współrzędne w układzie płaskim 2000 (w metrach)         / typ: float lub int
            
        '''
        
        def sigma(self,f):
            A0 = 1-(self.e2/4)-((3*self.e2**2)/64)-((5*self.e2**3)/256);
            A2 = (3/8) *(self.e2 + ((self.e2**2)/4) + ((15*(self.e2**3))/128));
            A4 = (15/256) *(self.e2**2 +((3*(self.e2**3))/4));
            A6 = (35*(self.e2**3))/3072;
            
            si = self.a*((A0*f)-(A2*m.sin(2*f))+(A4*m.sin(4*f))-(A6*m.sin(6*f)));
            return (si)
        
        b2 = (self.a**2) *(1-self.e2);
        ep2 = (self.a**2 - b2)/b2;
        t = m.tan(f);
        n2 = ep2 * ((m.cos(f))**2);
        dl = l - L0;
       
        def Np(self,f):
            N = self.a/m.sqrt(1-(self.e2)*m.sin(f)**2)
            return(N)
        
        N = Np(self,f)
        si = sigma(self,f)
        
        xgk = si + (((dl)**2)/2) *N*m.sin(f)*m.cos(f)*(1+(((dl)**2)/12) *(m.cos(f))**2 * (5-t**2 +9*n2 +4*(n2)**2)+(((dl)**4)/360) *(m.cos(f))**4 * (61-58*t**2 + t**4 + 270*n2 - 330* n2 *t**2));
        ygk = dl*N*m.cos(f) * (1+((dl**2)/6) *(m.cos(f))**2 *(1-t**2+n2) + ((dl**4)/120) *(m.cos(f))**4 *(5-18*t**2 +t**4 + 14*n2 - 58*n2*t**2)); 
        
        return(xgk,ygk)
    
    
        def u2000(self,xgk,ygk,L0):
            
            x2000 = xgk*0.999923
            y2000 = ygk*0.999923 + (L0/3)*1000000 + 500000
            return(x2000,y2000)         


    def  fl2u1992 (self,f,l):
            '''
            Funkcja przelicza współrzędne krzywoliniowe na xy na płaszczyźnie Gaussa-Krugera 
            i następnie do układu płaskiego 1992
            
            Argumenty:
                f,l - współrzędne krzywoliniowe (w radianach)       
                f(szerokosc geodezyjna)                          / typ: float lub int
                l(długosc geodezyjna)                            / typ: float lub int
                
                L0 - południk osiowy zależy od układu do jakiego bedziemy przeliczac lub z jakiego przeliczamy 
                     (w radianch)
                
                     Dla układu 1992: L0 = 19stopni
                 
            Wyniki:
                x92, y92 - współrzędne w układzie płaskim 1992 (w metrach)         / typ: float lub int
                
            '''
            L0 = 19*m.pi/180
            
            def sigma(self,f):
                A0 = 1-(self.e2/4)-((3*self.e2**2)/64)-((5*self.e2**3)/256);
                A2 = (3/8) *(self.e2 + ((self.e2**2)/4) + ((15*(self.e2**3))/128));
                A4 = (15/256) *(self.e2**2 +((3*(self.e2**3))/4));
                A6 = (35*(self.e2**3))/3072;
                
                si = self.a*((A0*f)-(A2*m.sin(2*f))+(A4*m.sin(4*f))-(A6*m.sin(6*f)));
                return (si)
            
            b2 = (self.a**2) *(1-self.e2);
            ep2 = (self.a**2 - b2)/b2;
            t = m.tan(f);
            n2 = ep2 * ((m.cos(f))**2);
            dl = l - L0;
           
            def Np(self,f):
                N = self.a/m.sqrt(1-(self.e2)*m.sin(f)**2)
                return(N)
            
            N = Np(self,f)
            si = sigma(self,f)
            
            xgk = si + (((dl)**2)/2) *N*m.sin(f)*m.cos(f)*(1+(((dl)**2)/12) *(m.cos(f))**2 * (5-t**2 +9*n2 +4*(n2)**2)+(((dl)**4)/360) *(m.cos(f))**4 * (61-58*t**2 + t**4 + 270*n2 - 330* n2 *t**2));
            ygk = dl*N*m.cos(f) * (1+((dl**2)/6) *(m.cos(f))**2 *(1-t**2+n2) + ((dl**4)/120) *(m.cos(f))**4 *(5-18*t**2 +t**4 + 14*n2 - 58*n2*t**2)); 
            
            return(xgk,ygk)
        
            def u1992 (self,xgk,ygk):
                
                x = xgk*0.9993 - 5300000
                y = ygk*0.9993 + 500000
                return(x92,y92) 
            
     
    def odl_2D(self,X,Y,X2,Y2):
        
        '''
        Funkcja oblicza odległosc 2D między punktami
        
        Argumenty:
             X,Y - współrzędne kartezjańskie jednego punktu (w metrach)        / typ: float lub int
             X2,Y2 - współrzędne kartezjańskie drugiego punktu (w metrach)        / typ: float lub int

        Wynik:
            odl - odległosc 2D miedzy punktami (w matrach)        / typ: float lub int
       
        '''
    
        odl = m.sqrt((X2-X)**2+(Y2-Y)**2)
        return(odl)
    
 
    def odl_3D(self,X,Y,Z,X2,Y2,Z2):
        
        '''
        Funkcja oblicza odległosc 3D między punktami
        
        Argumenty:
             X,Y,Z - współrzędne kartezjańskie jednego punktu (w metrach)            / typ: float lub int
             X2,Y2,Z2 - współrzędne kartezjańskie drugiwgo punktu (w metrach)        / typ: float lub int

        Wynik:
            odl - odległosc 3D miedzy punktami (w metrach)        / typ: float lub int
       
        '''
        
        odl = m.sqrt((X2-X)**2+(Y2-Y)**2+(Z2-Z)**2)
        return(odl)
    
    
    def kat_azymutu (self,f,l,f2,l2):
        
        '''
        Funkcja oblicza kąt azymutu
        
        Argumenty:
            f - szerokosc kątowa pierwszego punktu (pozycji satelity) {w radianach}   / typ: float lub int
            l - długość kątowa pierwszego punktu (pozycji satelity) {w radianach}     / typ: float lub int
            l2 - długość kątowa drugiego punktu (pozycji anteny) {w radianach}        / typ: float lub int
            f2 - szerokość kątowa drugiego punktu  (pozycji anteny) {w radianach}     / typ: float lub int
        
        Wynik:
            Az - kąt azymutu (w radianach)     / typ: float lub int
            
        '''
        
        Az = m.atan((m.tan(l2-l))/(m.sin(f2)))
        return(Az)
    
    def kat_elewacyjny (self,f,l,f2,l2):
        
        '''
        Funkcja oblicza kąt elewacyjny
        
        Argumenty:
            f - szerokosc kątowa pierwszego punktu (pozycji satelity) {w radianach}   / typ: float lub int
            l - długość kątowa pierwszego punktu (pozycji satelity) {w radianach}     / typ: float lub int
            l2 - długość kątowa drugiego punktu (pozycji anteny) {w radianach}        / typ: float lub int
            f2 - szerokość kątowa drugiego punktu  (pozycji anteny) {w radianach}     / typ: float lub int
        
        Wynik:
            El - kąt elewacyjny (w radianach)     / typ: float lub int
            
        '''
        x = m.acos((m.cos(l2-l)*m.cos(f2)))
        
        El = m.atan((m.cos(x)-0.1513)/m.sin(x))
        
        return(El)
    
    
    

if __name__=="__main__":
    #utworzenie obiektu
    geo = Transformacje(model = "grs80")
    #dane XYZ goecentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    X2 = 3664940.510; Y2 = 1409153.58; Z2 = 5009571.167
    
    f, l, h = geo.XYZ2flh(X, Y, Z)
    print('fi = ',f,'lam = ', l, 'h = ', h)
    
    X,Y,Z = geo.flh2XYZ(f,l,h)
    print('X = ',X,'Y = ',Y,'Z = ',Z)
    
    n,e,u = geo.fl2neu(f,l)
    print('n = ',n,'e = ', e, 'u = ', u)
    
    lam_deg = geo.rad2deg (l)
    #print('lam_deg = ',lam_deg )

    L0 = int(lam_deg) *m.pi/180
    
    x2000,y2000 = geo.fl2u2000(f,l,L0)
    print('x2000 = ',x2000, 'y2000 = ',y2000)
    
    x92,y92 = geo.fl2u1992(f,l)
    print('x1992 = ',x92, 'y1992 = ',y92)
    
    odl_2D = geo.odl_2D(X,Y,X2,Y2)
    print('odległsć_2D = ',odl_2D)
    
    odl_3D = geo.odl_3D(X,Y,Z,X2,Y2,Z2)
    print('odległsć_3D = ',odl_3D)
    
    f2, l2, h2 = geo.XYZ2flh(X2, Y2, Z2)
    print('f2 = ',f2,'l2 = ', l2, 'h2 = ', h2)
    
    kat_az= geo.kat_azymutu(f,l,f2,l2)
    print('kat_azymutu = ', kat_az)

    kat_el= geo.kat_elewacyjny(f,l,f2,l2)
    print('kat_elewacyjny = ', kat_el)
