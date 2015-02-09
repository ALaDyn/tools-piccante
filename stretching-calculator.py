#!/usr/bin/python
import math

def stretchingFunction( xi,  ialpha):
    return  ialpha*math.tan(xi / ialpha)

def inverseStretchingFunction( x, ialpha):
    return ialpha*math.atan( x / ialpha)

def derivativeStretchingFunction( csi, ialpha):
    buf =  math.cos(csi/ialpha)
    return 1.0/(buf**2)

def myHelpFunction( strFactor,   ialpha):
    return strFactor - stretchingFunction(1,  ialpha)

def computeAlphaFromStrFactor(strFactor):
    alphaMin = 1.001*2.0/math.pi
    alphaMax = 20. * alphaMin
    difference = 2. * strFactor
    guess = 0.0
    result = 0.0

    while (math.fabs(difference) > (strFactor*0.00001)):
        guess = 0.5*(alphaMin + alphaMax)
        result = myHelpFunction(strFactor, guess)
        if (result < 0):
            alphaMin = guess
        else:
            alphaMax = guess
            difference = result
        
    return guess

def myHelpFunction2( deriv,   ialpha):
    return deriv - derivativeStretchingFunction(1,  ialpha)

def computeAlphaFromDerivative(deriv):
    alphaMin = 1.0000001*2.0/math.pi
    alphaMax = 2000. * alphaMin
    difference = 2. * deriv
    guess = 0.0
    result = 0.0

    while (math.fabs(difference) > (deriv*0.0001)):
        guess = 0.5*(alphaMin + alphaMax)
        result = myHelpFunction2(deriv, guess)
        if (result < 0):
            alphaMin = guess
        else:
            alphaMax = guess
            difference = result
        
    return guess

alpha = 1.0
mode = input("which mode (0: stretching factor, 1: dxmax/dxmin, 2: complete)? ");
if (mode==0):
    strFactor = 1.0*input("stretching factor: ");
    print "strFactor = ", strFactor;
    #Nx = 1.0*input("Number of points stretched grid: ");
        
    alpha = computeAlphaFromStrFactor( strFactor );
    print "alpha = " , alpha
    print " dxmax/dxmin = {0:.3f}".format(derivativeStretchingFunction(1, alpha))
    
elif (mode==1):
    deriv = 1.0*input("dxmax/dxmin: ");
    print "dxmax/dxmin = ", deriv
    alpha = computeAlphaFromDerivative( deriv );
    strFactor = stretchingFunction( 1, alpha);
    print "strFactor = " , strFactor
   
elif (mode==2):
    xInput = 1.0*input("at what point do you want to fix the resolution?: ");
    deriv = 1.0*input("resolution at that point dxmax/dxmin: ");
    alpha = computeAlphaFromDerivative( deriv );
    Lx = 1.0*input("which size would you like?: ");
    strFactor = stretchingFunction(1, alpha);
    finalLxi = 1.0*(xInput/strFactor)*inverseStretchingFunction( strFactor*Lx/xInput, alpha)
    print "xi size = " , finalLxi
    print "strFactor = " , Lx / finalLxi
    print "final  dxmax/dxmin = {0:.3f}".format(derivativeStretchingFunction(finalLxi*strFactor/xInput, alpha))
    
if(mode!=2):
    print "do you want to increase the stretching area?"
    factor = input("percent factor (0: NO)? ");
    while (factor > 0):
        strFactor = stretchingFunction( 1.0 + factor*1.0/100, alpha);
        print "strFactor = " , strFactor
        print "final  dxmax/dxmin = {0:.3f}".format(derivativeStretchingFunction(1.0+factor/100.0, alpha))
        factor = input("continue? Tell me the new factor (0: NO)? ");
        
print "goodbye!"
