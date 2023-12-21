import math

def pr(x):
    p0 = 10.0
    a=1.0
    theta = 3*math.pi/4
    #res=1.5*p0*((a/x)*(a/x)-(a/x)*(a/x)*(a/x)*(a/x))
    #res = 0.5*p0*(1 - (a/x)*(a/x))+0.5*p0*(1 - 4*(a/x)*(a/x)+3*(a/x)*(a/x)*(a/x)*(a/x))*math.cos(2*theta)
    res = 0.5*p0*(1 - (a/x)*(a/x) + (1 -4*(a/x)*(a/x)+ 3*(a/x)*(a/x)*(a/x)*(a/x))*math.cos(2*theta))
    return res

def ptheta(x):
    p0 = 10.0
    a=1.0
    theta = 3*math.pi/4
    #res=0.5*p0*(2+(a/x)*(a/x)+3*(a/x)*(a/x)*(a/x)*(a/x))
    #res = 0.5*p0*(1 + (a/x)*(a/x))-0.5*p0*(1 + 3*(a/x)*(a/x)*(a/x)*(a/x))*math.cos(2*theta)
    res = 0.5*p0*(1 + (a/x)*(a/x) - (1 +3*(a/x)*(a/x)*(a/x)*(a/x))*math.cos(2*theta))
    return res

def ptay(x):
    #res=0
    p0 = 10.0
    a=1.0
    theta = 3*math.pi/4
    res = -0.5*p0*(1 +2*(a/x)*(a/x) - 3*(a/x)*(a/x)*(a/x)*(a/x))*math.sin(2*theta)
    return res

def p11(x):
    p0 = 10.0
    a=1.0
    theta = math.pi/4
    x = math.sqrt(x*x+x*x)
    #res=1.5*p0*((a/x)*(a/x)-(a/x)*(a/x)*(a/x)*(a/x))
    res = pr(x)*math.cos(theta)*math.cos(theta)+ptheta(x)*math.sin(theta)*math.sin(theta)-ptay(x)*math.sin(2*theta)
    return res

def p22(x):
    p0 = 10.0
    a=1.0
    theta = math.pi/4
    x = math.sqrt(x*x+x*x)
    #res=0.5*p0*(2+(a/x)*(a/x)+3*(a/x)*(a/x)*(a/x)*(a/x))
    res = pr(x)*math.sin(theta)*math.sin(theta)+ptheta(x)*math.cos(theta)*math.cos(theta)+ptay(x)*math.sin(2*theta)
    return res

def p12(x):
    #res=0
    p0 = 10.0
    a=1.0
    theta = math.pi/4
    x = math.sqrt(x*x+x*x)
    res = 0.5*(pr(x)-ptheta(x))*math.sin(2*theta)+ptay(x)*math.cos(2*theta)
    return res

xp = 20 / 500
yp = 20 / 500
for i in range(1,501):
	xtmp = xp * i
	ytmp = yp * i
	if xtmp < (1/math.sqrt(2)):
		print(str(math.sqrt(xtmp**2 + ytmp**2))+" "+str(0))
	else:
		print(str(math.sqrt(xtmp**2 + ytmp**2))+" "+str(p22(xtmp)))
