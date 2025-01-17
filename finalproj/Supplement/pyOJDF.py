import numpy as np
from math import *
from pylab import *
from time import *


N=4000
M=50
M1=M-1
M2=M-2
N1=N-1
CSCP,XL,XNO=1.32,4.32*.333,-1.92
TW,TCE=2.5,1.0
sig=1.38E-16
CP=111087600	
SC,PR=0.7,0.7
RO,XMU=0.34E-03,4.152E-04
EW,SB=50.3*30.0/100.0,2.32E+06
B=5.27E+7
#RO=AMBIENT DENSITY (G/CM3)
#XMU=AMBIENT VISCOSITY (G/CM-SEC)
#SB=PRE-EXPONENTIAL FACTOR, PYROLYSIS LAW (G/CM2-SEC)
#EW=ACTIVATION ENERGY ,     PYROLYSIS LAW (NONDIM.)


dt,dy=.0000001,.15/2
Q=69.0*.333
h=.001
E=45.3*30.0/100.0
epsons=1
epson=.9
QA=0.0
M1=M-1
M2=M-2
N1=N-1
ROS=1.180
tau=.1


class ReturnValue(object):
	#make an object to return: bundles all data together
	def __init__(self,T, F, YF, YO, U, W):
		self.T =T
		self.F =F
		self.YF =YF
		self.YO=YO
		self.U=U
		self.W=W
		
		

def main(A,YOE,time=0):
	#sets important variables based on input, Runs simulation
	#if time is set to 1, will run until steady state and find time to 
	#SS
	D=B/A
	dts=dt*CSCP*sqrt(2/(A*RO*XMU))*ROS*tau
	beta=SB/sqrt(RO*XMU*A)
	S2=epson*sig*(1000^3)*(1/CP)*sqrt(2/(RO*XMU*A))
	Xh=h/(CP*sqrt((RO*XMU*A)/2))
	
	T=np.zeros((N,M))
	F=np.zeros((N,M))
	YF=np.zeros((N,M))
	YO=np.zeros((N,M))
	U=np.zeros((N,M))
	W=np.zeros((N,M))
	YO[0, :]=YOE*np.ones(M)
	T[0, :]=np.ones(M)
	
	F[0,0]=0.0
	print dts
	
	U[0,:]=1
	
        F[0,1] =F[0,0]+dy*(U[0,1]+U[0,0])
	for J in range(0,M2):
        	F[0,J+2] =F[0,J] +2.*dy*(U[0,J+2] +4.*U[0,J+1] +U[0,J])/3.0
        i=0
	while i<N-1:
		#if F[i,0]<-.1:
		#	dtnew=.0001
		#else: 
		
		T[i+1,0],F[i+1,0]=solidheatup(T[i,:],F[i,:],dts,S2,Xh,beta)
		new=gasphase(T,F,YF,YO,U,W,i,dts,D,beta,YOE,A)
		if time==1:
			if new.W[i,:].max()>.8:
				print 'reaction has begun'
				print i
				print new.W[i,:].max()
				return
		T=new.T
		F=new.F
		YF=new.YF
		YO=new.YO
		U=new.U
		W=new.W
		i=i+1
	
	print 'done'
	return ReturnValue(T,F,YF,YO,U,W)
	
def solidheatup(T,F,dts,S2,Xh,beta):
	#the solid heat up computations
	Twnew=T[0]+((dts*(T[1]-T[0])/(dy*PR))-dts*S2*(2.0*(np.power(T[0],4.0))-(.33333**4)-1)-Xh*(T[0]-.33333)+dts*XL*F[0])
	
	Fwnew=-1*beta*exp(-EW/Twnew)
	return (Twnew,Fwnew)

def gasphase(T,F,YF,YO,U,W,k,dts,D,beta,YOE,A):
	#gas heat up computations.
	T[k+1,M-1]=1.0
        YO[k+1,M-1]=YOE
        YF[k+1,M-1]=0.0
        U[k+1,M-1]=1.0
        
        W[k,0] =D*YF[k,0]*YO[k,0]/exp(E/T[k,0])
        if W[k,0]<0:W[k,0]=0
        
        for Fi, I in zip(F[k,1:M1], range(1,M1)):
        	ACT=E/T[k,I]
        	ARH=exp(-1*ACT)
        	W[k,I]=D*YF[k,I]*YO[k,I]*ARH
        	if Fi >= 0:
        		TL =T[k,I-1]
        		TR =T[k,I]
        		YOL =YO[k,I-1]
        		YOR =YO[k,I]
        		YFL =YF[k,I-1]
        		YFR =YF[k,I]
        		UL =U[k,I-1]
        		UR =U[k,I]
        	else:
        		TL =T[k,I]
        		TR =T[k,I+1]
        		YOL =YO[k,I]
        		YOR =YO[k,I+1]
        		YFL =YF[k,I]
        		YFR =YF[k,I+1]
        		UL =U[k,I]
        		UR =U[k,I+1]

		T[k+1,I] =T[k,I] +(dts/(dy*dy))*((T[k,I+1]-2.*T[k,I]+T[k,I-1])/PR)+(dts/dy)*F[k,I]*(TR-TL)+dts*Q*W[k,I]
        	YF[k+1,I]=(YF[k,I]+(dts/(dy*dy))*(YF[k,I+1]-2.*YF[k,I]+YF[k,I-1])/SC +(dts/dy)*F[k,I]*(YFR-YFL))/(1.0+ dts*D*YO[k,I]*ARH)
        	YO[k+1,I]=(YO[k,I]+(dts/(dy*dy))*(YO[k,I+1]-2.*YO[k,I]+YO[k,I-1])/SC +(dts/dy)*F[k,I]*(YOR-YOL))/(1.0-dts*XNO*D*YF[k,I]*ARH)
        	if YF[k+1,I]<0 or YF[k+1,I]>1: YF[k+1,I]=YF[k,I]
        	if YO[k+1,I]<0 or YO[k+1,I]>1: YO[k+1,I]=YO[k,I]
		U[k+1,I] =U[k,I] +(dts/dy)*F[k,I]*(UR-UL) +(dts/(dy*dy))* (U[k,I+1]-2.*U[k,I]+U[k,I-1])+dts*(T[k,I]-U[k,I]*U[k,I])
        
	
	U[k+1,0]=0.
	ACT1=E/T[k+1,0]
	ARH1=1/exp(ACT1)
	YF0J=F[k,0]*2.0*dy*(YF[k,0]-1.0)*SC +YF[k,1]
	
        YF[k+1,0]=(YF[k,0]+(dts/(dy*dy))*(YF[k,1]-2.*YF[k,0]+YF0J)/SC+.5*(dts/dy)*(YF[k,1]-YF0J)*F[k,0]) /(1.0 +dts*D*YO[k,0]*ARH1)
        if YF[k+1,0]<0 or YF[k+1,0]>1: YF[k+1,0]=YF[k,0]
        YO0J=F[k,0]*2.0*dy*YO[k,0]*SC +YO[k,1]
        YO[k+1,0]=(YO[k,0]+(dts/(dy*dy))*(YO[k,1]-2.*YO[k,0]+YO0J)/SC+.5*(dts/dy)*(YO[k,1]-YO0J)*F[k,0])/(1.0-XNO*dts*D*YF[k,0]*ARH1)
        if YO[k+1,0]<0 or YO[k+1,0]>1: YO[k+1,0]=YO[k,I]
	W[k+1,0] = D*YF[k+1,0]*YO[k+1,0]*ARH1
        if W[k,I]<0:W[k,I]=0
        F[k+1,1] =F[k+1,0]+dy*(U[k+1,1]+U[k+1,0])
        for J in range(0,M2):
        	F[k+1,J+2] =F[k+1,J] +2.*dy*(U[k+1,J+2] +4.*U[k+1,J+1] +U[k+1,J])/3.0
        return ReturnValue(T,F,YF,YO,U,W)


def plotprofiles(dataobj, nstart):
	#plots an animation of the ignition process.
	ion()
	xsize=len(dataobj.T[0,:])
	xdata=np.zeros(xsize)
	for i in range(xsize):
		xdata[i]=i*dy
	Tydata=dataobj.T[0,:]
	YFydata=dataobj.YF[0,:]
	YOydata=dataobj.YO[0,:]
	Uydata=dataobj.U[0,:]
	Wydata=dataobj.W[0,:]
	line,=plot(xdata,Tydata)
	line2,=plot(xdata,YFydata)
	line3,=plot(xdata,YOydata)
	line4,=plot(xdata,Uydata)
	line5,=plot(xdata,Wydata)
	axis([0,5,-1,5])
	leg=legend(('T/Te','YF','YO','u/ue','w'), 'upper right')
	tsize=len(dataobj.T)
	for i in range(nstart,tsize):
		line.set_ydata(dataobj.T[i,:])
		line2.set_ydata(dataobj.YF[i,:])
		line3.set_ydata(dataobj.YO[i,:])
		line4.set_ydata(dataobj.U[i,:])
		line5.set_ydata(dataobj.W[i,:])
		draw()
		
		print i

	

def getcontour(xarr,yarr,zarr):
	#makes a contour plot for the temperature
	cplot=contourf(xarr,yarr,zarr,20,cmap=cm.jet)
	cbar=colorbar(cplot)
	cbar.ax.set_ylabel('T/Te')
	xlabel('stretch rate (a)')
	ylabel('YOE')
	title('Temperature map')
	show()

def findextinct():
	#get all temperature peaks at steady state for a range of yo 
	#and a values
	yoarr=.01+np.array(range(10))*.01
	aarr=np.array([.2,.4,.6,.8,1,2,4,6,8,10])
	Tpeak=np.zeros((10,10))
	etas=np.zeros((10,10))
	#get data... will take a very long time...
	for yo,i in zip(yoarr, range(len(yoarr))):
		for a,j in zip(aarr,range(len(aarr))):
			dataobj=main(a,yo)
			Tpeak[i,j]=dataobj.T[-1,:].max()
			n=np.where(dataobj.T[-1,:]==dataobj.T[-1,:].max())
			etas[i,j]=np.mean(n)*dy
	return (Tpeak, etas)

def timecompyo():
	#plot yo vs time to ignite- computations were solved seperatly and 
	#hard coded in
	yo=np.array([.1,.15,.2,.25,.3,.4,.5,.6,.7,.8,.9])
	a=10
	cdts=.000185397/(CSCP*sqrt(2/(a*RO*XMU))*ROS*tau)
	ntyo=cdts*np.array([2153, 1064,672,476,359,232,162,168,93,69,56])
	plot(yo,ntyo)
	xlabel('YO')
	ylabel('time')
	title('time to ignition: a=10')
	show()

def timecompa():
	a=np.array([1,2,3,3.5,4,4.5,5,6,7,8,9,10])
	cdts=np.zeros(len(a))
	for i in range(10):
		cdts[i]=.0000001
	for i in range(10,12):
		cdts[i]=.000001
	
	nta=cdts*np.array([40,89,182,950,1364,2682,3415,6100,8769,12329,1827,2272])
	plot(a,nta)
	xlabel('stretch rate')
	ylabel('time')
	title('time to ignition: YOE=0.3')
	show()
			
       
if __name__=='__main__':
        test=main(10,.3,time=1)
        plotprofiles(test,9000)
