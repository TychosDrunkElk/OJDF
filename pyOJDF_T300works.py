import numpy as np
from math import *
from pylab import *
from time import *

N=1000
M=25
M1=M-1
M2=M-2
N1=N-1
CSCP,XL,XNO=1.32,4.32,-1.92
TW,TCE=2.5,1.0
	
SC,PR=0.7,0.7
RO,XMU=1.176E-03,1.85E-04
EW,SB=50.3,2.32E+06
B=5.27E+7
#this should be user input, for testing just hardcoding the numbers
A=30.0

dt,dy=.001,.15
Q=69.0
YOE=.21
E=45.3
epsons=0.0
QA=1.0
M1=M-1
M2=M-2
N1=N-1
D=B/A
da=dt/dy
db=da/dy
ROS=1.180
tau=.1

dts=dt/CSCP*sqrt(2/(A*RO*XMU))*ROS*tau
beta=SB/sqrt(RO*XMU*A)
class ReturnValue(object):
	def __init__(self,T, F, YF, YO, U, W):
		
		self.T =T
		self.F =F
		self.YF =YF
		self.YO=YO
		self.U=U
		self.W=W

def main():
	T=np.zeros((N,25))
	F=np.zeros((N,25))
	YF=np.zeros((N,25))
	YO=np.zeros((N,25))
	U=np.zeros((N,25))
	W=np.zeros((N,25))
	YO[0, :]=YOE*np.ones(25)
	T[0, :]=np.ones(25)
	F[0,0]=-beta*exp(-EW)
	#U[0,:]=1.0
	
	F[0,0]=0.0
	print dts
	#for i in range(0,M-1):
	#	U[0,i]=i/24.0
	#F[0,:]=1
	for i in range(M1):
		U[0,i]=2*dy*i/M1
	
        F[0,1] =F[0,0]+dy*(U[0,1]+U[0,0])
	for J in range(0,M2):
        	F[0,J+2] =F[0,J] +2.*dy*(U[0,J+2] +4.*U[0,J+1] +U[0,J])/3.0
        i=0
	while i<N-1:
		#if F[i,0]<-.1:
		#	dtnew=.0001
		#else: 
		
		T[i+1,0],F[i+1,0]=solidheatup(T[i,:],F[i,:])
		new=gasphase(T,F,YF,YO,U,W,i)
		
		T=new.T
		F=new.F
		YF=new.YF
		YO=new.YO
		U=new.U
		W=new.W
		i=i+1
		
	#print F[0,:]
	#i=0
	#while i<N-1:
	#	dum=newstep(T,F,YF,YO,U,W,i)
	#	T=dum.T
	#	F=dum.F
	#	YF=dum.YF
	#	YO=dum.YO
	#	U=dum.U
	#	W=dum.W
	#	i=i+1		
		
	print 'done'
	return ReturnValue(T,F,YF,YO,U,W)
	
def solidheatup(T,F):
	Twnew=T[0]+((dts*(T[1]-T[0])/(dy*PR))+dts*XL*F[0])
	Fwnew=-beta*exp(-EW/Twnew)
	#print Twnew
	return (Twnew,Fwnew)

def gasphase(T,F,YF,YO,U,W,k):
	'''find steady state with new profile found with solidheatup.'''
	T[k+1,M-1]=1.0
        YO[k+1,M-1]=YOE
        YF[k+1,M-1]=0.0
        U[k+1,M-1]=1.0
        S=.2304*epsons/sqrt(A)
        W[k,0] =D*YF[k,0]*YO[k,0]/exp(E/T[k,0])
        
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

		T[k+1,I] =T[k,I] +db*((T[k,I+1]-2.*T[k,I]+T[k,I-1])/PR)+da*F[k,I]*(TR-TL)+dt*Q*W[k,I]
        	YF[k+1,I]=(YF[k,I]+db*(YF[k,I+1]-2.*YF[k,I]+YF[k,I-1])/SC +da*F[k,I]*(YFR-YFL))/(1.0+ dt*D*YO[k,I]*ARH)
        	YO[k+1,I]=(YO[k,I]+db*(YO[k,I+1]-2.*YO[k,I]+YO[k,I-1])/SC +da*F[k,I]*(YOR-YOL))/(1.0-dt*XNO*D*YF[k,I]*ARH)
        	U[k+1,I] =U[k,I] +da*F[k,I]*(UR-UL) +db* (U[k,I+1]-2.*U[k,I]+U[k,I-1])+dt*(T[k,I]-U[k,I]*U[k,I])
        
	
	U[k+1,0]=0.
	ACT1=E/T[k+1,0]
	ARH1=1/exp(ACT1)
	YF0J=F[k,0]*2.0*dy*(YF[k,0]-1.0)*SC +YF[k,1]
	
        YF[k+1,0]=(YF[k,0]+db*(YF[k,1]-2.*YF[k,0]+YF0J)/SC+.5*da*(YF[k,1]-YF0J)*F[k,0]) /(1.0 +dt*D*YO[k,0]*ARH1)
        #print YF[k+1,0]
        YO0J=F[k,0]*2.0*dy*YO[k,0]*SC +YO[k,1]
        YO[k+1,0]=(YO[k,0]+db*(YO[k,1]-2.*YO[k,0]+YO0J)/SC+.5*da*(YO[k,1]-YO0J)*F[k,0])/(1.0-XNO*dt*D*YF[k,0]*ARH1)
        W[k+1,0] = D*YF[k+1,0]*YO[k+1,0]*ARH1
        F[k+1,1] =F[k+1,0]+dy*(U[k+1,1]+U[k+1,0])
        for J in range(0,M2):
        	F[k+1,J+2] =F[k+1,J] +2.*dy*(U[k+1,J+2] +4.*U[k+1,J+1] +U[k+1,J])/3.0
        return ReturnValue(T,F,YF,YO,U,W)

def plotprofiles(datarray):
	ion()
	xdata=np.zeros(len(datarray[0,:]))
	for i in range(len(datarray[0,:])):
		xdata[i]=i*dy
	ydata=datarray[0,:]
	line,=plot(xdata,ydata)
	axis([0,5,-1,10])
	for i in range(len(datarray)):
		line.set_ydata(datarray[i,:])
		draw()
		
		print i
		
        	
        
if __name__=='__main__':
        test=main()
        plotprofiles(test.U)
