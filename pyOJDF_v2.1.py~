import numpy as np
from math import *
N=50
M=25
M1=M-1
M2=M-2
N1=N-1
CSCP,XL,XNO=1.32,4.32,-1.92
TW,TCE=2.5,1.0
	
SC,PR=0.7,0.7
RO,XMU=1.176E-03,1.85E-04
EW,SB=50.3,2.32E+06
B=5.27E+07
#this should be user input, for testing just hardcoding the numbers

dt,dy=.005,.15
Q=69.0
YOE=.2324
E=45.3
A=50.0
epsons=1.0
QA=0.0
M1=M-1
M2=M-2
N1=N-1
D=B/A
da=dt/dy
db=da/dy
class ReturnValue(object):
	def __init__(self,T, F, YF, YO, U, W):
		self.T =T
		self.F =F
		self.YF =YF
		self.YO=YO
		self.U=U
		self.W=W

def main():
	T=np.array([range(25),[2.5,3.3,4.1,4.9,5.8,6.6,7.1,6.7,5.4,4.2,3.2,2.5,2.0,1.6,1.4,1.3,1.2,1.1,1.0,1.0,1.0,1.0,1.0,1.0,1.0]])
	F=np.array([range(25),[-1.1,-.98,-.74,-.36,.13,.69,1.3,1.9,2.5,3.0,3.5,3.9,4.3,4.6,4.9,5.3,5.6,5.9,6.2,6.5,6.8,7.1,7.4,7.7,8.0]])
	YF=np.array([range(25),[.5125,.4552,.3920,.3240,.2536,.1850,.1243,.0807,.0565,.0402,.0279,.0190,.0126,.0082,.0052,.0033,.0020,.0012,.0007,.0004,.0002,.0001,.0001,.0,.0]])
	YO=np.array([range(25),[.0014,.0016,.0017,.0020,.0026,.0050,.0136,.0400,.0860,.1270,.1590,.1827,.1994,.2110,.2188,.2240,.2272,.2293,.2306,.2314,.2314,.2314,.2314,.2314,.2314]])
	U=np.zeros((2,25))
	W=np.zeros((2,25))
	
        #initial conditions
        U[-1,0]=0.0
        U[-1,M-1]=1.0
        efg='this is a test'
        for Fp,Fm,i in zip(F[-1,1:M-1], F[-1,0:M-2], range(1,M1)):
        	U[-1,i]=(Fp-Fm)/4./dy
        i=1
	while i<=N:
		dum=newstep(T,F,YF,YO,U,W)
		T=dum.T
		F=dum.F
		YF=dum.YF
		YO=dum.YO
		U=dum.U
		W=dum.W
		i=i+1
	
	print T[50,:]
		
		
	
	print 'done'
	return ReturnValue(T,F,YF,YO,U,W)
	
	
def newstep(T,F,YF,YO,U,W):
	dumT=T[-1,:]
	dumF=F[-1,:]
	dumYF=YF[-1,:]
	dumYO=YO[-1,:]
	dumU=U[-1,:]
	dumW=W[-1,:]
	TN=np.zeros(25)
        YON=np.zeros(25)
        YFN=np.zeros(25)
        UN=np.zeros(25)
        WN=np.zeros(25)
        FN=np.zeros(25)
        
        TN[M-1]=1.0
        YON[M-1]=YOE
        YFN[M-1]=0.0
        UN[M-1]=1.0
        S=.2304*epsons/sqrt(A)
        beta=SB/sqrt(RO*XMU*A)
        dumW[0] =D*dumYF[0]*dumYO[0]/exp(E/dumT[0])
        for Fi, I in zip(dumF[1:M1], range(1,M1)):
        	ACT=E/dumT[I]
        	ARH=exp(-1*ACT)
        	dumW[I]=D*dumYF[I]*dumYO[I]*ARH
        	if Fi >= 0:
        		TL =dumT[I-1]
        		TR =dumT[I]
        		YOL =dumYO[I-1]
        		YOR =dumYO[I]
        		YFL =dumYF[I-1]
        		YFR =dumYF[I]
        		UL =dumU[I-1]
        		UR =dumU[I]
        	else:
        		TL =dumT[I]
        		TR =dumT[I+1]
        		YOL =dumYO[I]
        		YOR =dumYO[I+1]
        		YFL =dumYF[I]
        		YFR =dumYF[I+1]
        		UL =dumU[I]
        		UR =dumU[I+1]
        	print ACT
        	TN[I] =dumT[I] +db*((dumT[I+1]-2.*dumT[I]+dumT[I-1])/PR)+da*dumF[I]*(TR-TL)+dt*Q*dumW[I]
        	YFN[I]=(dumYF[I]+db*(dumYF[I+1]-2.*dumYF[I]+dumYF[I-1])/SC +da*dumF[I]*(YFR-YFL))/(1.0+ dt*D*dumYO[I]*ARH)
        	YON[I]=(dumYO[I]+db*(dumYO[I+1]-2.*dumYO[I]+dumYO[I-1])/SC +da*dumF[I]*(YOR-YOL))/(1.0-dt*XNO*D*dumYF[I]*ARH)
        	UN[I] =dumU[I] +da*dumF[I]*(UR-UL) +db* (dumU[I+1]-2.*dumU[I]+dumU[I-1])+dt*(dumT[I]-dumU[I]*dumU[I])
        
	TN[0]=-EW/log(-1*dumF[0]/beta)
	UN[0]=0.
	ACT1=E/TN[0]
	ARH1=1/exp(ACT1)
	YF0J=dumF[0]*2.0*dy*(dumYF[0]-1.0)*SC +dumYF[1]
        YFN[0]=(dumYF[0]+db*(dumYF[1]-2.*dumYF[0]+YF0J)/SC+.5*da*(dumYF[1]-YF0J)*dumF[0]) /(1.0 +dt*D*dumYO[0]*ARH1)
        YO0J=dumF[0]*2.0*dy*dumYO[0]*SC +dumYO[1]
        YON[0]=(dumYO[0]+db*(dumYO[1]-2.*dumYO[0]+YO0J)/SC+.5*da*(dumYO[1]-YO0J)*dumF[0])/(1.0-XNO*dt*D*dumYF[0]*ARH1)
        WN[0] = D*YFN[0]*YON[0]*ARH1
        
        TG = TN[0] -CSCP*TCE +XL
        TIMAG = (-TN[1]+2.*TN[0] -.5*dumF[0]*TN[1]*PR*dy-Q*WN[0]*dy*dy*PR)/(1.-.5*dumF[0]*PR*dy)
        FN[0] = (((TN[1]-TIMAG)/(2.*dy)) - S*TN[0]**4.0+QA)/(-TG*PR)
        
        FN[1] =FN[0]+dy*(UN[1]+UN[0])
        for J in range(0,M2):
        	FN[J+2] =FN[J] +2.*dy*(UN[J+2] +4.*UN[J+1] +UN[J])/3.0
        
        T=np.vstack((T,TN))
        YF=np.vstack((YF,YFN))
        YO=np.vstack((YO,YON))
        U=np.vstack((U,UN))
        F=np.vstack((F,FN))
        W=np.vstack((W,WN))

        return ReturnValue(T,F,YF,YO,U,W)
        
        	
        
if __name__=='__main__':
        test=main()
        
        
