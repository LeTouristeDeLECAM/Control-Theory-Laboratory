import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output
from package_DBR import Delay_RT



#-----------------------------------        
def LL_RT(MV,Kp,Tlead,Tlag,Ts,PV,PVInit=0,method='EBD'):
    
    """
    The function "LL_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :T: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretisation method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
        
    
    The function "LL_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretisation method.
    """    
    
    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+(Tlead/Ts))*MV[-1]-(Tlead/Ts)*MV[-2]))#MV[-1])
            elif method == 'EFD':
                PV.append((1-K)*PV[-1] + K*Kp*((Tlead/Ts)*MV[-1]+(1-(Tlead/Ts))*MV[-2]))#MV[-2])
            #elif method == 'TRAP':
             #   PV.append((1/(2*T+Ts))*((2*T-Ts)*PV[-1] + Kp*Ts*(MV[-1] + MV[-2])))            
            else:
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+(Tlead/Ts))*MV[-1]-(Tlead/Ts)*MV[-2]))
    else:
        PV.append(Kp*MV[-1])

#-----------------------------------

def PID_RT(SP,PV,MAN,MV_MAN,MV_FF,K_C,T_I,T_D,alpha,Ts,MV_MAX,MV_MIN,MV,MV_P,MV_I,MV_D,E,MAN_FF=False,PVInit=0,method='EBD_EBD'): #,activateFF=activateFF
    
    """
This function should be imported in a loop 
SP: set point 
PV: process value
PVInit : PV initial value is applied when PID_RT start to run
MV : manipulated value 
MV_P : P-part manipulted value 
MV_I : I-part manipulated value
MV_D : D-part of the manipulated value
E : control error

K_C : controller gain
T_I : integration time constant [s]
T_D : derivative time constant [s]
T_FD :  derivative filter time constant [s] = alpha*T_D 
Ts : sampling time 

TRACK_ON : output tracking
TRACK_VAL : output te-racking value

MAN : manual mode
MAN_FF : feedforward active in manual mode
MV_FF: MV feedforward qui correspond à la partie qu'on anticipe de la disturbance D
MV_MIN : minimum value for MV (set the saturation value and anti wind-up)
MV_MAX :maximum value for MV ( set the saturation value)
MV_MAN : MV value in manual mode 


    """    
    
    #E-Error
    if len(PV)==0 : 
        E.append(SP[-1]-PVInit)
    
    else : 
        E.append(SP[-1] - PV[-1])
        
        
        
        
    if MAN[-1] == False : 
        

        
        #P-part
        MV_P.append(K_C*E[-1])
        
        
        #I-part
        if len(MV_I)==0 : 
            MV_I.append(( K_C * Ts/ T_I)* E[-1])
        else : 
            MV_I.append(MV_I[-1]+(K_C * Ts/ T_I)*E[-1])
   
        
        #D_part
        T_FD = alpha*T_D
    
        if len(MV_D)==0 and len(E)>1:
            MV_D.append((K_C * T_D/T_FD + Ts)*(E[-1] - E[-2]))
        elif len(E)==1:
            MV_D.append(0)
        else:
            MV_D.append((T_FD / (T_FD + Ts))*MV_D[-1]+(K_C * T_D/(T_FD + Ts))*(E[-1] - E[-2]))
            
        #windup 
        if MV_P[-1]+MV_I[]>MV_MAX:
            MV_I.append(MV_MAX-MV_P[-1])
        elif MV_P[-1]+MV_I[]<MV_MIN :
            MV_I.append(MV_MIN-MV_P[-1])
         
            
        #feedForward
                        
        
        #if len(MV_FF)== 0 or not activateFF :
        #    MV_feed_forward = 0 
        #else:
        #    MV_feed_forward = MV_FF[-1]
        
                        
        # MV_part                
        MV_SUM = MV_P[-1] + MV_I[-1] + MV_D[-1] #+ MV_feed_forward
        if MV_SUM > MV_MAX :
            MV_SUM =MV_MAX
        elif MV_SUM < MV_MIN : 
            MV_SUM = MV_MIN
        MV.append(MV_SUM)
        
    else:
        MV_I.append(0)
        MV_P.append(0)
        MV_D.append(0)
        if len(MV_MAN)==0:
            MV.append(0)
        else :
            MV.APPEND(MV_MAN[-1])

            
            
#----------------------------------------

def IMC_Tuning(Kp,T1,T2,theta,gamma,method='SOPDT'):
    """
    using the method of the course 
    """
    
    if method=='SOPDT':
        
        TauOLP = max(T1,T2)
        TauC = gamma*TauOLP
        Kc = (T1+T2)/(TauC)
        Kc = Kc/Kp
        TauI = T1+T2
        TauD = (T1 * T2)/(T1 + T2)
        alpha = 0.25
        return Kc, TauI, TauD, alpha
    
    
    elif method=='FOPDT':
        TauOLP = T1
        TauC = gamma*TauOLP
        
        Kc = (T1+(theta/2))/(TauC+(theta/2))
        Kc = Kc/Kp
        TauI = T1+(theta/2)
        TauD = T1*theta/(2*T1+theta)
        alpha = 0.25
        return Kc, TauI, TauD, alpha
        

#---------------------------------------------------------------------

def FF_RT(DV,Kd,Kp,Tlead1,Tlag1,Tlead2,Tlag2,thetaD,thetaP,Ts,FF_OUT,PV_LL1,PV_LL2,Dv0=0):
    diff = DV[-1] - Dv0
    KFF = Kd/Kp
    theta_FF = max(0,thetaD-thetaP)
    LL_RT(DV,KFF,Tlead1,Tlag1,Ts,PV_LL1,PVInit=0,method='EBD')
    LL_RT(PV_LL1,1,Tlead2,Tlag2,Ts,PV_LL2,PVInit=0,method='EBD')
    Delay_RT(PV_LL2,theta_FF,Ts,FF_OUT)
#   print(FF_OUT)
    
