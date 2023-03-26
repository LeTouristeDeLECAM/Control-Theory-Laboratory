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
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
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

def PID_RT(SP,PV,MAN,MV_MAN,MV_FF,K_C,T_I,T_D,alpha,Ts,MV_MAX,MV_MIN,MV,MV_P,MV_I,MV_D,E,MAN_FF=False,PVInit=0,method='EBD_EBD'  ): 
    
    """
    The function "PID_RT" needs to be included in a "for or while loop".
    
    :SP: SP (or SetPoint ) vector
    :PV: PV (or Process Value) vector
    :MAN : MAN (or Manual controller mode) vector (True or False)
    :MV_MAN: MV_MAN (or Manual value for MV) vector
    :MV_FF : MV_FF ( or Feedf ort;ard) vector
    :K_C: control ler gain
    :T_I : i ntegral time const ant [s]
    :T_D : derivative time constant [ s]
    :alpha: Tfd = alpha*Td where Tfd is t he derivative f ilter time const ant [s)
    :Ts: sampling period [s)
    :MV_MAX : maximum value for MV (used for saturation and anti wind-up)
    :MV_MIN: minimum value for MV (used for saturation and anti t,ind-up)
    :MV: MV (or Manipulated Value) vector
    :MV_P: MV_P (or Pr opotional part of MV) vector
    :MV_I: MV_I (or Integral part of MV} vect or
    :MV_D: MV_D (or Derivative part of MV) vector
    :E: E (or control Error) vector
    :Man_FF: Activated FF in manual mode (opt ional: default boolean value i s False)
    :PVInit: Initial value for PV (optional: default value is 0): used if PID_RT is ran f irst in t he squence and no value of PV is available yet.

    :method: discretisation method (optiona l : default value is 'EBD' )
        EBD-EBD: EBD for i ntegral action and EBD for derivative action
    
    
    The function "PID_RT" appends new values to the vectors "MV", "MV_P", "MV_I ", and "MV_D" .
    The appended values are based on the PID algorithm, the controller mode, and feedforward.
    Note that saturation of "MV" within the limits [MV_MIN MV_MAX) is impl emented with anti wlind-up.


    """    
    #print ('on excecute le PID')
    #E-Error
    if len(PV)==0 : 
        E.append(SP[-1]-PVInit)
    
    else : 
        E.append(SP[-1] - PV[-1])
        
        
        
        
    if MAN[-1] == False : 
        
        
            
        #Proportional-part
        MV_P.append(K_C*E[-1])
        
        
        #Itegral-part
        if len(MV_I)==0 : 
            MV_I.append(( K_C * Ts/ T_I)* E[-1])
        else : 
            MV_I.append(MV_I[-1]+(K_C * Ts/ T_I)*E[-1])
        
        
        #Derivative_part (Have to be filtred slide 194)
        T_FD = alpha*T_D 
        
        if len(MV_D)==0 and len(E)>1:
            MV_D.append((K_C * T_D/T_FD + Ts)*(E[-1] - E[-2]))
        elif len(E)==1:
            MV_D.append(0)
        else:
            MV_D.append((T_FD / (T_FD + Ts))*MV_D[-1]+(K_C * T_D/(T_FD + Ts))*(E[-1] - E[-2]))
            
        
        #feedForward
        if MAN_FF == True : 
            
            MV_feed_forward = MV_FF[-1]
            
        else:
            
            MV_feed_forward = 0

                  
            
            
        #windup 
        if MV_P[-1]+MV_I[-1]+MV_D[-1]>MV_MAX:
            MV_I[-1]=(MV_MAX-MV_P[-1]-MV_D[-1])
            
        elif MV_P[-1]+MV_I[-1]+MV_D[-1]<MV_MIN :
            MV_I[-1]=(MV_MIN-MV_P[-1]-MV_D[-1])
                        
        
        # MV_part
        MV_SUM = MV_P[-1] + MV_I[-1] + MV_D[-1] + MV_feed_forward
        
        MV.append(MV_SUM)
    
    else:
        MV_I.append(0)
        MV_P.append(0)
        MV_D.append(0)
        
        MV.append(MV_MAN[-1])
        
        
        
        # C'est faux !!!!!!!!!!!!!!!!!!!
        #if len(MV_MAN)==0:
        #    MV.append(0)
        #else :
        #    MV.append(MV_MAN[-1])
           
            
#----------------------------------------

def IMC_Tuning(Kp,T1,T2,theta,gamma,method='SOPDT'):
    """
    The function "IMC_Tuning" Calculate and return Kc, TauI, TauD, alpha optimized with the IMC_tuning grid 
    Kp : stactic gain
    T1 : First time constant
    T2 : Seconde time constant
    theta : Delay 
    gamma : agressivity [0,2...0,9] more close to 0,2 more its agressive
    method: Select which methode the IMC_tuning should use
        SOPDT:  Second order plus delay
        FOPDT: First order plus delay
    
    """
    
    if method=='SOPDT':
        
        # silde 186 Ligne B 
        TauOLP = max(T1,T2)
        TauC = gamma*TauOLP
        Kc = (T1+T2)/(TauC)
        Kc = Kc/Kp
        TauI = T1+T2
        TauD = (T1 * T2)/(T1 + T2)
        
        return Kc, TauI, TauD
    
    
    elif method=='FOPDT':
        #Slide 186 Ligne H
        TauOLP = T1
        TauC = gamma*TauOLP
        
        Kc = (T1+(theta/2))/(TauC+(theta/2))
        Kc = Kc/Kp
        TauI = T1+(theta/2)
        TauD = T1*theta/(2*T1+theta)
       
        return Kc, TauI, TauD

#---------------------------------------------------------------------



def FF_RT(DV,Kd,Kp,Tlead1,Tlag1,Tlead2,Tlag2,thetaD,thetaP,Ts,FF_OUT,PV_LL1,PV_LL2,Dv0=0):
          
    """
    The function "FF_RT" 
    DV:
    Kd:
    Kp:
    Tlead1:
    Tlag1:
    Tlead2:
    Tlag2:
    thetaD:
    thetaP:
    Ts: sampling period [s]
    FF_OUT:
    PV_LL1:
    PV_LL2:
    Dv0=0:
    
    """
    diff = DV[-1] - Dv0
    KFF = Kd/Kp
    theta_FF = max(0,thetaD-thetaP)

    LL_RT(DV,KFF,Tlead1,Tlag1,Ts,PV_LL1,PVInit=0,method='EBD')

    LL_RT(PV_LL1,1,Tlead2,Tlag2,Ts,PV_LL2,PVInit=0,method='EBD')

    Delay_RT(PV_LL2,theta_FF,Ts,FF_OUT)

    
    
