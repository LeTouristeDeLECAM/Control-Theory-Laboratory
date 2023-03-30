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
    :MV_FF : MV_FF ( or Feedfortward) vector
    :K_C: control ler gain
    :T_I : i ntegral time const ant [s]
    :T_D : derivative time constant [ s]
    :alpha: Tfd = alpha*Td where Tfd is t he derivative f ilter time const ant [s]
    :Ts: sampling period [s]
    :MV_MAX : maximum value for MV (used for saturation and anti wind-up)
    :MV_MIN: minimum value for MV (used for saturation and anti wind-up)
    :MV: MV (or Manipulated Value) vector
    :MV_P: MV_P (or Pr opotional part of MV) vector
    :MV_I: MV_I (or Integral part of MV) vector
    :MV_D: MV_D (or Derivative part of MV) vector
    :E: E (or control Error) vector
    :Man_FF: Activated FF in manual mode (optional: default boolean value is False)
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
        
        
    
            
    #Proportional-part
    MV_P.append(K_C*E[-1])
        
        
    #Itegral-part #slide 201
    if len(MV_I)<=2 : # peut etre changer la condition en len(E)<=2
        MV_I.append(( K_C * Ts/ T_I)* E[-1]) # attention veriefier je ne sais pas si c'est bien
    else : 
        MV_I.append(MV_I[-2]+(K_C * Ts/ T_I)*E[-1])
        
        
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

                  
            
            
    #windup slide 202
    if MV_P[-1]+MV_I[-1]+MV_D[-1]+MV_feed_forward>MV_MAX:
        MV_I[-1]=(MV_MAX-MV_P[-1]-MV_D[-1]-MV_feed_forward)
            
    elif MV_P[-1]+MV_I[-1]+MV_D[-1]+MV_feed_forward<MV_MIN :
        MV_I[-1]=(MV_MIN-MV_P[-1]-MV_D[-1]-MV_feed_forward)
                        
        
    if MAN[-1] == False :   #if we are not in manual mode  

        # MV_part
        MV_SUM = MV_P[-1] + MV_I[-1] + MV_D[-1] + MV_feed_forward
        
        
        #adding the value of MV in the MV vector (warning if the MV_SUM is< MV_MIN or >MV_MAX) 
        if (MV_SUM>=MV_MIN and MV_SUM<=MV_MAX):
            MV.append(MV_SUM)

        elif (MV_SUM<MV_MIN):
            MV.append(MV_MIN)

        elif (MV_SUM>MV_MAX) :
            MV.append(MV_MAX)

            
    
    else: #if we are in manual mode

        MV_I[-1]=MV_MAN[-1]- MV_P[-1] #slide 201 3eme ligne 

        #FeedForward in manual mode
        if MAN_FF == True : 
            MV_feed_forward = MV_FF[-1] 
        else:
            MV_feed_forward = 0
        
         # MV_part
        MV_SUM = MV_P[-1]+MV_I[-1] + MV_feed_forward #slide 201 4eme ligne
        
        #adding the value of MV in the MV vector (warning if the MV_SUM is< MV_MIN or >MV_MAX) 
        if (MV_SUM>=MV_MIN and MV_SUM<=MV_MAX):
            MV.append(MV_SUM)

        elif (MV_SUM<MV_MIN):
            MV.append(MV_MIN)

        elif (MV_SUM>MV_MAX) :
            MV.append(MV_MAX)
        
        
        
           
            
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
        
        # silde 186 Ligne I without T3 (T3=0)
        TauOLP = max(T1,T2)
        TauC = gamma*TauOLP
        Kc = (T1+T2)/(TauC+theta)
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



#def FF_RT(DV,Kd,Kp,Tlead1,Tlag1,Tlead2,Tlag2,thetad,thetap,Ts,FF_OUT,PV_LL1,PV_LL2,Dv0=0):
def FF_RT (DV , thetad, thetap, Ts, MVFFDelay, KFF, T1p, T1d, MVFFLL1, T2p, T2d, MV_FF, DV0=0):
          
    """
    The function "FF_RT" is feedforward realtime it's composed of one delay and two Lead-Lag
    DV: Disturbance vector
    thetad:
    thetap: 
    Ts: Sampling time
    MVFFDelay: MV Feedforward Delay vector
    KFF:rapport gain KFF = -Kd/Kp
    T1p: 
    T1d: 
    MVFFLL1: 
    T2p:
    T2d:
    MV_FF: MV_FF ( or Feedfortward) vector
    DV0=0:    
    """
    
    Delay_RT(DV - DV0 * np.ones_like(DV), np.max([thetad - thetap, 0]), Ts, MVFFDelay)
    LL_RT(MVFFDelay, KFF, T1p, T1d, Ts, MVFFLL1)
    LL_RT(MVFFLL1,1,T2p,T2d,Ts,MV_FF)
    
#---------------------------------------------------------------------  

class LS_Process:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kp'] = parameters['Kp'] if 'Kp' in parameters else 1.0
        self.parameters['theta'] = parameters['theta'] if 'theta' in parameters else 0.0
        self.parameters['T1'] = parameters['T1'] if 'T1' in parameters else 0.0
        self.parameters['T2'] = parameters['T2'] if 'T2' in parameters else 0.0    
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 0.0
        self.parameters['Kc/Ti'] = parameters['Kc/Ti'] if 'Ti' in parameters else 0.0
        self.parameters['KcTd/TfdS+1'] = parameters['KcTd/TfdS+1'] if 'Td' in parameters else 0.0
        
#---------------------------------------------------------------------  

class Controller:
    
    def __init__(self, parameters):
        
        self.parameters = parameters   
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 1.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 0.0
        self.parameters['Td'] = parameters['Td'] if 'Td' in parameters else 0.0
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 0.0


#---------------------------------------------------------------------  
    

def Margins(P,C,omega, Show = True):
    
    """
    
    L(s)=(Kp*Kc)/((T1*s+1)*(T2*s+1))*(1+1/(Ti*s)+(Td*s)/(Tfd*s+1))
        Ptheta=e^(-theta*s)
        PGain= Kp*Kc
        P1=1/(T1*s+1)
        P2=1/(T2*s+1)
        PID=(1+1/(Ti*s)+(Td*s)/(Tfd*s+1))
    
    :P: Process as defined by the class "Process".
        Use the following command to define the default process which is simply a unit gain process:
            P = Process({})
        
        Use the following commands for a SOPDT process:
            P.parameters['Kp'] = Kp
            P.parameters['Tlag1'] = T1_p
            P.parameters['Tlag2'] = T2_p
            P.parameters['theta'] = theta_p

    :C: Controller as defined by the class "Controler".
        Use the following command to define the default process which is simply a PID controller:
            C = Controller({})
        
        Use the following commands for the PID process:
            C.parameters['Kc'] = Kc
            C.parameters['Ti'] = Ti
            C.parameters['Td'] = Td   
            C.parameters['alpha']  = alpha 
        
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
    :Show: boolean value (optional: default value = True). If Show = True, the Bode diagram is shown. Otherwise Ls (P(j omega)) (vector of complex numbers) is returned.
    
    The function "LSBode" generates the Bode diagram of the process P
    """     
    
    s = 1j*omega
    
    Ptheta = np.exp(-P.parameters['theta']*s)
    PGain = (C.parameters['Kc']*P.parameters['Kp'])
    P1 = (1/((P.parameters['Tlag1']*s + 1)*(P.parameters['Tlag2']*s + 1)))
    #P2 = (1/(P.parameters['Tlag2']*s + 1))
    PID = (1+(C.parameters['Td']*s/(C.parameters['Td']*C.parameters['alpha']*s + 1)) + (1/(C.parameters['Ti']*s)))
    
    
    Ls = np.multiply(Ptheta,PGain)
    Ls = np.multiply(Ls,P1)
    #Ls = np.multiply(Ls,P2)
    Ls = np.multiply(Ls,PID)


    # Find the gain and phase margins

    gain = 20*np.log10(Ls).real
    phase = (180/np.pi)*np.unwrap(np.angle(Ls))
    
    indiceGain= Find_value(gain, 0)
    indicePhase=Find_value(phase, -180)
    

    print ('indiceGain: ', indiceGain)
    print ('indicePhase: ', indicePhase)

    #frequencyGain = omega[indicePhase-1]
    #frequencyPhase = omega[indiceGain-1]

    #gainMargin = 20*np.log10(np.abs(Ls[indiceGain]))
    #phaseMargin = 180 + (180/np.pi)*np.angle(Ls[indicePhase])
        

    
    fig, (ax_gain, ax_phase) = plt.subplots(2,1)
    fig.set_figheight(12)
    fig.set_figwidth(22)
    
    # Gain part
    ax_gain.semilogx(omega,20*np.log10(np.abs(Ls)),label='L(s)')
    gain_min = np.min(20*np.log10(np.abs(Ls)/5))
    gain_max = np.max(20*np.log10(np.abs(Ls)*5))
    ax_gain.set_xlim([np.min(omega), np.max(omega)])
    ax_gain.set_ylim([gain_min, gain_max])
    ax_gain.set_ylabel('Amplitude |L(s)| [db]')
    ax_gain.set_title('Bode plot of L(s) alpha: '+str(C.parameters['alpha']))
    ax_gain.legend(loc='best')

    # Phase part
    
    ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ls)),label='L(s)')   
    ax_phase.set_xlim([np.min(omega), np.max(omega)])
    ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ls))) - 10
    ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ls))) + 10
    ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
    ax_phase.set_ylabel(r'Phase $\angle L(s)$ [°]')
    ax_phase.legend(loc='best')


    # Draw the gain and phase crossover frequencies
    ax_gain.axhline(y=0, color='g', linestyle='--')
    ax_phase.axhline(y=-180, color='r', linestyle='--')

    
    try: 
        frequencyGain = omega[indicePhase-1]
        ax_gain.axvline(x=frequencyGain, color='y', linestyle='-')
        ax_phase.axvline(x=frequencyGain, color='y', linestyle='-')
        print ('Gain margin : ', gain[indicePhase-1] ,'dB at ', frequencyGain, 'rad/s ')
        
        
    except:
        print ('Error in the frequency gain computation')
        
    try :
        frequencyPhase = omega[indiceGain-1]        
        ax_gain.axvline(x=frequencyPhase, color='b', linestyle='-')
        ax_phase.axvline(x=frequencyPhase, color='b', linestyle='-')
        print ('Phase margin :',(180-phase[indiceGain-1]),'° at ',frequencyPhase,'rad/s' )
        
    except:
        print('Error in the frequency pahse computation ')

     
    
    

#---------------------------------------------------------------------  

def Find_value(vector, discriminat_Value):
    """
    vector: a vector
    discriminat_Value : discriminat Value
    return the indice below a discriminant value
    """
    
    
    for i in range(len(vector)):
        if vector[i] < discriminat_Value:
            return i
    
#---------------------------------------------------------------------  



    
    
    