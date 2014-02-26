import numpy as np
from math import *

# Water                                                                                                                                                                 
rho_0 = 998.2
nu_0  = 1.004e-6

# Air                                                                                                                                                                   
rho_1 = 1.205
nu_1  = 1.500e-5

# Surface tension                                                                                                                                                        
sigma_01 = 0.0

# Gravity                                                                                                                                                                
g = [0.0,0.0,-9.8]


# Domain
depthFactor=8.0#16.0  # ...TODO: remove depthFactor after debugging                                                                                                   
L = (20.0,
     0.25,
     depthFactor*0.61)

#wave/current properties                                                                                                                                                 
windspeed_u = 0.0
windspeed_v = 0.0
windspeed_w = 0.0

outflowHeight = 0.5*L[2]
outflowVelocity = (0.0,0.0,0.0)#not used for now                                                                                                                           

inflowHeightMean = 0.5*L[2]
inflowVelocityMean = (0.0,0.0,0.0)

regime = 30.0 # if regime > 25 ==> shallow water, < 4 ==> deep water, between ==> finite depth                                                                           
waveLength = regime*inflowHeightMean/depthFactor # xSponge                                                                                                                
k=(2.0*pi/waveLength,0.0,0.0) # if k[1]~0 ==> long waves (unidirectional) 

# NOTE: For Shallow Water Limit:  h < waveLength/25 ==> omega ~ sqrt(g*k^2*h) ~ 2*pi/period (no dispersion)                                                                #       For Deep Water Limit:     h > waveLength/4  ==> omega ~ sqrt(g*k)                                                                                                  #       For Finite Depth: waveLength/25 < h < waveLength/4 ==> omega = sqrt(g*k*tanh(k*h))                                                                               
if inflowHeightMean < (waveLength*25.0):
    omega = np.sqrt(-g[2]*inflowHeightMean)*k[0]
    df_dk = -g[2]*2.0*k[0]*inflowHeightMean
elif inflowHeightMean > (waveLength*4.0):
    omega = np.sqrt(-g[2]*k[0])
    df_dk = -g[2]
else:
    omega = np.sqrt(-g[2]*k[0]*np.tanh(k[0]*inflowHeightMean))
    df_dk = -g[2]*( np.tanh(k[0]*inflowHeightMean) + k[0]*inflowHeightMean/np.cosh(k[0]*inflowHeightMean)**2 )

#####################################################
# Setting desired level on nonlinearity for testing:                                                                                                                                  
# ... epsilon ~ 0.1 ==> weakly nonlinear, epsilon ~ 0.5 ==> highly nonlinear                                                                                              
epsilon = 0.05 # 0.01,0.02,0.05,0.1,0.15,0.2,0.4 # wave steepness                                                           
factor = epsilon*regime/(2*np.pi*depthFactor) # factor == amplitude/depth                                                                                                            
amplitude = inflowHeightMean*factor
period = 2.0*pi/omega

# Group Velocity ==> d/dk{omega} = f'(k)/(2*omega), where f'(k)=d/dk{omega(k)^2}                                                                                          
groupVelocity = df_dk / (2.0*omega)

# Add random phase                                                                                                                                                        
randomPhase = False

# Debugging 
print "regime is: ", regime
print "depthFactor is: ", depthFactor
print "wavelength is: ", waveLength
print "amplitude is: ", amplitude
