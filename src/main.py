import os
import sys

import numpy as np
import defaultWaveModules as wm
  
# Domain and mesh
depthFactor=4.0#16.0  # ...TODO: remove depthFactor after debugging
L = (10.0,
     0.0,#0.25,
     depthFactor*0.61)

# Water
rho_0 = 998.2
nu_0  = 1.004e-6

# Air
rho_1 = 1.205
nu_1  = 1.500e-5 

# Gravity
g = [0.0,0.0,-9.8]

inflowHeightMean = 0.5*L[2]
inflowVelocityMean = (0.0,0.0,0.0)

regime = 30.0 # if regime > 25 ==> shallow water, < 4 ==> deep water, between ==> finite depth
waveLength = regime*inflowHeightMean/depthFactor # xSponge
k=(2.0*pi/waveLength,0.0,0.0)
# NOTE: For Shallow Water Limit:  h < waveLength/25 ==> omega ~ sqrt(g*k^2*h) ~ 2*pi/period (no dispersion)
#       For Deep Water Limit:     h > waveLength/4  ==> omega ~ sqrt(g*k)
#       For Finite Depth: waveLength/25 < h < waveLength/4 ==> omega = sqrt(g*k*tanh(k*h))
if inflowHeightMean < (waveLength*25.0):
    omega = np.sqrt(-g[2]*inflowHeightMean)*k[0]
    df_dk = -g[2]*2.0*k[0]*inflowHeightMean
elif inflowHeightMean > (waveLength*4.0):
    omega = np.sqrt(-g[2]*k[0])
    df_dk = -g[2]
else:
    omega = np.sqrt(-g[2]*k[0]*np.tanh(k[0]*inflowHeightMean))
    df_dk = -g[2]*( np.tanh(k[0]*inflowHeightMean) + k[0]*inflowHeightMean/np.cosh(k[0]*inflowHeightMean)**2 )

# Setting desired level on nonlinearity: 
# ... epsilon ~ 0.1 ==> weakly nonlinear, epsilon ~ 0.5 ==> highly nonlinear
epsilon = 0.1 # 0.01,0.02,0.05,0.1,0.15,0.2,0.4 # wave steepness
factor = epsilon*regime/(2*np.pi*depthFactor) # factor == amplitude/depth 
amplitude = inflowHeightMean*factor
period = 2.0*pi/omega

# Group Velocity ==> d/dk{omega} = f'(k)/(2*omega), where f'(k)=d/dk{omega(k)^2}
groupVelocity = df_dk / (2.0*omega)

# Add random phase
randomPhase = False

################
# DEBUGGING
print "WAVE STEEPNESS IS: ", epsilon
print "AMPLITUDE IS: ", amplitude
print "GROUP SPPED IS: ", groupVelocity
print "WAVELENGTH IS: ", waveLength
#sys.exit()
##############

# Wave Field Object
waveField = wm.Linear2D(amplitude,omega,k,inflowHeightMean,rho_0,rho_1,randomPhase)
#waveField = wm.Linear3D(amplitude,omega,k,inflowHeightMean,rho_0,rho_1,randomPhase)
#waveField = wm.WaveGroup(amplitude,omega,k,inflowHeightMean,rho_0,rho_1,randomPhase)
#waveField = wm.Solitary(amplitude,omega,k,inflowHeightMean,rho_0,rho_1)
#waveField = wm.StokesWave(amplitude,omega,k,inflowHeightMean,rho_0,rho_1)
#waveField = wm.JONSWAP(amplitude,omega,k,inflowHeightMean,rho_0,rho_1)

#c_soliton = sqrt(fabs(g[2])*(inflowHeightMean+amplitude))

#cek debugging
#windVelocity = (0.0,0.0,0.0)
#inflowHeightMean = 2.64
#inflowVelocityMean = (0.0,0.0,0.0)
#period = 2.0
#omega = 2.0*math.pi/period
#waveheight = 0.4
#amplitude = waveheight/ 2.0
#wavelength = 6.2
#k = 2.0*math.pi/wavelength

#def theta(x,t):
#    return k[0]*x[0] - omega*t
#def z(x):
#    return x[2] - inflowHeightMean
#sigma = omega - k[0]*inflowVelocityMean[0]
#h = inflowHeightMean# - transect[0][1]
#print "h",h,"sigma",sigma,"inflowHeightMen",inflowHeightMean,"omega",omega,"k",k
#def waveHeight(x,t):
#    return inflowHeightMean + amplitude*cos(theta(x,t))
#def waveVelocity_u(x,t):
#    return sigma*amplitude*cosh(k[0]*(z(x)+h))*cos(theta(x,t))/sinh(k[0]*h)
#def waveVelocity_v(x,t):
#    return 0.0
#def waveVelocity_w(x,t):
#    return sigma*amplitude*sinh(k[0]*(z(x)+h))*sin(theta(x,t))/sinh(k[0]*h)
#cek debugging

def waveHeight(x,t):
     return inflowHeightMean + waveField.height(x,t)
#     # --------------------- CEK original -----------------------
#     #T = min(t,100) - 4.0
#     #return inflowHeightMean + amplitude/cosh(sqrt(3.0*amplitude/(4.0*inflowHeightMean**3)) * (x[0] - c_soliton*T))**2
#     #return inflowHeightMean + amplitude*sin(omega*t-k[0]*x[0])
#     # ----------------------------------------------------------

def waveVelocity_u(x,t):
     z = x[2] - inflowHeightMean
     return inflowVelocityMean[0] + waveField.velocity_u(x,t)
#     # CEK: return inflowVelocityMean[0] + omega*amplitude*cosh(k[0]*(z + inflowHeightMean))*sin(omega*t - k[0]*x[0])/sinh(k[0]*inflowHeightMean)
#     # CEK: return c_soliton*(waveHeight(x,t)-inflowHeightMean)/waveHeight(x,t)

def waveVelocity_v(x,t):
     z = x[2] - inflowHeightMean
     return inflowVelocityMean[1] + waveField.velocity_v(x,t)
#     # CEK:
#     #return 0.0

def waveVelocity_w(x,t):
     z = x[2] - inflowHeightMean
     return inflowVelocityMean[2] + waveField.velocity_w(x,t)
#     # CEK: return inflowVelocityMean[2] + omega*amplitude*sinh(k[0]*(z + inflowHeightMean))*cos(omega*t - k[0]*x[0])/sinh(k[0]*inflowHeightMean)
####

def wavePhi(x,t):
    return x[2] - waveHeight(x,t)

def wavePhi_init(x,t):
    return x[2] - inflowHeightMean # mean/flat initial surface profile
    # CEK original: return wavePhi(x,t) # interface is initialized at t=0 (not flat) 

# Computation Time for Wave(s) to return to wave maker (based on groupVelocity)
# ...TODO: remove debugFactor when done debugging 
T=2*period#3.00
print "Total Time of Computation is: ",T
