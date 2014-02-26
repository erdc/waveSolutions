import sys, math, random
from numpy import *

def sigma(omega, omega_peak):
    """ This function returns value of sigma in JONWSAP.
    Note that the inputs are acutally omega**2 and omega_peak**2,
    but it doesn't matter for the return value!

    This module (spectrum.py) is imported by the main JONSWAP.py module.
    It defines a JONSWAP spectrum and adds a desired degree of directional
    spreading with a random phase. 

    AUTHOR:  Matt Malej
    EMAIL:   matt.malej@erdc.dren.mil"""

    if omega <= omega_peak:
        sigma = 0.07
    else:
        sigma = 0.09
    # end if

    return sigma

def factorial_1(ns):
    fact = 1        # ~ default
    for factor in range(int(ns),0,-1):
        fact = fact * factor
    return fact

def factorial_2(ns):
    fact = 1        # ~ default
    for factor in range(int(ns),0,-2):
        fact = fact * factor
    return fact

def spreading(vkx, vky, ns, max_angle):
    temp1 = factorial_1(ns)      # user defined std. factorial function 
    temp2 = factorial_2(ns)      # user defined augmented factorial function
    
    # Normalization Constant
    D0 = temp1 * 2.0**(int(ns)-1) / (max_angle * temp2)
        
    if vkx == 0.0:
        spread = 0.0
    else:
        angle = arctan(vky/vkx)        
        if abs(angle) <= max_angle:
            spread = D0*( cos(pi/2.0*angle/max_angle)**(2*ns) )
        else:
            spread = 0.0

    return spread


def jonswap(nx, ny, Hs, gamma, fp, nspread, thetam, gv, kx, ky):
    # Constansts
    alpha = 5.0 / 16.0
    beta = 5.0 / 4.0        # *** different from Beta and theta_m !!!
    
    kx_min = kx[1]
    ky_min = ky[1]
    
    spectrum = zeros((nx,ny),complex)
    
    omega_peak = 2.0*pi*fp
    kp = omega_peak**2 / gv      # peak wavenumber (determinied from disper. rel. with h->infty)
    
    spread = spreading(kx[10], ky[10], nspread, thetam)
    
    # Define the spectrum over HALF the modes as the input is REAL and thus
    # the other half of the modes is just the complex conjugate of the first!
    for j in range(0, ny, 1):
        for i in range(0, nx/2+1, 1): # ~ loops over nx/2+1 (from  0^th to nx/2^th) terms
            kk = sqrt(kx[i]**2 + ky[j]**2)
            if kk == 0.0:
                spectrum[i,j] = 0.0 + 0.0j
            else:
                # sigma factor
                temp_sigma = sigma(gv*kk, gv*kp)  # ~ sending omegas SQUARED !
                             
                # spreading factor
                temp_spread = spreading(kx[i], ky[j], nspread, thetam)
                                
                # computing random phase (seed is apparently the system clock)
                phase = random.random()     #  where:  0 <= phase <= 1
                
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #    <<<<<<<<<< Computing JONSWAP Function >>>>>>>>>>
                #
                # Recall: the spectrum is composed of F(kk,theta)=S(kk)*D(theta)
                #         where S(kk)=S(kx,ky) is the JONSWAP spectrum and
                #         D(theta) is the directional spreading function!
                #
                temp1 = alpha * kp**2 * Hs**2 / (2.0 * kk**4)
                temp2 = exp(-beta * (kp/kk)**2)
                temp3 = gamma**( exp(- ( (sqrt(kk)-sqrt(kp))**2 / (2.0*(temp_sigma**2)*kp) ) ) )
                temp_JONSWAP = temp1 * temp2 * temp3
                
                # spectrum[i,j] ==> F(kk,theta) = S(kx,ky) * D(theta)
                spectrum[i,j] = sqrt(2.0 * temp_JONSWAP * temp_spread * kx_min * ky_min) * exp(2.0j * pi * phase)
                #  NOTE: we added sqrt(2.0) into the spectrum ******* !!!!!!
                #        ... if unclear  as to how the line above comes about
                #            see attached PDF for theory on relation between
                #            the above spectrum[i,j] and autocorrelation of the
                #            surface elevation (Weiner-Khinchin Theorem)
 
    # For real transforms the highest and lowest modes  are real!!!
    spectrum[0,0] = real(spectrum[0,0])
    spectrum[nx/2,0] = real(spectrum[nx/2,0])
    spectrum[0, ny/2] = real(spectrum[0, ny/2])
    spectrum[nx/2, ny/2] = real(spectrum[nx/2, ny/2])
        
 
    # Constructing the other half (complex conjugates) ~ linear vectors first
    spectrum[nx/2+1:nx-1, 0] = spectrum[nx/2-1:1:-1, 0].conj()          # along ky = 0
    spectrum[0, ny/2+1:ny-1] = spectrum[0, ny/2-1:1:-1].conj()          # along kx = 0
    spectrum[nx/2, ny/2+1:ny-1] = spectrum[nx/2, ny/2-1:1:-1].conj()    # along kx = Nx/2
    spectrum[nx/2+1:nx-1, ny/2] = spectrum[nx/2-1:1:-1, ny/2].conj()    # along ky = Ny/2
    
    # Now for the complex conjugates of the quadrants!    
    for j in range (1, ny/2, 1):
        for i in range (nx/2+1, nx, 1):            
            spectrum[i,j] = spectrum[nx-i, ny-j].conj()
    
    for j in range (ny/2+1, ny, 1):
        for i in range (nx/2+1, nx, 1):
            spectrum[i,j] = spectrum[nx-i, ny-j].conj()
     
    return spectrum
