
from math import pi
import numpy as np          # modules for computing FFT and other numerical stuff
import JONSWAP_p as JS
from spectrum import jonswap
from potential import velocityPotential

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import mpl_toolkits.mplot3d.axes3d as p3
except:
    pass


def JONSWAP():
    """Sets a wave field according to JONSWAP ocean wave spectrum.
A script for initializing data files of a Quasi-3D wave simulation via
the JONSWAP spectrum and directional distribution.

GENERAL NOTES
=============

MODULE:  JONSWAP.py (Python Script File)

AUTHOR:   Matt Malej (matt.malej@erdc.dren.mil)

PURPOSE:  This Python script is designed to set up and initalize the wave
          field via a JONSWAP spectrum and for both directionally confined and
          and directional varied distribution.
          Will serve as input for simulating linear and eventually weakly nonlinear
          random surface water waves in realistic settings. 

          The module imports JONSWAP_p.py (p = physics), where all the input
          parameters are specified. 

INPUT:    Phyical parameters from JONSWAP_p.py

Self NOTE: for FFT your transforming function needs to be periodic on your 
           domain and do NOT make your first and final point in the domain equal!
           ... i.e. for sin(x) do NOT make x=[0,2*pi], but instead x=[dx,2*pi]!

LAST UPDATE: January 14, 2014
"""

    # Some Global Variables
    omega_peak = 2.0*pi*JS.fp              # peak angular frequency
    kp = (2.0*pi*JS.fp)**2 / JS.gv      # peak wavenumber (determinied from disper. rel. with h-->infty)
    neqs = 2                               # number of equations

    # Discretization
    Nx = 2**JS.npw1                # number of Fourier modes
    Ny = 2**JS.npw2
    Lx = abs(JS.x1e - JS.x1b)         # domain length
    Ly = abs(JS.y1e - JS.y1b)
    kc_x = 2*pi/Lx
    kc_y = 2*pi/Ly              # std. wave factors if modes would be integers
    kc_x_modes = Nx*kc_x        
    kc_y_modes = Ny*kc_y        # wave factors for modes obtained from fftfreq()

    modes_x = np.fft.fftfreq(Nx)    # Fourier mode alignments
    modes_y = np.fft.fftfreq(Ny)

    kx = kc_x_modes * modes_x
    ky = kc_y_modes * modes_y

    # Generating mesh grid for 2D spectrum
    [kxx, kyy] = np.meshgrid(kx,ky)

    # ~ Call jonswap() in spectrum.py module *** Send 1D vectors kx and ky (not kxx nor kyy)!
    spectrum = jonswap(Nx, Ny, JS.Hs, JS.gamma, JS.fp, JS.nspread, JS.thetam, JS.gv, kx, ky)

    # (DISREGARD FOR NOW) - Imposing iniitial circular zero-pad filter
    if JS.filter == 1:
        spectrum[ 2*floor(kp/kc_x)+1 : (Nx-1)-2*floor(kp/kc_x), 2*floor(kp/kc_x)+1 : (Ny-1)-2*floor(kp/kc_x) ] = 0.0 + 0.0j 
    else:
        pass

    # *** Compute the surface elvation via iFFT2D of the spectrum  --> take REAL part
    #     just in case we have some left over (nonzero) imaginary part of 'spectrum'
    surface = np.real(np.fft.ifft2(spectrum) )

    # *** Compute the velocity potential from linear theory and horizontal velocity
    velPotential, velocity_u, velocity_v = velocityPotential(spectrum, JS.gv, kx, ky)

    try:
        #fig1 = plt.figure()
        #plt.clf()
        #plt.plot(kx[0:Nx/2-1,], abs(spectrum[0:Nx/2-1,0]), 'b', linewidth=2)
        #plt.show()

        # Attaching 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)
        
        x = np.linspace(JS.x1b, JS.x1e, Nx)
        y = np.linspace(JS.y1b,JS.y1e, Ny)
        [xx, yy] = np.meshgrid(x, y)
        
        #surf = ax.plot_surface(xx,yy,surface,rstride=2,cstride=2, cmap=cm.jet,
        #    linewidth=0.5, antialiased=False)
        ax.plot_wireframe(xx,yy,surface, rstride=4, cstride=4)
        plt.show()
    except:
        pass

    # Returning surface elvation along the line of wavemakers [x,y]
    loc = 10
    return surface #[loc,:]


if __name__=='__main__':
    JONSWAP()
