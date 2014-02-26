import numpy as np

def velocityPotential(spec, gv, kx, ky):
    """NOTE: 
This module defines the velocity potential based on linear 
theory from the free-surface spectrum. Higher order corrections
can be implemented if needed, but recall that the linear solution
for velocity potential of a Stokes wave is valid up to third
in wave steepness, hence capturing the required bound and free waves. """    
    [nx, ny] = spec.shape
    potentialFourier = np.zeros((nx,ny), complex)      # ~ for fft(vel potential)
       
    # ~ Define the spectrum over HALF the modes as the input is REAL and thus
    # the other half of the modes is just the complex conjugate of the first!
    for j in range(0, ny, 1):
        for i in range(0, nx/2+1, 1): # ~ loops over nx/2+1 (from  0^th to nx/2^th) terms
            kk = np.sqrt(kx[i]**2 + ky[j]**2)
            omega = np.sqrt(gv*kk)
            if omega == 0.0:
                potentialFourier[i,j] = 0.0
            else:
                potentialFourier[i,j] = -1.0j * gv * spec[i,j] / omega
                # Recall: spec = fft2(free surface)
            
    # Recall: for real transforms the highest and lowest modes are real!!!
    potentialFourier[0,0] = np.real(potentialFourier[0,0])
    potentialFourier[nx/2,0] = np.real(potentialFourier[nx/2,0])
    potentialFourier[0, ny/2] = np.real(potentialFourier[0, ny/2])
    potentialFourier[nx/2, ny/2] = np.real(potentialFourier[nx/2, ny/2])
        
    # Constructing the other half (complex conjugates) ~ linear vectors first
    potentialFourier[nx/2+1:nx-1, 0] = potentialFourier[nx/2-1:1:-1, 0].conj()        # along ky = 0
    potentialFourier[0, ny/2+1:ny-1] = potentialFourier[0, ny/2-1:1:-1].conj()        # along kx = 0
    potentialFourier[nx/2, ny/2+1:ny-1] = potentialFourier[nx/2, ny/2-1:1:-1].conj()  # along kx = Nx/2
    potentialFourier[nx/2+1:nx-1, ny/2] = potentialFourier[nx/2-1:1:-1, ny/2].conj()  # along ky = Ny/2
    
    # Now for the complex conjugates of the quadrants!    
    for j in range (1, ny/2, 1):
        for i in range (nx/2+1, nx, 1):            
            potentialFourier[i,j] = potentialFourier[nx-i, ny-j].conj()
    
    for j in range (ny/2+1, ny, 1):
        for i in range (nx/2+1, nx, 1):
            potentialFourier[i,j] = potentialFourier[nx-i, ny-j].conj()

    # Create a 2D meshgrid of wavenumbers
    [kxx, kyy] = np.meshgrid(kx,ky)

    # Horizontal velocityties (term-wise multiplication)
    uVelocity = 1.0j*kxx*potentialFourier
    vVelocity = 1.0j*kyy*potentialFourier
    
    # ~ Inverting back to physical space --> take REAL part just in case we have
    #   some left over (nonzero) imaginary part of 'potentialFourier' and u,v
    potential = np.real( np.fft.ifft2(potentialFourier) )
    u = np.real( np.fft.ifft2(uVelocity) )
    v = np.real( np.fft.ifft2(vVelocity) )
    
    return potential, u, v

