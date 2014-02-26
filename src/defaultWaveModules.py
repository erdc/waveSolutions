import numpy as np
from math import pi
import random
import JONSWAP_p as JS
""" Set of wave modules for initializing and
    driving the wavemaker field, that will be
    fed in and constrasted with wave flume 
    experimental data."""

class Linear2D:
    """
    A class for linearized solutions of 2D flow (1D + surface) for
    travelling progressive waves ~ exp(i(kx-wt))
    
    .. todo:: Finish the docs

    An equation
    
    .. math:: 
        
        \Phi_t = -g \zeta + ...

    More text, inline math :math:`x^3` 
    """
    #see this url for some useful restructured text directives
    # ==> http://matplotlib.sourceforge.net/devel/documenting_mpl.html#formatting-mpl-docs
    def __init__(self,amplitude,omega,k,depth,rho_0,rho_1,randomPhase):
        self.name = 'Linear2D'
        self.A = amplitude
        self.omega = omega
        self.k = k
        self.h = depth
        self.rho_0 = rho_0      # density of water
        self.rho_1 = rho_1      # density of air
        self.randomPhase = randomPhase
        if randomPhase:
                # computing random phase (random seed is the system clock)
                self.phase = 2*np.pi*random.random()   # where: 0< random.random()<= 1
        else:
                self.phase = 0.0    
    def theta(self,x,t):
        return self.k[0]*x[0] - self.omega*t + self.phase
    
    def height(self,x,t):
        """ Gives a linearized solution for the air-water interface to the 
            potential flow model in two dimensions (x,y,z=eta) for finite depth."""
        eta = self.A*np.exp(1j*self.theta(x,t))
        # NOTE: x[0] is normally exptected to be a vector here (in waveTank problem)!
        return np.real(eta)
        
    def velocity_u(self,x,t):
        """ Gives a linearized solution velocity in x-dir to the potential flow
            model in two dimensions (x,y,z) for finite depth, as well as, deep
            and shallow water limits.

            .. todo:: implement deep & shallow water limits."""
        g = (0.0,0.0,-9.81)         # gravity
        z = x[2] - self.h           # mean height would now =0
        
        # Finite Depth (0 < kh < infty)
        u = (g[2]*self.k[0]*self.A / self.omega) * np.cosh(self.k[0]*(z+self.h))/np.cosh(self.k[0]*self.h) * np.exp(1j*self.theta(x,t))  

        # Deep water (kh >> 1)
        # u = (g[2]*self.k[0]*self.A / self.omega) * np.exp(k[0]*z) * \
        #    np.exp(1j*(self.k[0]*x[0] - self.omega*t + self.phase))        
        # Shallow water (kh << 1)
        # u = (g[2]*self.k[0]*self.A / self.omega) * np.exp(1j*(self.k[0]*x[0] - self.omega*t + self.phase)) 
        return np.real(u)
 
    def velocity_v(self,x,t):
        v = 0.0 # 2D velocity only (x,z)
        return v

    def velocity_w(self,x,t):
        """ Gives a linearized solution velocity in x-dir to the potential flow
            model in two dimensions (x,y,z) for finite depth, as well as, deep
            and shallow water limits.

            .. todo:: implement deep & shallow water limits."""
        g = (0.0,0.0,-9.81)                          # gravity
        z = x[2] - self.h

        # Finite Depth (0 < kh < infty)                                 
        w = -1j*(g[2]*self.k[0]*self.A/self.omega) * np.sinh(self.k[0]*(z+self.h))/np.cosh(self.k[0]*self.h) * np.exp(1j*self.theta(x,t))
        # Deep Water (kh >> 1)
        # w = TODO
        # Shallow Water
        # ... TODO
        return np.real(w)
        
    def pressure(self,x,t):
        """ Gives linearized pressured with P_atm = 0 """
        g = (0.0,0.0,-9.81)
        z = x[2] - self.h
        p = self.rho_0*g[2]*self.A* np.cosh(self.k[0]*(z+self.h))/np.cosh(self.k[0]*self.h) \
            * np.exp(1j*self.theta(x,t))
        return np.real(p)


class Linear3D:
    """
    A class for linearized solutions of 2D flow (1D + surface) for
    travelling progressive waves ~ exp(i(kx+ky-wt))
    
    .. todo:: Finish the docs

    An equation
    
    .. math:: 
        
        \Phi_t = -g \zeta + ...

    More text, inline math :math:`x^3` 
    """
    #see this url for some useful restructured text directives
    #
    #http://matplotlib.sourceforge.net/devel/documenting_mpl.html#formatting-mpl-docs
    def __init__(self,amplitude,omega,k,depth,rho_0,rho_1,randomPhase):
        self.name = 'Linear3D'
        self.A = amplitude
        self.omega = omega
        self.k = k
        self.kappa = np.sqrt(k[0]**2+k[1]**2)
        self.h = depth
        self.rho_0 = rho_0      # density of water
        self.rho_1 = rho_1      # density of air
        self.randomPhase = randomPhase
        if randomPhase:
                # computing random phase (random seed is the system clock)
                self.phase = 2*np.pi*random.random()     #  where:  0 <= random.random() <= 1
        else:
                self.phase = 0.0    
       
    def theta(self,x,t):
        return (self.k[0]*x[0] + self.k[1]*x[1] - self.omega*t + self.phase)
    

    def height(self,x,t):
        """ Gives a linearized solution for the air-water interface to the 
            potential flow model in two dimensions (x,y,z=eta) for finite depth."""
        eta = self.A*np.exp(1j*self.theta(x,t))
        # NOTE: x[0] is normally exptected to be a vector here (in waveTank problem)!
        return np.real(eta)

    def velocity_u(self,x,t):
        """ Defines a linearized solution to the potential flow
            model in two dimensions (x,y,z) for finite depth,
            as well as, deep and shllow water limits.

            .. todo:: implement deep & shallow water limits."""
        g = (0.0,0.0,-9.81)         # gravity
        z = x[2] - self.h           # mean height would now =0

        # Finite Depth (0 < kh < infty)
        u = (g[2]*self.k[0]*self.A/self.omega) * np.cosh(self.kappa*(z+self.h))/np.cosh(self.kappa*self.h) * np.exp(1j*self.theta(x,t))  
        # Deep water (kh >> 1)
        # ... TODO                 
        # Shallow water (kh << 1)
        # ... TODO
        return np.real(u)

    def velocity_v(self,x,t):
        """ Gives a linearized solution velocity in y-dir to the potential flow
            model in two dimensions (x,y,z) for finite depth, as well as, deep
            and shallow water limits.

            .. todo:: implement deep & shallow water limits."""
        v = (g[2]*self.k[1]*self.A/self.omega) * np.cosh(self.kappa*(z+self.h))/np.cosh(self.kappa*self.h) * np.exp(1j*self.theta(x,t))
        return np.real(v)

    def velocity_w(self,x,t):
        """ Gives a linearized solution velocity in z-dir to the potential flow
            model in two dimensions (x,y,z) for finite depth, as well as, deep
            and shallow water limits.

            .. todo:: implement deep & shallow water limits."""
        g = (0.0,0.0,-9.81)                          # gravity
        z = x[2] - self.h

        # Finite Depth (0 < kh < infty)                                 
        w = -1j*(g[2]*self.kappa*self.A/self.omega) * np.sinh(self.kappa*(z+self.h))/np.cosh(self.kappa*self.h) * np.exp(1j*theta(x,t))
        # Deep Water (kh >> 1)
        # ... TODO
        # Shallow Water
        # ... TODO
        return np.real(w)

    def pressure(self,x,t):
        """ Gives linearized pressured with P_atm = 0 """
        g = (0.0,0.0,-9.81)
        z = x[2] - self.h        
        p = self.rho_0*(-g[2])*self.A* np.cosh(self.kappa*(z+self.h))/np.cosh(self.kappa*self.h) * np.exp(1j*self.theta(x,t))
        return np.real(p)


class WaveGroup:
    """ Class that defines a nearly monochromatic
        wave train/group of waves of same amplitude.
        
        .. todo:: Finish the docs. """

    def __init__(self,amplitude,omega,k,depth,rho_0,rho_1,randomPhase):
        self.name = 'WaveGroup'
        self.A = amplitude
        self.omega = omega
        self.k = k
        self.h = depth
        self.rho_0 = rho_0      # density of water                                      
        self.rho_1 = rho_1      # density of air
        self.N = 1              # number of nearly-monochromatic pairs
        self.diff = 0.05
        self.g = (0.0,0.0,-9.81)                          # gravity
        self.randomPhase = randomPhase
        if randomPhase:
                # computing random phase (random seed is the system clock)
                self.phase = 2*np.pi*random.random()     #  where:  0 <= random.random() <= 1
        else:
                self.phase = 0.0   

    def height(self,x,t):
        theta =  self.k[0]*x[0] - self.omega*t + self.phase# ~ NOTE: x[0] is a vector here!
        dtheta = self.diff*theta
        eta = np.zeros(x[0].shape) #self.A*np.sin(theta)

        for i in range(self.N):
            eta = eta + self.A*np.exp(1j*(i+1)*(theta+(i+1)*dtheta)) + self.A*np.exp(1j*(i+1)*(theta-(i+1)*dtheta))
            # NOTE: variations on the order of dtheta^2 are truncated (only up to order of group velocity kept)
            #       and inflowHeightMean already added in wavetank.py
        
        return np.real(eta)


    def velocity_u(self,x,t):
        """ Defines a linearized solution to the potential flow
            model in two dimensions (x,y,z) for finite depth,
            as well as, deep and shllow water limits, for slowly
            varying regular wavetrains.

            .. todo:: implement deep & shallow water limits. """
        theta =  self.k[0]*x[0] - self.omega*t + self.phase # ~ NOTE: x[0] is a vector here!
        dtheta = self.diff*theta
        z = x[2] - self.h
        u = np.zeros(x[0].shape)

        # Finite Depth (0 < kh < infty)
        for i in range(self.N):
            diffPos = (1+self.diff*(i+1))
            diffNeg = (1-self.diff*(i+1))
            u = u + (-self.g[2]*diffPos*self.k[0]*self.A / (diffPos*self.omega)) * \
                np.cosh(diffPos*self.k[0]*(z+self.h))/np.cosh(diffPos*self.k[0]*self.h) * \
                np.exp(1j*(i+1)*(theta+(i+1)*dtheta)) + \
                (-self.g[2]*diffNeg*self.k[0]*self.A / (diffNeg*self.omega)) * \
                np.cosh(diffNeg*self.k[0]*(z+self.h))/np.cosh(diffNeg*self.k[0]*self.h) * \
                np.exp(1j*(i+1)*(theta-(i+1)*dtheta))
        # Deep water (kh >> 1)
        # ... TODO
                 
        # Shallow water (kh << 1)
        # ... TODO

        return np.real(u)
        # NOTE: implemented based on linearized ideal flow

    def velocity_v(self,x,t):
        v = 0.0
        return v
        # NOTE: you can implement based on linearized ideal flow   


    def velocity_w(self,x,t):
        theta =  self.k[0]*x[0] - self.omega*t + self.phase # ~ NOTE: x[0] is a vector here!
        dtheta = self.diff*theta        
        z = x[2] - self.h        
        w = np.zeros(x[0].shape)
        
        for i in range(self.N):
            diffPos = (1+self.diff*(i+1))
            diffNeg = (1-self.diff*(i+1))        
            w = w + -1j * (-self.g[2]*diffPos*self.k[0]*self.A/(diffPos*self.omega)) * \
                np.sinh(diffPos*self.k[0]*(z+self.h))/np.cosh(diffPos*self.k[0]*self.h) * \
                np.exp(1j*(i+1)*(theta+(i+1)*dtheta)) + \
                -1j * (-self.g[2]*diffNeg*self.k[0]*self.A/self.omega) * \
                np.sinh(self.k[0]*(z+self.h))/np.cosh(self.k[0]*self.h) * \
                np.exp(1j*(i+1)*(theta-(i+1)*dtheta))
        
        return np.imag(w)

    def pressure(self,x,t):
        theta =  self.k[0]*x[0] - self.omega*t + self.phase # ~ NOTE: x[0] is a vector here!
        dtheta = self.diff*theta          
        z = x[2] - self.h
        p = np.zeros(x[0].shape)

        for i in range(self.N):
            diffPos = (1+self.diff*(i+1))
            diffNeg = (1-self.diff*(i+1))
            p = p + self.rho_0*(-self.g[2])*self.A* \
                np.cosh(diffPos*self.k[0]*(z+self.h))/np.cosh(diffPos*self.k[0]*self.h) * \
                np.exp(1j*(i+1)*(theta+(i+1)*dtheta)) + \
                self.rho_0*(-self.g[2])*self.A* np.cosh(diffNeg*self.k[0]*(z+self.h))/np.cosh(diffNeg*self.k[0]*self.h) * \
                np.exp(1j*(i+1)*(theta-(i+1)*dtheta))
                
        return np.real(p)
        # NOTE: also implemented on ideal flow via linearized Bernoulli eqn.


class Solitary:
    """ Class that defines a solitary wave profile
        of a constant initial amplitude.
        
        .. todo:: Finish the docs. """
    def __init__(self,amplitude,omega,k,depth,rho_0,rho_1):
        self.name = 'Solitary'
        self.g = (0.0,0.0,-9.8)
        self.amp = amplitude        # ~ disregarded for now
        self.omega = omega
        self.k = k
        self.C = self.omega/self.k[0]   # phase speed
        self.h = depth
        self.rho_0 = rho_0              # density of water
        self.rho_1 = rho_1              # density of air         
        self.alpha = -0.390             # from Nwogu '93 (Boussineqs alpha ~ -1/3)
        
        self.A = np.abs(self.C**2 - (-self.g[2]*self.h)) / self.C
        self.B = np.sqrt( np.abs(self.C**2-(-self.g[2]*self.h)) / (4.0* ((self.alpha+1.0/3)*(-self.g[2]*self.h**3) - \
                                                             self.alpha*self.h**2*self.C**2)**2 ) )
        self.A1 = np.abs(self.C**2 - (-self.g[2]*self.h)) / ( 3.0* ((self.alpha+1.0/3)*(-self.g[2]*self.h) - \
            self.alpha*self.C**2) ) * self.h
        self.A2 = - np.abs(self.C**2 - (-self.g[2]*self.h))/(2.0*(-self.g[2])*self.h*self.C**2) * \
            ((self.alpha+1.0/3)*(-self.g[2]*self.h)+2.0*self.alpha*self.C**2)/ \
            ((self.alpha+1.0/3)*(-self.g[2]*self.h)-self.alpha*self.C**2) * self.h
    def height(self,x,t):
        T = t - 2.0 
        xi = x[0] - self.C*T #offsetting start time
        eta = self.A1/(np.cosh(self.B*xi))**2 + self.A2/(np.cosh(self.B*xi))**4
        # inflowHeightMean already added in wavetank.py
        #        
        # NOTE: x[0] is a vector here for verification only!
        return eta

    def velocity_u(self,x,t):
        T = t-2.0
        xi = x[0] - self.C*T # offsetting start time
        eta = self.height(x,t)
        u = self.C * eta / (eta+self.h)        
        return u

    def velocity_v(self,x,t):
        v = 0.0        # no trasverse dependece 1D
        return v
    
    def velocity_w(self,x,t):
        w = 0.0         # negligable for small amplitude (of order eps**2 ==> (ka)**2)
        return w
        # NOTE: base it on ideal fluid flow

    def pressure(self,x,t):
        p = 0.0       # P_atm
        return
        # NOTE: define via Bernoulli eqn.

class StokesWave:
    """ Class that defines a 2nd order Stokes wave (2nd in wave steepness),
        which included all 3-wave interactions and is an appropriate
        model for gravity-capillary waves.

        For steep open ocean gravity waves, need to extend to third
        order to capture all 4-wave interactions.
        
        NOTE: Stokes waves are only defined in 2D (x,z) - no y dependence.
              This expansion if for kh = O(1), but the limits of kh-->infinity,
              or kh-->0 can be readily derived as well.
        
        .. todo:: add third order correction if desired. """
        
    def __init__(self,amplitude,omega,k,depth,rho_0,rho_1):
        self.name = 'StokesWave'
        self.g = (0.0,0.0,-9.8)        
        self.A = amplitude
        self.omega = omega  # linearized dispersion relation
        self.k = k
        self.depth = depth
        self.rho_0 = rho_0      # density of water
        self.rho_1 = rho_1      # density of air 
        self.alpha = 1.0/np.tanh(self.k[0]*self.depth)
        self.omega_NL = np.sqrt( -self.g[2]*self.k[0]*np.tanh(self.k[0]*self.depth) * \
            (1 + (self.k[0]*self.A)**2*((9.0/8)*(self.alpha**2-1)**2 + self.alpha**2)) )
        # omega_NL contains second order nonlinear (NL) correction to phase speed
        self.scaled_A = -self.g[2]*self.A/self.omega_NL
        self.potential = 0.0

    def height(self,x,t):
        theta = self.k[0]*x[0] - self.omega_NL*t
        eta_0 = self.A * np.cos(theta)
        eta_1 = self.alpha/4.0 * (3.0*self.alpha**2+1)*self.A**2*self.k[0] * np.cos(2.0*theta)
        eta_2 = -3.0/8 * (self.alpha**4-3*self.alpha**2+3)*(self.A**3*self.k[0]**2)*np.cos(theta) + \
            3.0/64*(8*self.alpha**6+(self.alpha**2-1)**2)*(self.A**3*self.k[0]**2)*np.cos(3.0*theta)
        eta = eta_0 + eta_1 + eta_2
        # inflowHeightMean already added in wavetank.py
        #
        # NOTE: x[0] is a vector here! ===> change if you want to pass both x,t as scalars!
        return eta

    def potential(self,x,t):
        """ Velocity Potential correct up to second order.
            
            NOTE: here it is setup such that the mean depth is at L[2]/2 and the
                  wavetank bottom is located at z=0, otherwise if the mean
                  water depth=0, then we need x[2]==>x[2]+self.depth chenged below!
        """
        theta = self.k[0]*x[0] - self.omega_NL*t
        phi_0 = self.scaled_A * np.cosh(self.k[0]*x[2])/np.cosh(self.k[0]*self.depth) * np.sin(theta)
        phi_1 = 3.0/8 * self.A*self.k[0]*self.scaled_A/self.alpha * (self.alpha**2-1) * \
            np.cosh(2.0*self.k[0]*x[2]) * np.sin(2.0*theta)
        phi_2 = 1.0/64 * (self.alpha**2-1)*(self.alpha**2+3)*(9.0*self.alpha**2-13) * \
            np.cosh(3.0*self.k[0]*x[2])/np.cosh(3.0*self.k[0]*self.depth) * \
            (self.A*self.k[0])**2 * self.scaled_A * np.sin(3.0*theta) 
    
        self.potential = phi_0 + phi_1 + phi_2

    def velocity_u(self,x,t):
        """ Defined via potential flow: u = d/dx{potential}. """
        theta = self.k[0]*x[0] - self.omega_NL*t        
        u0 = self.k[0]* self.scaled_A * np.cosh(self.k[0]*x[2])/np.cosh(self.k[0]*self.depth) \
            * np.cos(theta)
        u1 = 2.0*self.k[0] * 3.0/8 * self.A*self.k[0]*self.scaled_A/self.alpha* \
            (self.alpha**2-1)*np.cosh(2.0*self.k[0]*x[2]) * np.cos(2.0*theta)
        u2 = 3.0*self.k[0]* 1.0/64 * (self.alpha**2-1)*(self.alpha**2+3)*(9.0*self.alpha**2-13) * \
            np.cosh(3.0*self.k[0]*x[2])/np.cosh(3.0*self.k[0]*self.depth) * \
            (self.A*self.k[0])**2 * self.scaled_A * np.cos(3.0*theta)
        u = u0 + u1 + u2
        return u

    def velocity_v(self,x,t):
        """ Stokes wave exists in 2D only (x,z) - no y dependence. """
        v = 0.0
        return v
    
    def velocity_w(self,x,t):
        """ Defined via potential flow: w = d/dz{potential}. """
        theta = self.k[0]*x[0] - self.omega_NL*t
        
        w0 = self.k[0] * self.scaled_A * np.sinh(self.k[0]*x[2])/np.cosh(self.k[0]*self.depth) \
            * np.cos(theta)
        w1 = 2.0*self.k[0] * 3.0/8 * self.A*self.k[0]*self.scaled_A/self.alpha* \
            (self.alpha**2-1)*np.sinh(2.0*self.k[0]*x[2]) * np.cos(2.0*theta)
        w2 = 3.0*self.k[0] * 1.0/64 * (self.alpha**2-1)*(self.alpha**2+3)*(9.0*self.alpha**2-13) * \
            np.sinh(3.0*self.k[0]*x[2])/np.cosh(3.0*self.k[0]*self.depth) * \
            (self.A*self.k[0])**2 * self.scaled_A * np.cos(3.0*theta)
        w = w0 + w1 + w2
        return w
        # NOTE: base it on ideal fluid flow

    def pressure(self,x,t):
        """ Defined through first order pressure based on lineqrized equations
            but with a second order correction to velocity potential such that
            
                pressure = rho*omega_NL/k * potential

            with P_atm = 0.
            """
        p = self.rho_0 * self.omega_NL / self.k[0] * self.potential
        return p
        # NOTE: define via Bernoulli eqn.


class waveJONSWAP:
    """ Class that defines a wave field based on realistic wave spectrum.
        The JONSWAP wave spectrum with directional distribution and
        random phases.
        
        .. todo:: Finish the docs. """

    def __init__(self,amplitude,depth,L):
        # Based on Linearized Dispersion Relations
        self.name = 'waveJONSWAP'        
        self.omegaPeak = 2.0*pi*JS.fp                  # peak angular frequency
        self.kp = (2.0*pi*JS.fp)**2 /JS.gv             # peak wavenumber 
        self.h = depth
        self.A = amplitude
        
        # Discretization
        self.Nx = 2**JS.npw1                           # number of Fourier modes
        self.Ny = 2**JS.npw2
        self.Nz = np.arange(0,L[2],0.01)                 # vertical disretization
        
        self.Lx = abs(JS.x1e - JS.x1b)                   # domain length
        self.Ly = abs(JS.y1e - JS.y1b)
        self.kc_x = 2.0*pi/self.Lx
        self.kc_y = 2.0*pi/self.Ly                        
        self.kc_x_modes = self.Nx*self.kc_x        
        self.kc_y_modes = self.Ny*self.kc_y

        self.surface = np.zeros((self.Nx,self.Ny))                      # free surface
        self.velPotential = np.zeros((self.Nx,self.Ny))#,self.Nz))      # velocity potential
        self.u = np.zeros((self.Nx,self.Ny))#,self.Nz))                 # u velocity
        self.v = np.zeros((self.Nx,self.Ny))#,self.Nz))                 # v velocity
        self.w = np.zeros((self.Nx,self.Ny))#,self.Nz))                 # w velocity
                

    def height(self,x,t):
        self.JONSWAP(x,t)
        # ~ NOTE: x[0] is a vector here!
        return self.surface

    def velocity_u(self,x,t):
        """
        Returns the velocity at the free-surface (z=surface)
            
        NOTE: to extract velocity at any height z, write down the
            power series and take fft() of for every z
            ... little time consuming, might need a bit of Cython here!

        The velocity potential defined at the free-surace :math:`z=\zeta` is
        given by :math:`\Phi(x,y,t) \equiv \phi(x,y,z=\zeta,t)`
    
        .. math:: 
        
            \phi(x,y,z,t) = \Phi + (z-\zeta)W + \sum ... + \sum ...

        where :math:`W` is the vertical velocity defined at the free-surface
        given through a Dirichlet-to-Neumann operator relation.
        """
        return self.u
        # NOTE: base it on linearized ideal fluid flow

    def velocity_v(self,x,t):
        """
        Returns the velocity at the free-surface (z=surface)
            
        NOTE: to extract velocity at any height z, write down the
            power series and take fft() of for every z
            ... little time consuming, might need a bit of Cython here!

        The velocity potential defined at the free-surace :math:`z=\zeta` is
        given by :math:`\Phi(x,y,t) \equiv \phi(x,y,z=\zeta,t)`
    
        .. math:: 
        
            \phi(x,y,z,t) = \Phi + (z-\zeta)W + \sum ... + \sum ...

        where :math:`W` is the vertical velocity defined at the free-surface
        given through a Dirichlet-to-Neumann operator relation.
        """        
        return self.v
        # NOTE: base it on linearized ideal fluid flow

    def velocity_v(self,x,t):
        """
        Returns the velocity at the free-surface (z=surface)
            
        NOTE: to extract velocity at any height z, write down the
            power series and take fft() of for every z
            ... little time consuming, might need a bit of Cython here!

        The velocity potential defined at the free-surace :math:`z=\zeta` is
        given by :math:`\Phi(x,y,t) \equiv \phi(x,y,z=\zeta,t)`
    
        .. math:: 
        
            \phi(x,y,z,t) = \Phi + (z-\zeta)W + \sum ... + \sum ...

        where :math:`W` is the vertical velocity defined at the free-surface
        given through a Dirichlet-to-Neumann operator relation.
        """        
        return self.w
        # NOTE: base it on linearized ideal fluid flow
        
    def pressure(self,x,t):
        """ Pressure was defined via linearized theory (Bernoulli Eqn.) """
        p = self.rho_0 * self.omega/self.k[0] * self.velPotential      # P_atm=0 
        return p
    
    def JONSWAP(self,x,t):
        """Sets a wave field according to JONSWAP ocean wave spectrum."""
        from spectrum import jonswap
        from potential import velocityPotential
       
        modes_x = np.fft.fftfreq(self.Nx)    # Fourier mode alignments
        modes_y = np.fft.fftfreq(self.Ny)

        kx = self.kc_x_modes * modes_x
        ky = self.kc_y_modes * modes_y

        # Generating mesh grid for 2D spectrum
        [kxx, kyy] = np.meshgrid(kx,ky)
    
        # Call jonswap() in spectrum.py module *** Send 1D vectors kx and ky (not kxx or kyy)!
        spectrum = jonswap(self.Nx, self.Ny, JS.Hs, JS.gamma, JS.fp, JS.nspread, JS.thetam, JS.gv, kx, ky)

        # (DISREGARD FOR NOW) - Imposing iniitial circular zero-pad filter
        if JS.filter == 1:
            spectrum[ 2*floor(self.kp/self.kc_x)+1 : (self.Nx-1)-2*floor(self.kp/self.kc_x), \
                2*floor(selfkp/kc_x)+1 : (self.Ny-1)-2*floor(self.kp/self.kc_x) ] = 0.0 + 0.0j 
        else:
            pass

        # Compute the surface elvation via iFFT2D of the spectrum  --> take REAL part
        # just in case we have some left over (nonzero) imaginary part of 'spectrum'
        self.surface = np.real(np.fft.ifft2(spectrum))

        # Compute the velocity potential from linear theory and horizontal velocity
        [self.velPotential, self.u, self.v, self.w] = velocityPotential(spectrum, JS.gv, kx, ky)


