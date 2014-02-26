import sys
import os
import subprocess
from test_waveModules_p import *
import defaultWaveModules as wm
import JONSWAP_p as JS
import numpy as np
import argparse

try:
    import matplotlib.pyplot as plt
    from matplotlib import cm, animation
    import mpl_toolkits.mplot3d.axes3d as p3
except:
    print "[EXCEPTION RAISED] - Cannon IMPORT and initialize plotting with Matplotlib - EXIT"
    sys.exit()

""" Testing different analytic solution for wavetank benchmark case."""

#see nose docs for more complex testing

def test_Linear2D(args):
    """ Testing the Linearized 2D interface (phi) propagation. """
    # Collocation point
    N = 100
    NTime = 50
    tFinal = 5.0*period
    x = [np.linspace(0,L[0],N), 0.0, np.linspace(0,L[2],N)]
    time = np.linspace(0,tFinal,NTime)

    # Wave Field Object
    waveTest = wm.Linear2D(amplitude,omega,
                           k,inflowHeightMean,
                           rho_0,rho_1,randomPhase)

    # test flow control
    print "TEST CASE NAME is: ", waveTest.name
    saveFrames = args.save
    saveMovie = args.movie
    playMovie = args.play
    frameRate = NTime/tFinal # for saving movie with ffmpeg (below)
    if ~os.path.exists(waveTest.name):
        subprocess.call(["mkdir",waveTest.name])
    
    # monitoring the progress
    sys.stdout.write('Processing.')
    # get free-surface and velocity vector vel=(u,0,w)
    for [n,T] in enumerate(time):
        sys.stdout.write('.')
        sys.stdout.flush() #flush the stdout buffer
        surface = waveTest.h + waveTest.height(x,T) 
        [uu,ww] = np.meshgrid(np.zeros(x[0].size), np.zeros(x[2].size))
        speed = np.zeros([x[0].size,x[2].size])

        # velocity streamlines setup
        for [j,xj] in enumerate(x[0]):
            for [i,zi] in enumerate(x[2]):
                XZ_location = [xj,0.0,zi]
                if zi < surface[j]:
                    U = waveTest.velocity_u(XZ_location,T) # x is expected as 3D vector x=>[x,y,z]
                    W = waveTest.velocity_w(XZ_location,T)
                    uu[i,j] = U 
                    ww[i,j] = W 
                    speed[i,j] = np.sqrt(U**2 + W**2)
                else:
                    uu[i,j] = 0.0 # note the row/column alignemnt here
                    ww[i,j] = 0.0
                    speed[i,j] = 0.0
    
        # plotting surface and velocity streamlines
        fig = plt.figure()
        plt.clf() #clear figure
        surfacePlot = plt.plot(x[0],surface)
        speedPlot = plt.streamplot(x[0],x[2],uu,ww,
                               density=(1,1),
                               color=speed,
                               linewidth=2.0*speed/speed.max())
        plt.axis((0.0, L[0], 0.0, L[2]))
            
        if saveFrames:
            plt.savefig(waveTest.name+'/'+'frame%4.4d.png' % n)
        else:
            plt.show() # user wants to view/plot every frame (time consuming)
    print "[DONE]"

    if saveFrames and saveMovie:
        try:
            subprocess.call(r"ffmpeg -r "+str(frameRate)+" -y -i "+waveTest.name+"/frame%4d.png -vcodec libx264 -sameq "+waveTest.name+"/"+waveTest.name+".mp4",shell=True)
        except:
            print "[EXEPTION RAISED] - Cannot save "+waveTest.name+".mp4 with ffmpeg"
            sys.exit()

    if saveMovie and playMovie:
        try:
            subprocess.call(r"open "+waveTest.name+"/"+waveTest.name+".mp4", shell=True)
        except:
            print "[EXCEPTION RAISED] - Cannot play "+waveTest.name+".mp4 movie"
            sys.exit()


def test_Linear3D(args):
    """ Testing the Linearized 2D interface (phi) propagation. """

    # Collocation points
    N = 100
    Nx = N; Ny = N
    NTime = N
    tFinal = 5.0*period
    x = [np.linspace(0,L[0],Nx), np.linspace(0,L[1],Ny), np.linspace(0,L[2],N)]
    time = np.linspace(0,tFinal,NTime)
    [xx,yy,zz] = np.meshgrid(x[0],x[1],x[2])
    #xx = x[0]
    #yy = x[1].reshape(-1,1)

    # Wave Field Object
    waveTest = wm.Linear3D(amplitude,omega,
                           k,inflowHeightMean,
                           rho_0,rho_1, randomPhase)

    # test flow control
    print "TEST CASE NAME is: ", waveTest.name
    saveFrames = args.save
    saveMovie = args.movie
    playMovie = args.play
    frameRate = NTime/tFinal # for saving movie with ffmpeg (below)
    if ~os.path.exists(waveTest.name):
        subprocess.call(["mkdir",waveTest.name])
    


        #
        #  LEFT OFF HERE 1/31/2014 
        #


    # monitoring the progress
    sys.stdout.write('Processing.')
    # get free-surface and velocity vector vel=(u,v,w)
    for [n,T] in enumerate(time):
        sys.stdout.write('.')
        sys.stdout.flush() #flush the stdout buffer
        surface = waveTest.h + waveTest.height(x,T) 
        [uu,ww] = np.meshgrid(np.zeros(x[0].size), np.zeros(x[2].size))
        speed = np.zeros([x[0].size,x[2].size])

        # velocity streamlines setup
        #for [j,xj] in enumerate(x[0]):
        #    for [i,zi] in enumerate(x[2]):
        #        XZ_location = [xj,0.0,zi]
        #        if zi < surface[j]:
        #            U = waveTest.velocity_u(XZ_location,T) # x is expected as 3D vector x=>[x,y,z]
        #            W = waveTest.velocity_w(XZ_location,T)
        #            uu[i,j] = U 
        #            ww[i,j] = W 
        #            speed[i,j] = np.sqrt(U**2 + W**2)
        #        else:
        #            uu[i,j] = 0.0 # note the row/column alignemnt here
        #            ww[i,j] = 0.0
        #            speed[i,j] = 0.0
    
        # plotting surface and velocity streamlines
        fig = plt.figure()
        ax = fig.gcd(projection='3d')
        plt.clf() #clear figure
        surface_temp = np.sin(xx+yy)
        surfacePlot = ax.plot_trisurf(xx,yy,surface_temp,
                                      cmap=cm.jet, linewidth=0.1)
        #speedPlot = plt.streamplot(x[0],x[2],uu,ww,
        #                       density=(1,1),
        #                       color=speed,
        #                       linewidth=2.0*speed/speed.max())
        #plt.axis((0.0, L[0], 0.0, L[2]))
        plt.show()

        sys.exit()

        if saveFrames:
            plt.savefig(waveTest.name+'/'+'frame%4.4d.png' % n)
        else:
            plt.show() # user wants to view/plot every frame (time consuming)
    print "[DONE]"

    if saveFrames and saveMovie:
        try:
            subprocess.call(r"ffmpeg -r "+str(frameRate)+" -y -i "+waveTest.name+"/frame%4d.png -vcodec libx264 -sameq "+waveTest.name+"/"+waveTest.name+".mp4",shell=True)
        except:
            print "[EXEPTION RAISED] - Cannot save "+waveTest.name+".mp4 with ffmpeg"
            sys.exit()

    if saveMovie and playMovie:
        try:
            subprocess.call(r"open "+waveTest.name+"/"+waveTest.name+".mp4", shell=True)
        except:
            print "[EXCEPTION RAISED] - Cannot play "+waveTest.name+".mp4 movie"
            sys.exit()


    # Plot result if appropriate
    #assert correctResult.all() == result.all(), "Linear2D.height returned %f should be %f" % (result,correctResult)

    fig = plt.figure()
    #y = waveTest.height(x,t[0]) 
    #line, = plt.plot(x[0],y)
    #plt.axis((0.0, L[0], 0, L[2]))

    def f(xx,yy):
        return np.sin(xx)  + np.cos(yy) #waveTest.height(x,t[count])

    im = plt.imshow(f(xx,yy), cmap=plt.get_cmap('jet'))

    def updatefig(*args):
        global xx,yy#x,t
        
        im.set_array(f(x,t,count))

    #def animate(i):
    #    line.set_ydata(waveTest.height(x,t[i]))  # update the data
    #    return line,

    #def init():
    #    line.set_ydata(np.ma.array(x[0], mask=True))
    #    return line,

    # Get initial wave field for all x,y at t=0
    #result = np.transpose(waveTest.height(x,t[0]))

    #try:
        # Attaching 3D axis to the figure
        #fig = plt.figure()
        #ax = p3.Axes3D(fig)
        
        #surf = ax.plot_surface(xx,yy,result,rstride=2,cstride=2, cmap=cm.jet,
        #        linewidth=0.5, antialiased=False)
        #ax.plot_wireframe(xx,yy,surface, rstride=4, cstride=4)

    #    anim = animation.FuncAnimation(fig, animate, np.arange(1, 200), init_func=init,
    #        interval=25, blit=True)    

        #anim.save('monochromatic_wave_BC_animation.mov', fps=20)
     #   plt.show()
    #except:
     #   pass


def test_WaveGroup(args):
    """ Testing the Linearized 2D interface (phi) propagation. """

    # Collocation points
    N = 100
    x = [np.linspace(0,L[0],N), 0.0, 0.0]
    t = np.linspace(0,10*period,N)
    [xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    waveTest = wm.WaveGroup(amplitude,omega,k,inflowHeightMean,
                            rho_0,rho_1,randomPhase)

    # Plot result if appropriate
    #assert correctResult.all() == result.all(), "Linear2D.height returned %f should be %f" % (result,correctResult)

    fig = plt.figure()
    y = waveTest.height(x,t[0])
    line, = plt.plot(x[0],y)
    plt.axis((0.0, L[0], 0.0, L[2]))

    def animate(i):
        line.set_ydata(waveTest.height(x,t[i]))  # update the data
        return line,

    def init():
        line.set_ydata(np.ma.array(x[0], mask=True))
        return line,
        
    try:
        anim = animation.FuncAnimation(fig, animate, np.arange(1, 1000), init_func=init,
            interval=2, blit=True)    

        #anim.save('modulated_waveGroup_BC_animation.mov', fps=30)
        plt.show()
    except:
        pass


def test_Solitary(showPlots=False):
    """ Testing the Linearized 2D interface (phi) propagation. """

    # Collocation points
    N = 100
    x = [np.linspace(0,L[0],N), 0.0, 0.0]
    t = np.linspace(0,10*period,N)
    [xx, tt] = np.meshgrid(x,t) 
    
    # Wave Field Object
    waveTest = wm.Solitary(2.0*A,omega,k,h,rho_0,rho_1)

    fig = plt.figure()
    y = waveTest.height(x,t[0])
    line, = plt.plot(x[0],y)
    plt.axis((0.0, L[0], 0, L[2]))

    def animate(i):
        line.set_ydata(waveTest.height(x,t[i]))  # update the data
        return line,

    def init():
        line.set_ydata(np.ma.array(x[0], mask=True))
        return line,
        
    try:
        anim = animation.FuncAnimation(fig, animate, np.arange(1, 500), init_func=init,
            interval=20, blit=True)    

        #anim.save('solitary_wave_BC_animation.mov', fps=30)
        plt.show()
    except:
        pass


def test_StokesWave(showPlots=False):
    """ Testing the 2nd order nonlinear Stokes Wave solution. """

    # Collocation points
    N = 100
    x = [np.linspace(0.0,L[0],N), 0.0, np.linspace(0.0,L[2],N)]
    time = np.linspace(0.0,10*period,N)
    
    # Wave Field Object
    waveTest = wm.StokesWave(amplitude,omega,k,inflowHeightMean,rho_0,rho_1)

    # check initial surface and velocity
    T = time[0] # initial time for now => no time loop yet
    surface = waveTest.depth + waveTest.height(x,T) # add depth waveTest.h if desired later for plotting!
    [uu,ww] = np.meshgrid(np.zeros(x[0].size), np.zeros(x[2].size))
    speed = np.zeros([x[0].size,x[2].size])

    # velocity streamline setup
    for [j,xj] in enumerate(x[0]):
        for [i,zi] in enumerate(x[2]):
            XZ_location = [xj,0.0,zi]
            #surface = waveTest.height(x,t[0]) 
            if zi < surface[j]:
                U = waveTest.velocity_u(XZ_location,T) # x is expected as 3D vector x=>[x,y,z]
                W = waveTest.velocity_w(XZ_location,T)
                uu[i,j] = U # *** check row/column alignment here!!!
                ww[i,j] = W 
                speed[j,i] = np.sqrt(U**2 + W**2)
            else:
                uu[i,j] = 0.0 #check row/column alignment
                ww[i,j] = 0.0
                speed[i,j] = 0.0
    
    # plotting surface and velocity streamlines
    fig = plt.figure()
    surfacePlot = plt.plot(x[0],surface)
    speedPlot = plt.streamplot(x[0],x[2],uu,ww,
                               density=(1,1),
                               color=speed,
                               linewidth=5.0*speed/speed.max())
    plt.axis((0.0, L[0]/2, 0.0, L[2]))
    plt.show() 

def test_waveJONSWAP(args):
    """ Testing the initialization via JONWAP wave spectrum."""
    # Wave Field Object
    waveTest = wm.waveJONSWAP(amplitude,h,L)

    # Discretization
    Nx = 2**JS.npw1                # number of Fourier modes
    Ny = 2**JS.npw2
    kc_x = 2*np.pi/L[0]
    kc_y = 2*np.pi/L[1]              # std. wave factors if modes would be integers
    kc_x_modes = Nx*kc_x        
    kc_y_modes = Ny*kc_y        # wave factors for modes obtained from fftfreq()

    # Collocation point
    N = 100
    x = [np.linspace(0,L[0],Nx), np.linspace(0,L[1],Ny), 0.0]
    t = np.linspace(0,10*period,N)
    [xx, yy] = np.meshgrid(x[0],x[1])

    # Get initial wave field for all x,y at t=0
    result = np.transpose(waveTest.height(x,t[0]))

    try:
        # Attaching 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)
        
        surf = ax.plot_surface(xx,yy,result,rstride=2,cstride=2, cmap=cm.jet,
                linewidth=0.5, antialiased=False)
        #ax.plot_wireframe(xx,yy,surface, rstride=4, cstride=4)
        plt.show()
    except:
        pass
    

if __name__ == '__main__':
    # parsing the input arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-s','--save', action='store_true',
                         help='save the PNG frames (default: no saving if flag not set)')
    parser.add_argument('-m','--movie', action='store_true', 
                        help='create a movie from frames with ffmpeg (default: no MP4 movie created in flag not set)')
    parser.add_argument('-p','--play', action='store_true',
                        help='play a movie created with ffmpeg (default: movie not played if flag not set)')
    args = parser.parse_args()
    print 'Input arguments are: ', args

    #test_Linear2D(args)
    test_Linear3D(args) #---> NEED TO FINISH
    #test_WaveGroup(args)
    #test_Solitary(args)
    #test_StokesWave(args)
    #test_waveJONSWAP(args)
