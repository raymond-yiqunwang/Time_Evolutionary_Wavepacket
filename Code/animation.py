import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure(num=None, figsize=(8,8), dpi=100, facecolor='w', edgecolor='k')

ax = plt.axes(xlim=(-1.8, 1.8),ylim=(0,1))
line, = ax.plot([], [], lw=2, color = 'r')
ttl = ax.text(.70, .95, '', transform = ax.transAxes, va='center')
dt = 0.01
#matplotlib.animation.Animation._blit_draw = _blit_draw

# initialization function: plot the background of each frame
def init():
    ttl.set_text('')
    line.set_data([], [])
    return line,ttl

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(-1.5,1.5,100)
    nt = i+1
    os.system("sed -i '' '55s/.*/    do it = 1, %d /g' animation.f90" %nt)
    os.system("mpif90 -I ~/Documents/program/Fortran/Library/ ~/Documents/program/Fortran/Library/*.o /usr/local/lib/libfftw3.a -o av animation.f90")
    os.system("./av > data")
    y = np.loadtxt("dens_t.txt")
    ttl.set_text("Time = %f" %(i*dt))
    line.set_data(x, y)
    
    return line,ttl

# call the animator.  'blit=True' means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=100, interval=20, blit=True,repeat=False)
#writer = animation.writers['ffmpeg'](fps=30)
#anim.save('Grid.mp4',writer=writer,dpi=100)
# call our new function to display the animation
#matplotlib.display_animation(anim)a
plt.show()
