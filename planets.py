from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


class Planet():
    _registry = []

    def __init__(self, name, mass, xdistance, ydistance, xvelocity, yvelocity, size, colour):
        self._registry.append(self)
        self.name = name
        self.mass = mass
        self.xdistance = xdistance
        self.ydistance = ydistance
        self.xvelocity = xvelocity
        self.yvelocity = yvelocity
        self.size = size
        self.colour = colour
        self.xpos = []
        self.ypos = []
        self.forcex = 0
        self.forcey = 0
        self.old_forcex = 0
        self.old_forcey = 0
        self.point_pot = 0
        self.epot = []
        self.ekin = []
        self.sundist = []

G = 4 * np.pi*np.pi

def convert_to_solarmass(m):
    """ Celestial body mass in kg and returns it in terms of solar mass """
    new_mass = m/(1.9885*10**30)
    return new_mass


def fx(x,y,x_next,y_next,m,m_next):
    """ Return force on current celestial body from next one, x direction """
    G = 4 * np.pi*np.pi

    r = np.sqrt((x_next-x)**2 + (y_next-y)**2)

    if r == 0:
        value = 0

    else:
        value = G*m_next*m*(x_next-x)/(r*r*r)

    return value

def fy(x,y,x_next,y_next,m,m_next):
    """ Return force on current celestial body from next one, y direction """
    G = 4 * np.pi*np.pi

    r = np.sqrt((x_next-x)**2 + (y_next-y)**2)

    if r == 0:
        value = 0
    else:
        value = G*m_next*m*(y_next-y)/(r*r*r)

    return value

def calculate_potential(x,y,x_next,y_next,m, m_next):
    """ Return potential on current celestial body from next one """
    G = 4 * np.pi*np.pi

    r = np.sqrt((x_next-x)**2 + (y_next-y)**2)

    if r == 0:
        value = 0
    else:
        value = -G*m_next*m/r

    return value

def r(x,y,x_next,y_next):
    """ Calculates normed distance between two celestial bodies """
    value = np.sqrt((x_next-x)**2 + (y_next-y)**2)
    return value

m_asteroid = 1.9e27

sunmass = convert_to_solarmass(1.9885*10**30)
mercmass = convert_to_solarmass(3.30e23)
venusmass = convert_to_solarmass(4.87e24)
earthmass = convert_to_solarmass(5.97e24)
marsmass = convert_to_solarmass(6.42e23)
jupitermass = convert_to_solarmass(1.90e27)
saturnmass = convert_to_solarmass(5.68e26)
uranusmass = convert_to_solarmass(8.68e25)
neptunemass = convert_to_solarmass(1.02e26)
halleysmass = convert_to_solarmass(2.2e14)
asteroidmass = convert_to_solarmass(m_asteroid)


                # Name       Mass        x      y   vx vy    size  colour
sun     = Planet("Sun",     sunmass,     0,     0,  0, 0,     5,   "gold")
mercury = Planet("Mercury", mercmass,    0.387, 0,  0, 10.09, 3,   "grey")
venus   = Planet("Venus",   venusmass,   0.723, 0,  0, 7.386, 5,   "darkorange")
earth   = Planet("Earth",   earthmass,   1,     0,  0, 6.28, 6,    "steelblue")
mars    = Planet("Mars",    marsmass,    1.52,  0,  0, 5.08, 5,    "peru")
jupiter = Planet("Jupiter", jupitermass, 5.20,  0,  0, 2.746, 12,  "sandybrown")
saturn  = Planet("Saturn",  saturnmass,  9.58,  0,  0, 2.047, 10,  "wheat")
uranus  = Planet("Uranus",  uranusmass,  19.2,  0,  0, 1.4413, 9,  "lightsteelblue")
neptune = Planet("Neptune", neptunemass, 30.05, 0,  0, 1.153, 8,   "royalblue")


# asteroid = Planet("Halley's Comet", halleysmass, 15, 15, -0.8, -1, 4, "lightskyblue" ) #Good set of initial conditions
asteroid = Planet("Asteroid", asteroidmass, -10, 15, 6.23, -4.998, 4, "dimgrey" ) #Close flyby of Jupiter

#asteroid = Planet("Asteroid", asteroidmass, 100, 100, 6.23, -5, 4, "dimgrey" )



# parameter controlling plot interval: do not plot too often not to slow down the animation
t = 0.	             # start time
tmax = 200.	         # final time
dt = 0.002           # time step
counter = 0

stepsperframe = 25
numframes     = int(tmax/(stepsperframe*dt))

plotSize = 35       # Set x- and y-limits for animation plot

###### Set font sizes for legends and axis #####
plt.style.use('seaborn')
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['axes.labelsize'] = 19
plt.rcParams['axes.titlesize'] = 19
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['figure.titlesize'] = 19

##### Initialize figure for animation
fig = plt.figure(figsize=(10,8))
ax  = plt.subplot(xlim=(-plotSize, plotSize), ylim=(-plotSize, plotSize))
plt.title("Asteroid mass = " + str(m_asteroid) + " kg")
plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
plt.axhline(y=0, c = "steelblue")    # draw a default hline at y=1 that spans the xrange
plt.axvline(x=0, c = "steelblue")    # draw a default vline at x=1 that spans the yrange

##### Set size and colour of every body
sun_plot = ax.plot([],[], 'o', c = sun.colour, markersize = sun.size )[0]
mercury_plot = ax.plot([],[], 'o', c = mercury.colour, markersize = mercury.size )[0]
venus_plot = ax.plot([],[], 'o', c = venus.colour, markersize = venus.size )[0]
earth_plot = ax.plot([],[], 'o', c = earth.colour, markersize = earth.size )[0]
mars_plot = ax.plot([],[], 'o', c = mars.colour, markersize = mars.size )[0]
jupiter_plot = ax.plot([],[], 'o', c = jupiter.colour, markersize = jupiter.size )[0]
saturn_plot = ax.plot([],[], 'o', c = saturn.colour, markersize = saturn.size )[0]
uranus_plot = ax.plot([],[], 'o', c = uranus.colour, markersize = uranus.size )[0]
neptune_plot = ax.plot([],[], 'o', c = neptune.colour, markersize = neptune.size )[0]
asteroid_plot = ax.plot([],[], 'o', c = asteroid.colour, markersize = asteroid.size )[0]

years_text = ax.text(10,30, "initialize string", fontsize = 12)  # Text to display number of earth years passed in animation


def init():
    """ Initialize animation for all bodies in system """"
    sun_plot.set_data([], [])
    mercury_plot.set_data([], [])
    venus_plot.set_data([], [])
    earth_plot.set_data([], [])
    mars_plot.set_data([], [])
    jupiter_plot.set_data([], [])
    saturn_plot.set_data([], [])
    uranus_plot.set_data([], [])
    neptune_plot.set_data([], [])
    asteroid_plot.set_data([], [])
    years_text.set_text("initialize string")

    return sun_plot, mercury_plot, venus_plot, earth_plot, mars_plot,
    jupiter_plot, saturn_plot, uranus_plot, neptune_plot, years_text, asteroid_plot

def velocity_verlet():
    """ Performs one step of velocity Verlet integrator of entire solar system """
    global t, counter

    for planet in Planet._registry:

        planet.xdistance += planet.xvelocity*dt + 0.5*planet.forcex/planet.mass*dt**2
        planet.ydistance += planet.yvelocity*dt + 0.5*planet.forcey/planet.mass*dt**2

        planet.old_forcex = planet.forcex
        planet.old_forcey = planet.forcey

        planet.forcex = 0
        planet.forcey = 0


    for planet in Planet._registry:

        foo_energy = 0

        for planet2 in Planet._registry:

            planet.forcex += fx(
                planet.xdistance, planet.ydistance, planet2.xdistance,
                planet2.ydistance, planet.mass, planet2.mass)

            planet.forcey += fy(
                planet.xdistance, planet.ydistance, planet2.xdistance,
                planet2.ydistance, planet.mass, planet2.mass)

            foo_energy += calculate_potential(
                planet.xdistance, planet.ydistance, planet2.xdistance,
                planet2.ydistance, planet.mass, planet2.mass)

        planet.epot.append(foo_energy)

    for planet in Planet._registry:

        planet.xvelocity += 0.5*(planet.forcex + planet.old_forcex)/planet.mass*dt
        planet.yvelocity += 0.5*(planet.forcey + planet.old_forcey)/planet.mass*dt

        v2  = planet.xvelocity**2 + planet.yvelocity**2

        planet.ekin.append( 0.5*planet.mass*v2 )

        # Append position in position list
        planet.xpos.append(planet.xdistance)
        planet.ypos.append(planet.ydistance)


    if counter % round(1/dt) == 0:
        earth.sundist.append(r(earth.xdistance, earth.ydistance, sun.xdistance, sun.ydistance))

    counter += 1
    t += dt


def euler():
    """ Performs one step of the Euler-Cromer integrator of entire solar system """
    global t, counter


    for planet in Planet._registry: # Update force on every celestial body

        planet.forcex = 0
        planet.forcey = 0

        for planet2 in Planet._registry:

            planet.forcex += fx(
                planet.xdistance, planet.ydistance, planet2.xdistance,
                planet2.ydistance, planet.mass, planet2.mass)

            planet.forcey += fy(
                planet.xdistance, planet.ydistance,  planet2.xdistance,
                planet2.ydistance, planet.mass, planet2.mass)

    for planet in Planet._registry: # Update position and velocity on all celestial bodies

        planet.xvelocity += dt*planet.forcex/planet.mass
        planet.yvelocity += dt*planet.forcey/planet.mass

        planet.xdistance += planet.xvelocity*dt
        planet.ydistance += planet.yvelocity*dt

        # Calculate kinetic energy and append it to kinetic energy list
        v2  = planet.xvelocity**2 + planet.yvelocity**2
        planet.ekin.append( 0.5*planet.mass*v2 )

        # Append position in position list
        planet.xpos.append(planet.xdistance)
        planet.ypos.append(planet.ydistance)


    for planet in Planet._registry: # Calculate potential energy on every celestial body
        foo_energy = 0

        # Calculate total potential energy on every body in system
        for planet2 in Planet._registry:
            foo_energy += calculate_potential(
                planet.xdistance, planet.ydistance, planet2.xdistance,
                planet2.ydistance, planet.mass, planet2.mass)

        planet.epot.append(foo_energy)

    if counter % round(1/dt) == 0:
        earth.sundist.append(r(earth.xdistance, earth.ydistance, sun.xdistance, sun.ydistance))


    counter += 1
    t += dt


def animate(framenr):
    """ Animation function. Plots position of every planet continously"""
    for it in range(stepsperframe):
        velocity_verlet()
        #euler()

    sun_plot.set_data(sun.xdistance, sun.ydistance)
    mercury_plot.set_data(mercury.xdistance, mercury.ydistance)
    venus_plot.set_data(venus.xdistance, venus.ydistance)
    earth_plot.set_data(earth.xdistance, earth.ydistance)
    mars_plot.set_data(mars.xdistance, mars.ydistance)
    jupiter_plot.set_data(jupiter.xdistance, jupiter.ydistance)
    saturn_plot.set_data(saturn.xdistance, saturn.ydistance)
    uranus_plot.set_data(uranus.xdistance, uranus.ydistance)
    neptune_plot.set_data(neptune.xdistance, neptune.ydistance)
    asteroid_plot.set_data(asteroid.xdistance, asteroid.ydistance)

    # Prints number of earth years passed
    years_text.set_text("Number of years passed: " + str(round(framenr*dt*stepsperframe,2)))

    return sun_plot, mercury_plot, venus_plot, earth_plot, mars_plot,
    jupiter_plot, saturn_plot, uranus_plot, neptune_plot, asteroid_plot, years_text

#Call the animator, blit=True means only re-draw parts that have changed
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=numframes, interval=50, blit=True, repeat=False)

#Loop for long calculations, should be used without animation
# for it in range(stepsperframe*numframes):
#     print(str(it/(stepsperframe*numframes)*100) + " %")
#     #euler()
#     velocity_verlet()

plt.show()
fig = plt.figure(figsize=(10,8))

# Loop through and plot all bodies in Planet class
for planet in Planet._registry:
    plt.plot(planet.xpos, planet.ypos, label = planet.name, c = planet.colour)

plt.legend(loc='best')  # Sets location of legend

years = 0
for i in range(len(earth.xpos)): # Calculate the number of earth years passed
    if earth.ypos[i-1] < 0 and earth.ypos[i] >= 0 :
        years += 1

plt.xlim(-plotSize, plotSize)
plt.ylim(-plotSize, plotSize)
plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
plt.title("Trajectory of asteroid. " + str(years) + " earth years completed. dt = "  + str(dt))
plt.show()
