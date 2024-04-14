#                           MD Simulation Engine for a 2DLJ gas                           #
from tkinter import *   # needed to create the GUi window, the canvas and the particles
import numpy as np
import time
import random
import math

#  constants - taken either from the Gromacs Manual as indicated or recalculated, so they are correctly incorporated
#  into the simulation. There are constants for the widget and for the particles in the simulation.
WIDTH = 800
HEIGHT = 600
dt = 0.001 * pow(10, -12)
charge = 1.6 * pow(10, -19)
mass = 1.66 * pow(10, -27)
sigma = 1
epsilon = 1
g = 0.1
k = 1 / (4 * math.pi * epsilon)
particles_number = 20
# -----------------------------------------------------#

class Particle:
    """
    Particle class, here everything is done - creating the particle, calculating the forces,
    the collisions between the particles, their acceleration and updating their velocities all the time
    """
    def __init__(self, canvas, number):
        self.canvas = canvas
        self.number = number
        self.place = np.random.rand(number, 2) * HEIGHT
        self.g = g
        self.coords, self.xvel, self.yvel, self.charge = self.create()
        self.move()

    def forces(self):
        """
        This function calculates the Lennard-Jones force and the Coulomb force every two particles in the simulation.
        After getting the total force (Coulomb's + Lennard-Jones') it is split into x- and y-direction parts and from there.
        Newton's second law is then applied to calculate the acceleration in x- and y-direction which is then used to
        update the velocities of the particles.
        """
        for i, particle1 in enumerate(self.coords):
            charge_1 = self.charge[i]   # getting the charge of the particle
            for j, particle2 in enumerate(self.coords):
                charge_2 = self.charge[j]   # getting the charge of the particle

                # making sure that we do not calculate the same particle,
                # otherwise the calculations are not going to be possible
                if particle1 == particle2:
                    continue
                else:
                    # getting the coordinates of both particles
                    coords1 = canvas.coords(particle1)
                    coords2 = canvas.coords(particle2)
                    # calculating the where the centers of the atoms are
                    x_1 = (coords1[0] + coords1[2]) / 2
                    x_2 = (coords2[0] + coords2[2]) / 2
                    y_1 = (coords1[1] + coords1[3]) / 2
                    y_2 = (coords2[1] + coords2[3]) / 2
                    # calculating the distance connecting the two particles and the distance in x and y plane
                    dist = math.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2) * pow(10, -9)
                    dx = ((x_1 - x_2) / dist) * pow(10, -9)
                    dy = ((y_1 - y_2) / dist) * pow(10, -9)
                lennard = 48 * epsilon * ((sigma * pow(abs(dist), 12)) - 0.5 * (sigma * pow(abs(dist), 6))) / pow(abs(dist), 2)
                coulomb = (k * charge_1 * charge_2) / (pow(abs(dist), 3))
                # splitting the total force to an x- and y-component
                force_x = (lennard + coulomb) * dx
                force_y = (lennard + coulomb) * dy
                # updating the acceleration and the velocities
                a_x = force_x / mass
                a_y = force_y / mass
                self.xvel[i] += a_x * dt
                self.xvel[j] -= a_x * dt
                self.yvel[i] += (a_y + g) * dt * 0.5
                self.yvel[j] -= (a_y + g) * dt * 0.5
# -----------------------------------------------------#

    def collisions(self):
        """
        This function keeps track of whether two particles bump into each other, for simplicity the collisions are
        perfectly elastic. If the particles collide (the distance between them is equal or less a particle's diameter),
        their velocity is going to invert and thus will move in different direction than before the collision.
        """
        for i, particle1 in enumerate(self.coords):
            for j, particle2 in enumerate(self.coords):
                # making sure we are not calculating the same particle twice
                if particle1 == particle2: continue
                else:
                    # calculating the distance between the two particles using their coordinates
                    coords1 = self.canvas.coords(particle1)
                    coords2 = self.canvas.coords(particle2)
                    # calculating where the centers of the circles are
                    x_1 = (coords1[0] + coords1[2]) / 2
                    x_2 = (coords2[0] + coords2[2]) / 2
                    y_1 = (coords1[1] + coords1[3]) / 2
                    y_2 = (coords2[1] + coords2[3]) / 2
                    dist = math.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2)
                    #   if the distance is smaller than the diameter of the particles, collision takes place
                    if dist <= 10:  # 10 is the diameter of the particles, set on line 115 when creating the particles
                        self.xvel[i], self.xvel[j] = self.xvel[j], self.xvel[i]
                        self.yvel[i], self.yvel[j] = self.yvel[j], self.yvel[i]

# -----------------------------------------------------#

    def create(self):
        """
        This function creates the particles at random locations across the canvas with random charge, random initial
        for the x- and y-direction and stores them into lists used later in the code
        """
        cords_atoms = []
        velx_atoms = []
        vely_atoms = []
        charges = []
        for i, x_y in enumerate(self.place):
            x, y = x_y
            self.atom = self.canvas.create_oval(x - 5, y - 5, x + 5, y + 5, fill='blue')
            cords_atoms.append(self.atom)
            velx_atoms.append(random.randint(0,1))
            vely_atoms.append(random.randint(0, 1))
            charges.append(random.choice([-1, 1]) * charge)
        return cords_atoms, velx_atoms, vely_atoms, charges

    # -----------------------------------------------------#
    def move(self):
        """
        This function moves the particles across the canvas. This is done by updating the velocities according to the
        force acting upon the particles, the collisions between them and the periodic boundaries conditions are set.
        The PBCs are controlled by if statements in a for loop, checking every particle
        """

        for i, particle in enumerate(self.coords):
            coordinates = self.canvas.coords(particle)
            # for bringing the particle back into the frame when they leave from the right
            if coordinates[2] > self.canvas.winfo_width():
                self.canvas.moveto(particle, 0, (coordinates[1] + coordinates[3]) / 2)
            # for bringing the particle back into the frame when they leave from the left
            if coordinates[0] < 0:
                self.canvas.moveto(particle, 795, (coordinates[1] + coordinates[3]) / 2)
            # for bringing the particle back into the frame when they leave from the bottom
            if coordinates[3] > self.canvas.winfo_height():
                self.canvas.moveto(particle, (coordinates[0] + coordinates[2]) / 2, 0)
            # for bringing the particle back into the frame when they leave from the top
            if coordinates[1] < 0:
                self.canvas.moveto(particle, (coordinates[0] + coordinates[2]) / 2, 595)
            self.canvas.move(particle, self.xvel[i], self.yvel[i])
        self.forces()
        self.collisions()

# -----------------------------------------------------#

# initiating the window and the canvas
window = Tk()
window.title('MD simulation of LJ gas')
canvas = Canvas(window, width=WIDTH, height=HEIGHT, background='black')
canvas.pack()
atoms = Particle(canvas, particles_number)
# -----------------------------------------------------#

# starting the program
while True:
    atoms.move()    # calling the functions of the particles and also to create them
    window.update()
    time.sleep(0.0167)

window.mainloop()
# -----------------------------------------------------#