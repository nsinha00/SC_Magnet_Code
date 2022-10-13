# NOTE: Modification for a sloped bobbin, with increasing no. of turns in each layer (this kind of shape: \____/, as opposed to |____|) -- Edit by Nish
# strategy: increase L iteratively for each layer
# want length and field value
# strategy: fix slope and vary field, then fix field and vary length

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 14:56:09 2021

@author: alexlindgrenruby
"""
import math
import numpy as np
import matplotlib.pyplot as plt

#Script contains functions to calculate components of vec field B for solenoid

#Experimental parameters
I = 20 #Amps
N = 1780 #no. turns
R = 28e-3 #m, radius of coil
L = 25e-3 #m, length of coil # Note: will change iteratively
n = N/L #1/m, turn density # Note: will change iteratively

#winding parameters
wire_d = 0.279e-3 #m, wire diameter
n_l = np.floor(L/wire_d) # wires per 2
layers = np.ceil(N/n_l) # number of layers
remain = N%n_l #number of turns on the outer layer
coil_outer_rad = R+(layers*wire_d)

# Code added by Nish
#Parameters specifically for Slanted Bobbin
L_min = 25e-3 #m, min length of coil #! may want to vary this and keep L_max at our current value (25e-3 m)
L_max = 30e-3 #m, max length of coil #* Get correct value
height = math.sqrt(2)*20*wire_d #m, height of coil
total_wire_len = 342 #m, total lenth of wire (for slanted bobbin functions) #* Get correct value
grade = height / ( (L_max - L_min)/2 ) #grade of bobbin slope
#! Fix layer count, vary bobbin dimensions
#! Plot field w/r/t angle

mu = 1.25663706 * 10 ** -6 #m*kg*s**-2*A**-2, permitivity of free space


#Computational parameters
arcs = 100 #no of divisions for the angle summations
sub_angle = 2*math.pi / arcs #subtended angle corrosponding to each arc segment 
Pr = 0 #only consider points along the z axis
phi = 1 #dummy variable, not used (unless off-axis points are considered)

#I'm not convinced that this function does what it is supposed to
def fieldxy(Pz):
    #B field x component
    Bxy = 0 #initialize x component
    for m in range(N):
        for theta in range(arcs):
            zmth = -0.5*L + m/n
            theta = theta * sub_angle
            Bxy += ((I * R * (Pz - zmth) * sub_angle) / 
                   (Pr**2  + (Pz - zmth)**2 + R**2 - 2*Pr*R*math.cos(theta-phi))**(3/2))
    Bxy = Bxy * mu / (4 * math.pi)
    #why is the loop variable theta being used here
    Bx = Bxy * math.cos(theta)
    By = Bxy * math.sin(theta)
    return Bx, By
    
def fieldz(Pz):
    #B field z component 
    Bz= 0 #initialize z component
    for m in range(N):
        for theta in range(arcs):
            zmth = -0.5*L + m/n
            theta = theta * sub_angle
            Bz += (-(I * R * (Pr * math.cos(theta - phi) - R) * sub_angle) / 
                   (Pr**2  + (Pz - zmth)**2 + R**2 - 2*Pr*R*math.cos(theta-phi))**(3/2))
    Bz = Bz * mu / (4 * math.pi)
    return Bz

def field(Pz):
    #compound function - returns vector field at given position
    Bx, By = fieldxy(Pz)
    Bz = fieldz(Pz)
    B = [Bx, By, Bz]
    return B

#code added by Nolan
def fieldz2(Pz):
    #Check the Bz field component, since getting significant nonzero Bx field
    Bz= 0 #initialize z component
    for m in range(N):
        zmth = -0.5*L + m/n
        Bz += (I * R**2) / ((Pz - zmth)**2 + R**2)**(3/2)
    Bz = Bz * mu / 2.0
    return Bz
    
def wire_len(rad, sol_len, layer):
    #find the number of turns per layer
    N_layer = int(np.floor(sol_len/wire_d))
    #calculates the length of superconducting wire needed to make the magnet
    total_len = 0
    for i in range(int(layer)):
        circumf = 2*math.pi*(rad+i*wire_d)
        total_len += N_layer*circumf
    return total_len

def field_dist_plotter(final_dist):
    #plots the magenetic field strength as a function of distance from the end
    # of the magnet assembly
    
    #calculate
    distances = np.linspace(0,final_dist,50)
    field_v_dist = np.array([fieldz(dist+L*0.5+2.5e-3) for dist in distances])
    
    #plot
    plt.figure()
    plt.plot(distances*1e2, field_v_dist)
    plt.xlabel("Distance from Magnet (cm)")
    plt.ylabel("Magentic Field (T)")
    plt.show()

#parameterized field function
def fieldz_param(Pz_offset, rad, sol_len, layer):
    #find the number of turns per layer
    N_layer = int(np.floor(sol_len/wire_d))
    
    #fix the distance from the end of the magnet
    Pz = Pz_offset+0.5*sol_len+2.5e-3
    
    #B field z component 
    Bz= 0 #initialize z component
    for i in range(layer):
        layer_rad = rad+(i*wire_d)
        for m in range(N_layer):
            for theta in range(arcs):
                zmth = -0.5*sol_len + m*wire_d
                theta = theta * sub_angle
                Bz += (-(I * layer_rad * (Pr * math.cos(theta - phi) - layer_rad) * sub_angle) / 
                       (Pr**2  + (Pz - zmth)**2 + layer_rad**2 - 2*Pr*layer_rad*math.cos(theta-phi))**(3/2))
    Bz = Bz * mu / (4 * math.pi)
    return Bz


#parameterized field function -- SLANTED BOBBIN -- code added by Nish
def fieldz_per_layer(Pz_offset, rad, layer_num):
    L_layer = L_min + 2*layer_num*wire_d / grade
    #find the number of turns per layer
    N_layer = int(np.floor(L_layer/wire_d))
    
    #fix the distance from the end of the magnet
    Pz = Pz_offset+0.5*L_layer+2.5e-3 #? ignore for now?
    
    #B field z component 
    Bz= 0 #initialize z component
    layer_rad = rad+(layer_num*wire_d)
    for m in range(N_layer):
        for theta in range(arcs):
            zmth = -0.5*L_layer + m*wire_d
            theta = theta * sub_angle
            Bz += (-(I * layer_rad * (Pr * math.cos(theta - phi) - layer_rad) * sub_angle) / 
                    (Pr**2  + (Pz - zmth)**2 + layer_rad**2 - 2*Pr*layer_rad*math.cos(theta-phi))**(3/2))
    Bz = Bz * mu / (4 * math.pi)
    return Bz

def len_per_layer(rad, layer_num):
    L_layer = L_min + 2*layer_num*wire_d / grade
    num_turns = np.floor(L_layer / wire_d)
    return 2*math.pi*(rad + layer * wire_d) * num_turns

# calculate total Bz
used_len = 0
layer = 1
Bz = 0
while used_len < total_wire_len:
    #add B field at each layer
    Bz += fieldz_per_layer(0, R, layer)
    #iterate to the next layer
    used_len += len_per_layer(R, layer)
    layer += 1

print(layer)
print(Bz)


#desired function calls
#field_dist_plotter(5e-2)
print("Wire Length Required: {}[m]".format(wire_len(R, L, layers)))

#do some parameter optimization
#param_size = 9
#rads = np.linspace(2e-2,3e-2, param_size)
#lens = np.linspace(1e-2,3e-2, param_size)
#field_array = np.zeros([param_size, param_size])
#length_array = np.zeros([param_size, param_size])
#for i in range(param_size):
#    for j in range(param_size):
#        field_array[i,j] = fieldz_param(2e-2, rads[i], lens[j], 20)
#        length_array[i,j] = wire_len(rads[i], lens[j], 20)
#
#good_params = (field_array>0.2)&(length_array<350)