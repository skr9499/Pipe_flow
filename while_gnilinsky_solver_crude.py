#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 09:07:59 2023

@author: sreekrishnaravishankar
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as cp

#=====================================================================================================#
def colebrook_equation(Re, epsilon, D):
    f = 0.05                                                                    # Initial guess for friction factor
    error = 1                                                                   # Error threshold
    max_iterations = 10
    iterations = 0

    while (error > 1e-6) and (iterations < max_iterations):
        f_new = (-2 * math.log10((epsilon / (3.7 * D)) + (2.51 / (Re * math.sqrt(f))))) ** -2
        error = abs(f_new - f)
        f = f_new
        iterations = iterations + 1    
    friction_factor = f
    return friction_factor

def Reynolds_number(G, dx, viscosity):
    Re = (G * dx) / viscosity
    return Re

def Cp_fluid(P_x, rho_x):
    Cp = cp('C', 'P', P_x, 'D', rho_x, 'Hydrogen')
    return Cp

def Thermal_conducitvity(P_x, rho_x):
    k = cp('L', 'P', P_x, 'D', rho_x, 'Hydrogen')
    return k

def Viscosity(P_x, rho_x):
    Viscosity = cp('V', 'P', P_x, 'D', rho_x, 'Hydrogen')
    return Viscosity

def Prandtl_number(cp, mu, k):
    Prandtl_number = (cp * mu) / k
    return Prandtl_number

def Gnielinski(f_x, Re_x, Pr):
    Nu = ((f_x / 8) * (Re_x - 1000) * Pr) / (1 + 12.7 * ((f_x / 8) ** 0.5) * ((Pr ** (2 / 3)) - 1))
    return Nu

def Momentum_Eqn(P_x, rho_x, G, f_x, D, dx, P_new):
    v_x = 1 / rho_x
    spec_vol_new = (P_x - P_new + ((G ** 2) * v_x) * (1 - ((f_x * dx) / (2 * D)))) / (G ** 2)
    rho_new = 1 / spec_vol_new
    return rho_new

def Eqn_of_State(P_new, rho_new):
    T_new = cp('T', 'P', P_new, 'D', rho_new, 'Hydrogen')
    return T_new

T_surface = 1200
P_x = 45e5
T_x = 25
samples = 10000
mdot = 37.962
D = 1e-2
Epsilon = 0.3853e-6
A_pipe = math.pi * (D ** 2) / 4
G = mdot / A_pipe
rho_x = cp('D', 'P', P_x, 'T', T_x, 'Hydrogen')
vel_x = G / rho_x
dx = 0.01
L_pipe = 0.5
pipe_elements = L_pipe / dx
segment = 0.00
segment_array = np.zeros(int(pipe_elements))
P_x_array = np.zeros(int(pipe_elements))
T_x_array = np.zeros(int(pipe_elements))
vel_x_array = np.zeros(int(pipe_elements))
rho_x_array = np.zeros(int(pipe_elements))
error_final_array = np.zeros(int(pipe_elements))
current_index = 0

while (segment < pipe_elements) and (current_index<pipe_elements):
    P_start = P_x - 1e5
    P_new = np.linspace(P_x, P_start, samples)
    rho_new = np.zeros(samples)
    T_new = np.zeros(samples)
    error = np.zeros(samples)

    Cp = Cp_fluid(P_x, rho_x)
    k = Thermal_conducitvity(P_x, rho_x)
    mu = Viscosity(P_x, rho_x)
    Pr = Prandtl_number(Cp, mu, k)
    Re_x = Reynolds_number(G, dx, mu)
    f_x = colebrook_equation(Re_x, Epsilon, D)
    Nu = Gnielinski(f_x, Re_x, Pr)
    h = k * Nu / D
    A_pipe_element = (math.pi * D) * dx
    T_out = T_surface - ((T_surface - T_x) * math.exp(-(h * A_pipe_element) / (mdot * Cp)))

    for i in range(len(P_new)):                                                 #Pnew = keep vel head const
        rho_new[i] = Momentum_Eqn(P_x, rho_x, G, f_x, D, dx, P_new[i])
        T_new[i] = Eqn_of_State(P_new[i], rho_new[i])
        error[i] = abs(T_new[i] - T_out)

    for i in range(len(error)):
        index = np.argmin(error)
    error_final = error[index]
    P_final = P_new[index]
    rho_final = rho_new[index]
    vel_final = G / rho_final
        
    segment_array[current_index] = segment
    P_x_array[current_index] = P_x
    T_x_array[current_index] = T_x
    vel_x_array[current_index] = vel_x
    rho_x_array[current_index] = rho_x
    error_final_array[current_index] = error_final

    P_x = P_final
    rho_x = rho_final
    vel_x = vel_final
    T_x = T_out
    segment = segment + dx
    current_index = current_index+1

        
plt.figure(1)

plt.plot(segment_array, P_x_array, label='P_x')
plt.xlabel('Segment')
plt.ylabel('P_x')
plt.legend()

plt.figure(2)
plt.plot(segment_array, T_x_array, label='T_x')
plt.xlabel('Segment')
plt.ylabel('T_x')
plt.legend()

plt.figure(3)
plt.plot(segment_array, vel_x_array, label='vel_x')
plt.xlabel('Segment')
plt.ylabel('vel_x')
plt.legend()

plt.figure(4)
plt.plot(segment_array, rho_x_array, label='rho_x')
plt.xlabel('Segment')
plt.ylabel('rho_x')
plt.legend()


# Plotting the error_final along the segment variable
plt.figure(5)
plt.plot(segment_array, error_final_array, label='error_final') 
plt.xlabel('Segment')
plt.ylabel('Error Final')
plt.grid()
plt.legend()
plt.show()

