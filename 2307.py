# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 00:36:09 2023

@author: S335830
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI as cp

#===================================Functions==================================

def Cp_fluid(P_x, T_x):
    Cp = cp('C', 'P', P_x, 'T', T_x, 'Hydrogen')
    return Cp

def Thermal_conducitvity(P_x, T_x):
    k = cp('L', 'P', P_x, 'T', T_x, 'Hydrogen')
    return k

def Prandtl_number(cp, mu, k):
    Prandtl_number = (cp * mu) / k
    return Prandtl_number

def Viscosity(P_x, T_x):
    Viscosity = cp('V', 'P', P_x, 'T', T_x, 'Hydrogen')
    return Viscosity

def Reynolds_number(G, D, viscosity):
    Re = (G * D) / viscosity
    return Re

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

def Nusselt_No(f_x, Re_x, Pr):
    Nu = ((f_x / 8) * (Re_x - 1000) * Pr) / (1 + 12.7 * ((f_x / 8) ** 0.5) * ((Pr ** (2 / 3)) - 1))
    return Nu

def Asymptotic_T_increase (Nu, D, k, dx, T_surface, mdot, Cp):
    h = k * Nu / D
    A_pipe_element = (math.pi * D) * dx
    T_out = T_surface - ((T_surface - T_x) * math.exp(-(h * A_pipe_element) / (mdot * Cp)))
    return T_out

def Approx_Pressure_drop(G,dx, rho_x, f, D):
    v_x = 1/rho_x
    drop = 2*(((G**2)*f*dx*v_x)/(2*D))
    return drop

#================================Pipe Dimensions===============================
mdot = 37
D = 5e-2
A_pipe = math.pi * (D ** 2) / 4
G = mdot / A_pipe
Epsilon = 0.3853e-6
dx = 0.01
L_pipe = 0.5  

Pipe_elements = L_pipe / dx
segment = 0.00
Pipe_Length = np.zeros(int(Pipe_elements))
P_field = np.zeros(int(Pipe_elements))
T_field = np.zeros(int(Pipe_elements))
Vel_field = np.zeros(int(Pipe_elements))
Rho_field = np.zeros(int(Pipe_elements))
Error_array = np.zeros(int(Pipe_elements))
current_index = 0
samples = 10000
 
#===========================Pipe Inlet Thermodynamics==========================
P_x = 25e5
T_x = 25
rho_x = cp('D', 'P', P_x, 'T', T_x, 'Hydrogen')
vel_x = G / rho_x                                             
T_surface = 1200

#=================================Computation==================================
mu = Viscosity(P_x, T_x)
Cp = Cp_fluid(P_x, T_x)
k = Thermal_conducitvity(P_x, T_x)
mu = Viscosity(P_x, T_x)
Pr = Prandtl_number(Cp, mu, k)
Re_x = Reynolds_number(G, D, mu)
f_x = colebrook_equation(Re_x, Epsilon, D)
Nu = Nusselt_No(f_x, Re_x, Pr)
T_out = Asymptotic_T_increase(Nu, D, k, dx, T_surface, mdot, Cp)
#=========================Solver Initialisation================================
P_drop_Approx = Approx_Pressure_drop(G, dx, rho_x, f_x, D)
Probable_P_start = P_x - P_drop_Approx
Probable_P = np.linspace(P_x,Probable_P_start,samples)
