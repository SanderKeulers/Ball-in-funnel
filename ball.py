#-*- coding: utf-8 -*-
"""
Created on Wed Oct  6 22:29:43 2021

@author: sande
"""

import numpy as np
import math as M
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def Trechter(omega):
    #%% Constants
    alpha   = 1/np.sqrt(4)          # Half angle of the funnel [rad]
   
    g       = 9.81                  # Gravitational acceleration [m s**-2]
    m       = 0.01                  # Mass of ball [kg]
    R       = 0.01                  # Radius of ball [m]
    I       = (2/5)*m*R**2          # Moment of inertia[kg m**2]
    c       = 0
   
    #%% Functions
    def d2theta_dt(t,r,dr_dt,dtheta_dt):
        """
        Returns the change in angular velocity as function of angular velocity,
        change in radius/height and radius/height
        """
        return M.exp(-c*t)*((7/5)*-2*dtheta_dt*dr_dt/r)
   
    def d2dr_dt(t,r,dtheta_dt,dr_dt):
        """
        Returns the change in vertical velocity as a function of the angular velocity
        and the radius/height
        """
        return M.exp(c*t)*(((7/5)*r*dtheta_dt**2 * np.sin(alpha)**2 - g*np.cos(alpha)*np.sin(alpha)))
   
   
    #%% Initial conditions
    r       = 0.5           # Initial radius of cone of the position of marble
    h       = 2*r/M.tan(alpha) # Initial height corresponding to initial radius
    dr_dt   = 0             # Initial vertical velocity
    theta   = 0             # Initial angular position
    dtheta_dt = omega      # Initial angular velocity
   
    #%% Time definition and creating emtpy lists
    t      = 0
    dt     = 0.0001
   
    height = []         # Height of ball = h
    x      = []         # x-position of ball = cos(alpha)*h
    y      = []         # y-position of ball = cos(alpha)*h
   
    time   = []
   
    while h>0.03:
        t = t+dt
        time.append(t)
        k1 = dt*d2theta_dt(t,r,dr_dt,dtheta_dt)
        l1 = dt*d2dr_dt(t,r,dr_dt,dtheta_dt)
       
        k2 = dt*d2theta_dt(t+0.5*dt,r+0.5*dt,dr_dt+k1/2,dtheta_dt+l1/2)
        l2 = dt*d2dr_dt(t+0.5*dt,r+0.5*dt,dr_dt+k1/2,dtheta_dt+l1/2)
       
        k3 = dt*d2theta_dt(t+0.5*dt,r+0.5*dt,dr_dt+k2/2,dtheta_dt+l2/2)
        l3 = dt*d2dr_dt(t+0.5*dt,r+0.5*dt,dr_dt+k2/2,dtheta_dt+l2/2)
       
        k4 = dt*d2theta_dt(t+0.5*dt,r+0.5*dt,dr_dt+k3,dtheta_dt+l3)
        l4 = dt*d2dr_dt(t+0.5*dt,r+0.5*dt,dr_dt+k3,dtheta_dt+l3)
       
        k = 1/6 * (k1+2*k2+2*k3+k4)
        l = 1/6 * (l1+2*l2+2*l3+l4)
         
        dtheta_dt = dtheta_dt + k
        dr_dt = dr_dt + l
       
        theta = theta + dtheta_dt*dt
        r = r + dr_dt * dt
       
        if theta > 2*M.pi:
            theta = 0 
           
        h = r/M.sin(alpha)
        height.append(r/M.sin(alpha))
        x.append(M.cos(theta)*r)
        y.append(M.sin(theta)*r)
       
        if h<0.03:
            alpha = M.pi
   

    #%%
   
   
   
    ax.plot(x,y,height,label='$\omega = $'+str(omega))
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(0,2)
    plt.legend()

w = np.arange(1*M.pi, 8*M.pi, M.pi)

for i in w:
    Trechter(i)
