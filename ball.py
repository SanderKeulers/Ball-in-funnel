# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 07:28:38 2023

@author: Sander.Keulers
"""

from matplotlib.widgets import Slider
import matplotlib as mpl
import matplotlib.pyplot as plt
import math as M
import numpy as np

def Funnel(omega, alpha, mu):
    
    """
    Main function, dependent on initial angular velocity (omega) and half angle of the funnel (alpha)
    """
    
    ## Constants
    
    g       = 9.81                  # Gravitational acceleration [m s**-2]
    m       = 0.05                  # Mass of ball [kg]
    R       = 0.01                  # Radius of ball [m]
    I       = (2/5)*m*R**2          # Moment of inertia[kg m**2]

    c       = mu*m*g
    
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
    
    ## Initial conditions
    
    r       = 0.5           # Initial radius of cone of the position of marble
    h       = 3             # Initial height 
    dr_dt   = 0             # Initial vertical velocity
    theta   = 0             # Initial angular position
    dtheta_dt = omega      # Initial angular velocity
    
    ## Time definition and creating empty lists
    
    t      = 0
    dt     = 0.001
    
    height = []         # Height of ball = h
    x      = []         # x-position of ball = cos(alpha)*h
    y      = []         # y-position of ball = cos(alpha)*h
    
    time   = []
    
    
    ## Main while loop following fourth-order Runga-Kutta scheme
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
            
    return x,y,height 
            
            
#%%


fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2, top=0.75)


ax_alpha = fig.add_axes([0.3, 0.85, 0.4, 0.05])
ax_alpha.spines['top'].set_visible(True)
ax_alpha.spines['right'].set_visible(True)

ax_omega = fig.add_axes([0.3, 0.92, 0.4, 0.05])
ax_omega.spines['top'].set_visible(True)
ax_omega.spines['right'].set_visible(True)

ax_mu = fig.add_axes([0.3, 0.78, 0.4, 0.05])
ax_mu.spines['top'].set_visible(True)
ax_mu.spines['right'].set_visible(True)


slider_alpha = Slider(ax=ax_alpha,valmin=1/np.sqrt(8), valmax=1/np.sqrt(2), label='alpha',
              valfmt=' %1.1f ')

slider_omega = Slider(ax=ax_omega,valmin=M.pi, valmax=14*M.pi, label='omega',
              valfmt=' %1.1f ')


slider_mu = Slider(ax=ax_mu,valmin=0.0, valmax=0.6, label='mu',
              valfmt=' %0.01f ')


x,y,height = Funnel(8*M.pi,1/np.sqrt(4),0.0)



f_funnel, = ax.plot(x,height)
ax.set_xlim(min(x)-np.std(x),max(x)+np.std(x))
ax.set_ylim(min(height)-np.std(height),max(height)+np.std(height))
ax.set_xlabel('x-axis [m]')
ax.set_ylabel('Height [m]')


def update(val):
    omega = slider_omega.val
    alpha = slider_alpha.val
    mu = slider_mu.val
    x,y,height = Funnel(omega,alpha,mu)
    f_funnel.set_data(x, height)
    ax.set_xlim(min(x)-np.std(x),max(x)+np.std(x))
    ax.set_ylim(min(height)-np.std(height),max(height)+np.std(height))
    fig.canvas.draw_idle()
    
slider_alpha.on_changed(update)
slider_omega.on_changed(update)
slider_mu.on_changed(update)







