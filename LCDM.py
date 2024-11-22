
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets 
from module import *

# All equations dimensionless

# Constants and initial parameter values
a0 = 0
omega_m = 0.3
omega_r = 9e-5
omega_lambda = 1 - omega_m
H0 = 0.7/3000
G = 6.67e-11
eta = np.linspace(0, 1.8e4, 1000)


# Change values according to those selected by the user
def sliderCallbackLCDM(val):
    
    # Store the values chosen by the user
    line = radioHandle.value_selected
    omega_m = sliderHandleO_m.val
    omega_lambda = 1 - omega_m
    
    # Find the analytic and numerical solutions using the new omega values
    a = solveLCDM(modelLCDM, eta, 0, 1.8e4, a0, omega_m, omega_r, 
                  omega_lambda, H0)
    a_analytic = analytical(eta, omega_m, omega_r, H0, G, eta)
    
    # Reset the y data using the new solutions
    axesHandle.set_ydata(a.y[0]) 
    axesHandle2.set_ydata(a_analytic)
      
    # Show only the chosen plots
    if line == 'Both':
        axesHandle2.set_visible(True)
        axesHandle.set_visible(True) 
     
    elif line == 'Numerical':
        axesHandle2.set_visible(False)
        axesHandle.set_visible(True) 
        
    elif line == 'Analytical':
        axesHandle2.set_visible(True)
        axesHandle.set_visible(False)
           
    # Change the value of omega lambda in the title
    ax.set_title(r'$\Omega_{\Lambda}$ = ' + str(round(omega_lambda, 4)))
    # Redraw the axes 
    plt.draw()  
       
# Find numerical solution initially
a = solveLCDM(modelLCDM, eta, 0, 1.8e4, a0, omega_m, omega_r, omega_lambda, H0)

# Find analytic solution initially
a_analytic = analytical(eta, omega_m, omega_r, H0, G, eta)

# Create figure and axes
fig = plt.figure('Scale factor against conformal time for LCDM cosmology', 
                 figsize=(9, 7))  
ax = plt.axes([0.1, 0.32, 0.5, 0.6]) 

# Plotting analytical and numerical (python) scale factor solutions to the 
# LCDM Freidmann equation against conformal time
axesHandle, = plt.loglog(eta, a.y[0], label="Numerical")
axesHandle2, = plt.loglog(eta, a_analytic, 
                          label="Analytical, $\Omega_{\Lambda} = 0$")
ax.set_title(r'$\Omega_{\Lambda}$ = ' + str(omega_lambda))
ax.set_xlabel(r'$\eta (Mpc)$'), ax.set_ylabel(r'$a$')
plt.legend(loc="upper left")

# Add the slider to control the omega_m (and so omega_lambda) value. 
sax = plt.axes([0.1, 0.15, 0.5, 0.05]) 
# Set slider values from 0.0 to 0.34 (the behaviour is too exponential
# after this), and set initial value of slider to omega_m
sliderHandleO_m = widgets.Slider(sax, r'$\Omega_{m}$', 0.0001, 0.34, 
                                 valinit=omega_m)
sliderHandleO_m.on_changed(sliderCallbackLCDM)

# Add radio buttons to plot both or analytical/numerical solutions only
rax = plt.axes([0.72, 0.54, 0.2, 0.2])
radioHandle = widgets.RadioButtons(rax, ('Both', 
                                         'Numerical', 
                                         'Analytical'))
radioHandle.on_clicked(sliderCallbackLCDM)

# Add button to close figure window
bax = plt.axes([0.82, 0.08, 0.1, 0.05]) # Add new axes to the figure
buttonHandle = widgets.Button(bax, 'Close')
buttonHandle.on_clicked(closeCallback)


plt.show()

