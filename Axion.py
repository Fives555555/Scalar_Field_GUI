
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import matplotlib.widgets as widgets
from module import *

plt.close("all")

# Parameters
G = 6.67e-11
omega_m = 0.3
omega_r = 8.5e-5
omega_lambda = 1 - omega_m
H0 = 0.7/3000
mu = 6.65
f = 0.1

# geomspace returns numbers spaced evenly on log scale
eta = np.geomspace(1, 1.8e4, 10000)

#Initial conditions
a0 = H0 * np.sqrt(omega_r) * eta[0]
psi0 = 0
phi0 = 0.2075 # decoupled model
y0 = np.array([a0, psi0, phi0, a0])


# Change values according to those selected by the user
def sliderCallback(val):
    
    # Store the values chosen by the user
    line_a = radioHandle.value_selected
    line_lcdm = radioHandlel.value_selected
    line_axion = radioHandleax.value_selected
    omega_m = sliderHandleO_m.val
    omega_lambda = 1 - omega_m
    omega_r = sliderHandlerad.val
    H0 = sliderHandleh.val/3000
    G = sliderHandleG.val * 1e-11
    mu = sliderHandlemu.val
    f = sliderHandlef.val
    phi0 = sliderHandlep.val
    y0 = np.array([a0, psi0, phi0, a0])
    
    # Find the solutions using the new user chosen values
    a = solveLCDM(modelLCDM, eta, 0, 1.8e4, a0, omega_m, omega_r, 
                  omega_lambda, H0)
    a_analytic = analytical(eta, omega_m, omega_r, H0, G, eta)
    solutions = solveAxion(modelAxion, eta, 1, 1.8e4, y0, mu, f, omega_m, 
                           omega_r, omega_lambda, H0)
    xs, cs_scalar, cs_lcdm = HubbleRate(solutions, mu, f, H0, omega_m, 
                                        omega_r, omega_lambda)
    m, r, lam = LCDMEnergy(solutions, omega_m, omega_r, omega_lambda)
    
    V_array = mu * (1 - np.cos(solutions.y[2]/f))
    
    m1, r1, lam1, scal = AxionEnergy(solutions, V_array, H0, omega_m, omega_r, 
                                     omega_lambda)
    
    # Reset each plot's data using the new solutions
    # ax1
    axion.set_xdata(solutions.t)
    axion.set_ydata(solutions.y[0])
    numerical.set_ydata(a.y[0]) 
    analytic.set_ydata(a_analytic)
    
    # ax2
    v.set_xdata(solutions.t)
    v.set_ydata(solutions.y[1])
    
    # ax3
    p.set_xdata(solutions.y[0])
    p.set_ydata(solutions.y[2])
    
    # ax4
    H.set_xdata(1/xs - 1)
    H.set_ydata(cs_scalar(xs)/cs_lcdm(xs))
    
    # ax5
    Om.set_xdata(1/solutions.y[3]- 1)
    Om.set_ydata(m)
    Or.set_xdata(1/solutions.y[3]- 1)
    Or.set_ydata(r)
    OL.set_xdata(1/solutions.y[3]- 1)
    OL.set_ydata(lam)
      
    # ax6
    Om1.set_xdata(1/solutions.y[0] - 1)
    Om1.set_ydata(m1)
    Or1.set_xdata(1/solutions.y[0] - 1)
    Or1.set_ydata(r1)
    OL1.set_xdata(1/solutions.y[0] - 1)
    OL1.set_ydata(lam1)
    Scal.set_xdata(1/solutions.y[0] - 1)
    Scal.set_ydata(scal)    
    
    # Show only the chosen plots for scale factor and conformal time plot
    if line_a == 'All':
        axion.set_visible(True)
        numerical.set_visible(True)
        analytic.set_visible(True)
     
    elif line_a =='Axion + LCDM':
        axion.set_visible(True)
        numerical.set_visible(False)
        analytic.set_visible(False)
        
    elif line_a == 'Numerical LCDM':
        axion.set_visible(False)
        numerical.set_visible(True)
        analytic.set_visible(False) 
        
    elif line_a == 'Analytical LCDM':
        axion.set_visible(False)
        numerical.set_visible(False)
        analytic.set_visible(True)           
    
    # Show only chosen plots for LCDM energy densities
    if line_lcdm == 'All':
        Om.set_visible(True)
        Or.set_visible(True)
        OL.set_visible(True)
    elif line_lcdm == r'$\Omega_{m}$':
        Om.set_visible(True)
        Or.set_visible(False)
        OL.set_visible(False)
    elif line_lcdm == r'$\Omega_{r}$':
        Om.set_visible(False)
        Or.set_visible(True)
        OL.set_visible(False)
    elif line_lcdm == r'$\Omega_{\Lambda}$':
        Om.set_visible(False)
        Or.set_visible(False)
        OL.set_visible(True)
    
    # SHow only chosen plots for LCDM + single axion field energy densities
    if line_axion == 'All':
        Om1.set_visible(True)
        Or1.set_visible(True)
        OL1.set_visible(True)
        Scal.set_visible(True)
    elif line_axion == r'$\Omega_{m,\phi}$':
        Om1.set_visible(True)
        Or1.set_visible(False)
        OL1.set_visible(False)
        Scal.set_visible(False)
    elif line_axion == r'$\Omega_{r,\phi}$':
        Om1.set_visible(False)
        Or1.set_visible(True)
        OL1.set_visible(False)
        Scal.set_visible(False)
    elif line_axion == r'$\Omega_{\Lambda, \phi}$':
        Om1.set_visible(False)
        Or1.set_visible(False)
        OL1.set_visible(True)
        Scal.set_visible(False)
    elif line_axion == r'$\Omega_{\phi}$':
        Om1.set_visible(False)
        Or1.set_visible(False)
        OL1.set_visible(False)
        Scal.set_visible(True)
    
    # Change the value of omega lambda in plot 1 title
    ax1.set_title(r'$\Omega_{\Lambda}$ = ' + str(round(omega_lambda, 4)))    
    # Redraw the axes 
    plt.draw()  


# Solve single axion + LCDM cosmology
solutions = solveAxion(modelAxion, eta, 1, 1.8e4, y0, mu, f, omega_m, omega_r, 
                       omega_lambda, H0)

# solve numerical LCDM and analytic LCDM to plot on first figure
a = solveLCDM(modelLCDM, eta, 0, 1.8e4, a0, omega_m, omega_r, omega_lambda, H0)
a_analytic = analytical(eta, omega_m, omega_r, H0, G, eta)


# Create figure window
fig = plt.figure('', figsize=(17, 9))

# Plot axion + LCDM cosmology scale factor solution against 
# conformal time with numerical LCDM and analytic LCDM
ax1 = fig.add_axes([0.07, 0.67, 0.2, 0.3])
axion, = ax1.loglog(solutions.t, solutions.y[0], label='Axion')
numerical, = ax1.loglog(eta, a.y[0], label="Numerical")
analytic, = ax1.loglog(eta, a_analytic, 
                       label="Analytical, $\Omega_{\Lambda} = 0$")
ax1.set_title(r'$\Omega_{\Lambda}$ = ' + str(omega_lambda))
ax1.set_xlabel(r'$\eta$')
ax1.set_ylabel(r'$a$')
ax1.legend()

# Velocity of axion scalar field evolution against conformal time
ax2 = fig.add_axes([0.33, 0.3, 0.2, 0.3])
v, = ax2.semilogx(solutions.t, solutions.y[1])
ax2.set_xlabel(r'$\eta$')
ax2.set_ylabel(r'$\psi$')
ax2.tick_params(axis='y', rotation=60)

# Axion scalar field phi evolution against scale factor
ax3 = fig.add_axes([0.33, 0.67, 0.2, 0.3])
p, = ax3.semilogx(solutions.y[0], solutions.y[2])
ax3.set_xlabel(r'$a$')
ax3.set_ylabel(r'$\phi$')

# Plot of Hubble rate for scalar field + LCDM compared to LCDM value
xs, cs_scalar, cs_lcdm = HubbleRate(solutions, mu, f, H0, omega_m, omega_r, 
                                    omega_lambda)
ax4 = fig.add_axes([0.07, 0.3, 0.2, 0.3])
H, = ax4.semilogx(1/xs - 1, cs_scalar(xs)/cs_lcdm(xs), label="Cubic Spline")
ax4.set_xlabel(r'$z$')
ax4.set_ylabel(r'$H/H_{\Lambda CDM}$')
ax4.grid()
ax4.grid(which='minor')
ax4.invert_xaxis()
ax4.set_ylim([0.99, 1.025])

# Plot of LCDM energy densities
m, r, lam = LCDMEnergy(solutions, omega_m, omega_r, omega_lambda)
ax5 = fig.add_axes([0.58, 0.67, 0.2, 0.3])
Om, = ax5.semilogx(1/solutions.y[3]- 1, m, label=r'$\Omega_{m}$')
Or, = ax5.semilogx(1/solutions.y[3]- 1, r, label=r'$\Omega_{r}$')
OL, = ax5.semilogx(1/solutions.y[3]- 1, lam, label=r'$\Omega_{\Lambda}$')
ax5.set_xlabel('z')
ax5.legend(loc="center left")
ax5.grid()
ax5.grid(which='minor')
ax5.set_xlim([10e-3, 10e5])
ax5.invert_xaxis()

# Plot of axion + LCDM energy densities
V_array = mu * (1 - np.cos(solutions.y[2]/f))
m1, r1, lam1, scal = AxionEnergy(solutions, V_array, H0, omega_m, omega_r, 
                                 omega_lambda)
ax6 = fig.add_axes([0.58, 0.3, 0.2, 0.3])
Om1, = ax6.semilogx(1/solutions.y[0] - 1, m1, label=r'$\Omega_{m,\phi}$')
Or1, = ax6.semilogx(1/solutions.y[0] - 1, r1, label=r'$\Omega_{r,\phi}$')
OL1, = ax6.semilogx(1/solutions.y[0] - 1, lam1, 
                    label=r'$\Omega_{\Lambda, \phi}$')
Scal, = ax6.semilogx(1/solutions.y[0] - 1, scal, label=r'$\Omega_{\phi}$')
ax6.set_xlabel('z')
ax6.legend(loc="center left", fontsize=9)
ax6.grid()
ax6.grid(which='minor')
ax6.set_xlim([10e-3, 10e5])
ax6.invert_xaxis()


# Add the slider to control the omega_m (and so omega_lambda) value. 
sax = plt.axes([0.12, 0.2, 0.2, 0.03]) 
# Set slider values from 0.0 to 0.34 (the behaviour is too exponential 
# after this), and set initial value of slider to omega_m
sliderHandleO_m = widgets.Slider(sax, r'$\Omega_{m}$', 0.0001, 0.34, 
                                 valinit=omega_m)
sliderHandleO_m.on_changed(sliderCallback)

# Add slider to control the value of omega_r
radax = plt.axes([0.12, 0.15, 0.2, 0.03])
sliderHandlerad = widgets.Slider(radax, r'$\Omega_{r}$', 1e-5, 8e-4, 
                                 valinit=omega_r)
sliderHandlerad.on_changed(sliderCallback)

# Add slider to control the value of H0
hax = plt.axes([0.12, 0.1, 0.2, 0.03])
sliderHandleh = widgets.Slider(hax, r'$h_{0}$', 0.1, 0.7, valinit=0.7)
sliderHandleh.on_changed(sliderCallback)

# Add slider to control value of G
gax = plt.axes([0.12, 0.05, 0.2, 0.03])
sliderHandleG = widgets.Slider(gax, r'$G (x10^{-11})$', 0.1, 10, valinit=6.67)
sliderHandleG.on_changed(sliderCallback)

# Add slider to control value of mu
muax = plt.axes([0.46, 0.2, 0.2, 0.03])
sliderHandlemu = widgets.Slider(muax, r'$\mu$', 0, 20.0, valinit=mu)
sliderHandlemu.on_changed(sliderCallback)

# Add slider to control value of f
fax = plt.axes([0.46, 0.15, 0.2, 0.03])
sliderHandlef = widgets.Slider(fax, r'$f$', 0.0, 2.0, valinit=f)
sliderHandlef.on_changed(sliderCallback)

# Add slider to control value of phi0
pax = plt.axes([0.46, 0.1, 0.2, 0.03])
sliderHandlep = widgets.Slider(pax, r'$\phi_{0}$', 0.005, 0.3, valinit=phi0)
sliderHandlep.on_changed(sliderCallback)


# Add radio buttons to change shown scale factor solutions for the models
rax = plt.axes([0.82, 0.8, 0.12, 0.15])
radioHandle = widgets.RadioButtons(rax, ('All', 
                                         'Axion + LCDM', 
                                         'Numerical LCDM', 
                                         'Analytical LCDM'))
radioHandle.on_clicked(sliderCallback)

# Add radio buttons to change shown LCDM energy density curves 
lax = plt.axes([0.82, 0.54, 0.12, 0.15])
radioHandlel = widgets.RadioButtons(lax, ('All', 
                                          r'$\Omega_{m}$', 
                                          r'$\Omega_{r}$', 
                                          r'$\Omega_{\Lambda}$'))
radioHandlel.on_clicked(sliderCallback)

# Add radio buttons to change shown LCDM + axion energy density curves 
axax = plt.axes([0.82, 0.28, 0.12, 0.15])
radioHandleax = widgets.RadioButtons(axax, ('All', r'$\Omega_{m,\phi}$', 
                                            r'$\Omega_{r,\phi}$', 
                                            r'$\Omega_{\Lambda, \phi}$',
                                            r'$\Omega_{\phi}$'))
radioHandleax.on_clicked(sliderCallback)


# Add button to close figure window
bax = plt.axes([0.87, 0.07, 0.07, 0.05]) # Add new axes to the figure
buttonHandle = widgets.Button(bax, 'Close')
buttonHandle.on_clicked(closeCallback)


plt.show()


