
# Functions used for LCDM and Axion files

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import matplotlib.widgets as widgets 


def modelLCDM(t, a, omega_m, omega_r, H0, omega_lambda=0.0): 
    """Differential equation for LCDM cosmology (k=0) to be solved 
    numerically"""
    dadeta = H0 * (omega_m * a + omega_r + omega_lambda * a**4)**0.5
    
    return dadeta
 
    
def analytical(t, omega_m, omega_r, H0, G, eta):
    """Analytic solution for LCDM, where k=0 and lambda=0"""
    a_eq = omega_r/omega_m
    rho_c = (3 * H0**2)/(8 * np.pi * G)
    rho_eq = 2 * omega_m * rho_c * a_eq**(-3)
    eta_star = ((np.pi * G * rho_eq * a_eq**2)/3)**(-0.5)
    a_analytic = a_eq * ((eta/eta_star)**2 + 2 * (eta/eta_star))
    
    return a_analytic


def solveLCDM(model, t_array, t_initial, t_final, y_initial, omega_m, omega_r, 
              omega_lambda, H0):
    """Function to make calling solve_ivp easier to understand"""
    solutions = solve_ivp(model, [t_initial, t_final], 
                          [y_initial, y_initial], t_eval=t_array, 
                          args=(omega_m, omega_r, H0, omega_lambda))
    
    return solutions


def modelAxion(t, y, mu, f, omega_m, omega_r, omega_lambda, H0):
    """Differential equations for single axion + LCDM cosmology to be solved"""
    dphi_deta = y[1]   
    V = mu * (1 - np.cos(y[2]/f))
    dV_dphi = (mu/f) * np.sin(y[2]/f)
    da_deta = np.sqrt((1/3) * ((V * (y[0]**4) * (H0**2)) + (((y[0]**2)/2) \
        * (dphi_deta)**2)) + (omega_m * y[0] + omega_r + omega_lambda * \
                             y[0]**4) * H0**2 )        
    dpsi_deta = ((-2/y[0]) * da_deta * y[1]) - ((y[0]**2) * (H0**2) * dV_dphi)  
    da_deta_lcdm = H0 * (omega_m * y[3] + omega_r + omega_lambda * y[3]**4)**0.5
   
    return np.array([da_deta, dpsi_deta, dphi_deta, da_deta_lcdm])


def solveAxion(model, t_array, t_initial, t_final, y_initial, mu, f, 
               omega_m, omega_r, omega_lambda, H0):
    """Makes calling solutions with solve_ivp easier"""
    solutions = solve_ivp(model, [t_initial, t_final], y_initial, 
                          t_eval=t_array, args=(mu, f, omega_m, 
                                                omega_r, omega_lambda, H0), 
                          method='LSODA')
    
    return solutions


def HubbleRate(solutions, mu, f, H0, omega_m, omega_r, omega_lambda):
    """Function for calculating specific range of values of 
    Hubble rate to be plotted"""
    V_array = mu * (1 - np.cos(solutions.y[2]/f))

    H_scalar_field = (np.sqrt((1/3) * ((V_array * (solutions.y[0]**4) * \
                                        (H0**2)) + (((solutions.y[0]**2)/2) \
        * (solutions.y[1])**2)) +  (omega_m * solutions.y[0] + omega_r + \
                                    omega_lambda * solutions.y[0]**4) * \
                                    H0**2))/solutions.y[0] 
        
    H_LCDM = (H0 * (omega_m * solutions.y[3] + omega_r + omega_lambda * \
                    solutions.y[3]**4)**0.5)/solutions.y[3]

    cs_scalar = CubicSpline(solutions.y[0], H_scalar_field)
    cs_lcdm = CubicSpline(solutions.y[3], H_LCDM)
    xs = np.arange(0.01, 1, 0.00001) # scale factor in chosen range
    
    return xs, cs_scalar, cs_lcdm


def LCDMEnergy(solutions, omega_m, omega_r, omega_lambda):
    """Calculates the energy densities for LCDM cosmology"""
    LCDM_densities = omega_m * solutions.y[3] + omega_r + \
        omega_lambda * solutions.y[3]**4
    m = (omega_m * solutions.y[3]) / (LCDM_densities)
    r = (omega_r) / (LCDM_densities)
    lam = (omega_lambda * solutions.y[3]**4) / (LCDM_densities)
    
    return m, r, lam


def AxionEnergy(solutions, V_array, H0 , omega_m, omega_r, omega_lambda):
    """Calculates energy densities for single axion + LCDM cosmology"""
    LCDM_s_densities = omega_m * solutions.y[0] + omega_r + \
        omega_lambda * solutions.y[0]**4  
    s_density = (1/3) * (V_array * (solutions.y[0]**4) + 1/(H0**2) * \
                         (((solutions.y[0]**2)/2) * (solutions.y[1])**2))                                                                             
    m1 = (omega_m * solutions.y[0]) / (s_density + LCDM_s_densities)
    r1 = (omega_r) / (s_density + LCDM_s_densities)
    lam1 = (omega_lambda * solutions.y[0]**4) / (s_density + LCDM_s_densities)
    scal = (s_density) / (s_density + LCDM_s_densities)
    
    return m1, r1, lam1, scal


def closeCallback(event):
    """Close all open figure windows from button press"""
    plt.close('all') 


