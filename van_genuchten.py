#!/usr/bin/env python

"""
van Genuchten model to estimate soil water content

A closed-form equation for predicting the hydraulic conductivity of
unsaturated soils. 1980, Soil Sci. Soc. Am., 44, 892-898

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (13.08.2012)"
__email__ = "mdekauwe@gmail.com"

import numpy as np
import matplotlib.pyplot as plt

def van_genuchten_func(alpha, n, h, theta_fc, theta_wp):
    """
    Soil water content (theta) as a function of the pressure head
    """
    m = 1.0 - 1.0 / n
    theta = theta_r + (theta_s - theta_r) / (1.0 + (alpha * np.fabs(h))**n)**m
    theta = np.where(h>=0.0, theta_s, theta)
    beta = (theta - theta_s) / (theta_s - theta_r)

    return (theta, beta)

if __name__ == "__main__":

    # random example, Silt Loam G.E. 3, Table 1
    cm_2_m = 100.0
    h = np.linspace(-5, 2) * cm_2_m # cm-1
    theta_s = 0.396 # saturated water content (cm3 cm-3)
    theta_r = 0.131 # residual water content (cm3 cm-3)
    alpha = 0.00423 # approximates the inverse of the air-entry value, (m-1)
    n = 2.06        # shape param, measure of the pore-size distribution

    (theta, beta) = van_genuchten_func(alpha, n, h, theta_s, theta_r)

    plt.subplot(211)
    plt.plot(h/100., theta, "k-")
    plt.ylabel("Volumetric Water Content (-)")
    plt.xlim(2, -5)

    plt.subplot(212)
    plt.plot(h/100., beta, "k-")
    plt.ylabel("Beta")
    plt.xlabel("Pressure Head, h (m)")
    plt.xlim(2, -5)
    plt.show()
