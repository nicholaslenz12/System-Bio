# -*- coding: utf-8 -*-

# Plots phase portraits and nullclines for model 2
# -----------------------------------------------------------------------------


#  ----------------------------------------------------------------------------
#  MODULES
#  ----------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from math import log
from math import sqrt


# -----------------------------------------------------------------------------
# GLOBALS/PARAMETERS
# -----------------------------------------------------------------------------
r1 = log(2)/20
r2 = log(2)/40
B1 = 1.5
B2 = 1.5
K1 = .1
K2 = 3
gwt = .01*r1
gc = .001*r2
A = 0
S = 30


# -----------------------------------------------------------------------------
# MATHEMATICAL FUNCTIONS
# -----------------------------------------------------------------------------
def DR(R,Pc):
    """
    Time derivative of R
    """
    return r1*R*(K1**B1/(K1**B1 + (A/R)**B1))*(S/(S + R*Pc + Pc)) \
    - gwt*A - r2*R*(K2**B2/(K2**B2 + (A)**B2))*(S/(S + R*Pc + Pc)) + R*gc

vDR = np.vectorize(DR)

def DPc(R,Pc):
    """
    Time derivative of Pc
    """
    return r2*(K2**B2/(K2**B2 + (A)**B2))*(S/(S + R*Pc + Pc))*Pc - gc*Pc

vDPc = np.vectorize(DPc)

def mag(x,y):
    return sqrt(x**2 + y**2)

vmag = np.vectorize(mag)

def PcNull(R):
    """
    Explicit definition of the \dot{P_C}=0 nullcline as a function of R
    """
    return (r2*(K2**B2/(K2**B2 + A**B2))*S - gc*S)/((1+R)*gc)

vPcNull = np.vectorize(PcNull)

def RNull(R):
    """
    Explicit definition of the \dot{R}=0 nullcline as a function of R
    """
    return -S*R*(r1*(K1**B1/(K1**B1 + (A/R)**B1)) - r2*(K2**B2/(K2**B2 + (A)**B2))) \
        /((-gwt*A + R*gc)*(R+1)) + S/(R+1)

vRNull = np.vectorize(RNull)


# -----------------------------------------------------------------------------
# PLOT SETUP
# -----------------------------------------------------------------------------
small_float = 0.0000001
granularity = 30 # Mesh spacing

# Axis limits
R_max = 400
Pc_max = 600

# Nullcline colormap
cmap = plt.get_cmap("tab10")

# Grid for plotting nullclines
Rrange = np.linspace(small_float, R_max, 10000, endpoint=False)
Rrange_asymptote = np.linspace(small_float, A*gwt/gc, 10000, endpoint=False)

# Mesh for vector field
Rsteps,Pcsteps = np.meshgrid(
    np.linspace(small_float, R_max, granularity, endpoint=False),
    np.linspace(small_float, Pc_max, granularity, endpoint=False))

# Vector components
U = vDR(Rsteps, Pcsteps)
V = vDPc(Rsteps, Pcsteps)

# Length of vector per mesh point
N = np.sqrt(U**2+V**2)
U2, V2 = U/N, V/N # Normalize (Potential divide by zero)



# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
# TeX engine
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Plot vector field
plt.quiver(Rsteps, Pcsteps, U2, V2, width=.002, headwidth=2, headlength=2,
    color='gray')

# Plot nullclines
plt.plot(Rrange, vPcNull(Rrange), label=r'$\dot{R} = 0$', color=cmap(0))
plt.plot(Rrange_asymptote, vRNull(Rrange_asymptote), label=r'$\dot{P_C} = 0$',
    color=cmap(1))
plt.plot([0, R_max],[0, 0], color=cmap(1))

# Set plot style
plt.grid(linestyle='-', linewidth=0.5)
plt.xlim(0, R_max)
plt.ylim(-Pc_max/100, Pc_max)
plt.legend()
plt.xlabel(r'$R$')
plt.ylabel(r'$P_C$')
plt.show()
