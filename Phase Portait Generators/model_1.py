# -*- coding: utf-8 -*-

# Plots phase portraits and nullclines for model 2
# -----------------------------------------------------------------------------


#  ----------------------------------------------------------------------------
#  MODULES
#  ----------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.optimize import root_scalar
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


# -----------------------------------------------------------------------------
# MATHEMATICAL FUNCTIONS
# -----------------------------------------------------------------------------
def DR(R):
    """
    Time derivative of R
    """
    return r1*R*(K1**B1/(K1**B1 + (A/R)**B1)) \
    - gwt*A - r2*R*(K2**B2/(K2**B2 + (A)**B2)) + R*gc

vDR = np.vectorize(DR)

def DPc(Pc):
    """
    Time derivative of Pc
    """
    return r2*(K2**B2/(K2**B2 + (A)**B2))*Pc - gc*Pc

vDPc = np.vectorize(DPc)

def mag(x,y):
    return sqrt(x**2 + y**2)

vmag = np.vectorize(mag)

def PcNull(R):
    """
    Explicit definition of the \dot{P_C}=0 nullcline as a function of R
    """
    pass

vPcNull = np.vectorize(PcNull)

def RNull(R):
    """
    Explicit definition of the \dot{R}=0 nullcline as a function of R
    """
    pass

vRNull = np.vectorize(RNull)


# -----------------------------------------------------------------------------
# PLOT SETUP
# -----------------------------------------------------------------------------
small_float = 0.0000001
granularity = 30 # Mesh spacing

# Axis limits
R_max = 20
Pc_max = 30

# Nullcline colormap
cmap = plt.get_cmap("tab10")

# Grid for plotting nullclines
Rrange = np.linspace(small_float, R_max, 10000, endpoint=False)

# Mesh for vector field
Rsteps,Pcsteps = np.meshgrid(
    np.linspace(small_float, R_max, granularity, endpoint=False),
    np.linspace(small_float, Pc_max, granularity, endpoint=False))

# Vector components
U = vDR(Rsteps)
V = vDPc(Pcsteps)

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
plt.plot([0, R_max], [0, 0], label=r'$\dot{P_C} = 0$', color=cmap(1))
if A > 0:
    sol = root_scalar(DR, method='bisect', bracket=[small_float, 20],
        x0=small_float)
    R_root = sol.root
    plt.plot([R_root, R_root],[0, Pc_max], label=r'$\dot{R} = 0$', color=cmap(0))
else:
    plt.plot(0, 0, label=r'$\dot{R} = 0$', color=cmap(0))

# Set plot style
plt.grid(linestyle='-', linewidth=0.5)
plt.xlim(0, R_max)
plt.ylim(-1, Pc_max)
plt.legend()
plt.xlabel(r'$R$')
plt.ylabel(r'$P_C$')
plt.show()
