from math import pi
""" 
Input parameters for simulating both Linear and
weakly (as well as fully) nonlinear directional 
water waves. """    

# Initial Discretization of the grid is: 2^npw1 x 2^npw1 
npw1 = 7
npw2 = 7

# Field dimensions (b=beginning, e=end)
x1b = 0.0
x1e = 20.0
y1b = 0.0
y1e = 0.25

# Order of nonlinearity
eqorder = 3
# NOTE: set eqorder to 1 if you want directional LINEAR waves

# Gravity
gv = 9.81

# Water depth
depth = 1.0#5.0

# Significant wave height
# ... todo: add integral formulations for slope
Hs = 0.0326

# Peak frequency
fp = 1.7

# Peak enhancement
gamma = 3.3

# Degree of directional distribution/spreading (nspread~1000, then nearly unidirectional)
nspread = 2.0
# (if =1 then thetam in range [0,pi/2],  or if =2,4,8,16,32,... then keep thetam=pi/2)

# Angular Spread  
thetam = pi/2.0

# Initial filter (keep it =0 as default, or =1 to set as Trulsent et al. for mNLS)
filter = 0
# .... disregard for now and keep it =0
