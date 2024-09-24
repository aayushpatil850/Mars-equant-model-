#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
import numpy as np
import datetime
from tqdm import tqdm
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize
get_ipython().run_line_magic('matplotlib', 'inline')

# Load data
marsData = pd.read_csv('/home/aayushjeevan/Downloads/01_data_mars_opposition_updated.csv')




# In[8]:


# Extract heliocentric longitude values
marsHeliocentric_longitude = marsData.iloc[:, 5:9].values
print(marsHeliocentric_longitude)

# Convert heliocentric longitude to degrees
marsHeliocentric_longitude_in_degrees = (
    marsData['ZodiacIndex'] * 30 + 
    marsData['Degree'] + 
    marsData['Minute.1'] / 60.0 + 
    marsData['Second'] / 3600.0
)

# Convert degrees to radians
marsHeliocentric_longitude_in_radians = np.deg2rad(marsHeliocentric_longitude_in_degrees)

# Compute time intervals between consecutive observations
times = [0]
for i in range(1, len(marsData)):
    prev_date = datetime.datetime(
        marsData['Year'][i-1], marsData['Month'][i-1],
        marsData['Day'][i-1], marsData['Hour'][i-1],
        marsData['Minute'][i-1]
    )
    
    curr_date = datetime.datetime(
        marsData['Year'][i], marsData['Month'][i],
        marsData['Day'][i], marsData['Hour'][i],
        marsData['Minute'][i]
    )
    
    duration = curr_date - prev_date
    days_elapsed = duration.days + duration.seconds / (60 * 60 * 24)
    times.append(days_elapsed)
times = np.array(times)
#print(times)


# In[9]:


# Stack times and heliocentric longitudes
oppositions = np.column_stack((times, marsHeliocentric_longitude_in_degrees))


# In[10]:


def getIntersectionPoint(h, k, theta, r, c):
    
    cos_theta = np.cos(np.radians(theta))
    sin_theta = np.sin(np.radians(theta))

    cos_c = np.cos(np.radians(c))
    sin_c = np.sin(np.radians(c))

    b = 2 * ((h * cos_theta) + (k * sin_theta) - (cos_c * cos_theta) - (sin_c * sin_theta))
    c1 = h**2 + k**2 + 1 - 2 * (h * cos_c) - 2 * (k * sin_c) - r**2
    l1 = -b / 2

    try:
        l2 = np.sqrt(b**2 - 4 * c1) / 2
    except:
        l2 = 0
        
    root1 = l1 + l2
    root2 = l1 - l2
    
    ell_value = root1 if root1 > 0 else root2
    
    return (h + ell_value * cos_theta), (k + ell_value * sin_theta)

def MarsEquantModel(c, r, e1, e2, z, s, oppositions):
    
    errors = []
    xpos = []
    ypos = []
    
    h = e1 * np.cos(np.radians(e2 + z))
    k = e1 * np.sin(np.radians(e2 + z))

    theta_new = z
    for i in range(12):
        theta = (s * times[i]) + theta_new
        x, y = getIntersectionPoint(h, k, theta, r, c)
        xpos.append(x)
        ypos.append(y)
        angle = np.degrees(np.arctan2(y, x)) % 360
        errors.append(abs(oppositions[i][1] - angle))
        theta_new = theta
    max_error = max(errors)
    return errors, max_error


# In[11]:



errors, maxError = MarsEquantModel(149.0000001, 8.59999999999999, 1.60000000001, 93.19999999, 55.800000001, 0.524093, oppositions)
print("12 spokes errors-:", errors)
print("---------------------------------------------------------")
print('maximum error:-', maxError)

def bestOrbitInnerParams(r, s, oppositions):
    """
    Find the optimal parameters for the Mars orbit model by minimizing the maximum error.
    """
    min_error = 1e20
    for c in tqdm(np.arange(149, 149.5, 0.01)):
        for e2 in np.arange(93, 94, 0.1):
            for z in np.arange(55, 56, 0.1):
                for e1 in np.arange(1.4, 1.6, 0.01):
                    errors, max_error = MarsEquantModel(c, r, e1, e2, z, s, oppositions)
                    if min_error > max_error:
                        min_error = round(max_error, 4)
                        best_c = c
                        best_e1 = round(e1, 4)
                        best_e2 = e2
                        best_z = z
                        best_errors = errors
                        print("C value:-", best_c, " E1 value:-", best_e1, " E2 value:-", best_e2, " z value:-", best_z, " Maximum Error:-", min_error)
    return best_c, best_e1, best_e2, best_z, best_errors, min_error

# Find the best orbit parameters
c, e1, e2, z, errors, maxError = bestOrbitInnerParams(9, 0.524, oppositions)
print("------------------------------------------------------------------------------")
print("C value:-", c, " E1 value:-", e1, " E2 value:-", e2, " z value:-", z, " Error List:-", errors, " Maximum Error:-", maxError)


# In[12]:


def bestS(r, oppositions):
    """
    Determine the optimal value of `s` for the Mars orbit model by minimizing the maximum error.
    """
    least_error = 1e20
    optimal_s = None
    best_errors = None
    
    try:
        for s in np.arange(680, 689, 0.2):
            print(f"Testing S value: {360 / s}")  # Debug print
            best_c, best_e1, best_e2, best_z, errors, max_error = bestOrbitInnerParams(r, 360 / s, oppositions)
            if least_error > max_error:
                optimal_s = 360 / s
                best_errors = errors
                least_error = round(max_error, 4)
                print(f"S value: {optimal_s} Error List: {best_errors} Maximum Error: {max_error}")
    except Exception as e:
        print(f"Error occurred in bestS function: {e}")
    
    return optimal_s, best_errors, least_error

def bestR(s, oppositions):
    """
    Determine the optimal radius `r` for the Mars orbit model by minimizing the maximum error.
    """
    least_error = 1e20
    best_r = None
    best_errors = None
    
    try:
        for r in np.arange(5, 10, 0.1):
            print(f"Testing R value: {r}")  # Debug print
            best_c, best_e1, best_e2, best_z, errors, max_error = bestOrbitInnerParams(r, s, oppositions)
            if least_error > max_error:
                best_r = r
                best_errors = errors
                least_error = round(max_error, 4)
                print(f"R value: {best_r} C value: {best_c} E1 value: {best_e1} E2 value: {best_e2} Z value: {best_z} Maximum Error: {least_error}")
    except Exception as e:
        print(f"Error occurred in bestR function: {e}")
    
    return best_r, best_errors, least_error


s, errors, maxError = bestS(7, oppositions)
print("------------------------------------------------------------------------------")
print(f"S value: {s} Error List: {errors} Maximum Error: {maxError}")

r, errors, maxError = bestR(0.524, oppositions)
print("------------------------------------------------------------------------------")
print(f"Best R value: {r} Error List: {errors} Maximum Error: {maxError}")


# In[13]:


def plotMarsOrbit(c,r,e1,e2,z,oppositions):
    
    
    figure, axes = plt.subplots(1,1,figsize=(10,10))
    #change default range 
#     axes.figsize(5,5)
    axes.set_xlim((-10, 10))
    axes.set_ylim((-10, 10))
    
    
    
    CentreXPos = math.cos(math.radians(c))
    CentreYPos = math.sin(math.radians(c))
    h = e1 * math.cos(math.radians(e2+z))
    k = e1 * math.sin(math.radians(e2+z))
    xpos=[]
    ypos=[]
    
    thetaNew = z
    for i in range(12):
        theta = (s * oppositions[i][0]) + thetaNew
        x,y = getIntersectionPoint(h,k,theta,r,c)
        plt.scatter(x,y,c = 'k')
        plt.plot([h,x],[k,y],linestyle = 'dashed')
        xpos.append(x)
        ypos.append(y)
        thetaNew = theta
        
    orbit = plt.Circle((CentreXPos, CentreYPos), r, color='blue', fill = False)

    
    for i in range(0,12):
        xpos = math.cos(math.radians(oppositions[i][1]))
        ypos = math.sin(math.radians(oppositions[i][1]))
        
        x = [0,xpos*20]
        y = [0,ypos*20]


        plt.plot(x,y)
    

    
    axes.annotate('centre', xy =(CentreXPos, CentreYPos), fontsize = 14)
    axes.annotate('sun', xy =(0, 0),fontsize = 14)
    axes.annotate('equant', xy =(h, k),fontsize = 14)
    
    axes.scatter([0,CentreXPos,h],[0,CentreYPos,k])
    
    
    axes.add_artist(orbit)
    plt.title( 'Mars Predicted Orbit' )
    plt.show()


# In[14]:


plotMarsOrbit(149,8.599999,1.6,93.199999999,55.800000001,oppositions)


# In[15]:


errors, maxError = MarsEquantModel(149.0000001,8.59999999999999,1.60000000001,93.19999999,55.800000001,0.524093,oppositions)


# In[16]:


errors, maxError


# In[17]:


def bestMarsOrbitParams(oppositions):
    
    # Initialize the best error with a very large number
    min_error = 1e24

    # Iterate over possible values of `r`
    for r in tqdm(np.arange(8, 9, 0.1)):
        # Iterate over possible values of `s`
        for s in np.arange(686, 687, 0.1):
            # Call the function to get errors and other parameters
            c, e1, e2, z, error_list, max_error = bestOrbitInnerParams(r, 360 / s, oppositions)
            
            # Check if the current maximum error is less than the best error found
            if min_error > max_error:
                # Update the best parameters if the current configuration is better
                optimum_r = r
                optimum_s = 360 / s
                optimum_c = c
                optimum_e1 = e1
                optimum_e2 = e2
                optimum_z = z
                best_error_list = error_list
                min_error = max_error
    
    return (optimum_r, optimum_s, optimum_c, optimum_e1, optimum_e2, optimum_z, 
            best_error_list, min_error)


# In[19]:


# Ensure `oppositions` is defined in your environment before calling this function
opt_r, opt_s, opt_c, opt_e1, opt_e2, opt_z, err_list, err = bestMarsOrbitParams(oppositions)
print(f"Optimal R: {opt_r}")
print(f"Optimal S: {opt_s}")
print(f"Optimal C: {opt_c}")
print(f"Optimal e1: {opt_e1}")
print(f"Optimal e2: {opt_e2}")
print(f"Optimal Z: {opt_z}")
print(f"Best Error: {err}")
print(f"error list: {err_list}")

