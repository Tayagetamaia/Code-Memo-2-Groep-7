# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define physical parameters for the system
k = 84  # Spring constant in N/m   (42)
m = 1.9372e-6  # Mass of the object in kg
y = 0.0235  # Damping coefficient in kg/s     (0.045101) 0.0235

epsilon_0 = 8.854e-12  # Permittivity of free space in F/m
A = 1e-6  # Area of the capacitor plates in m^2
d = 2e-6  # Initial distance between plates in meters
V_s = 1e-6  # Voltage in volts

def calculate_voltage(x):
    C1 = epsilon_0 * A / (d - x)
    C2 = epsilon_0 * A / (d + x)
    V1 = -V_s * C1 / (C1 + C2)
    V2 = V_s * C2 / (C1 + C2)
    return V1 + V2


peak_treshold = 14
rest_treshold = 0.005
search_depth = 100

# Read in the position data from a text file, assuming whitespace delimiter remove or add the # to choose which position file to use
#data = pd.read_csv('Posities\posities_1_Team_07.txt', sep='\s+')    

data = pd.read_csv('Posities\posities_2_Team_07.txt', sep='\s+')   

# Calculate velocity by finite difference of position 'x' with respect to time 't'
# Shift(1) moves the data down one row for subtraction, division by time interval gives velocity
data['v (m/s)'] = (data['x'] - data['x'].shift(1)) / (data['t'] - data['t'].shift(1))

# Calculate acceleration by finite difference of velocity with respect to time
# Resulting NaN values from shift are replaced with 0
data['a (m/s^2)'] = (data['v (m/s)'] - data['v (m/s)'].shift(1)) / (data['t'] - data['t'].shift(1))
data = data.fillna(0)

# Convert time column to a numpy array for numerical operations
t = np.array(data['t'].values)

# Convert the acceleration data for use in numerical integration and initialize arrays
xf = np.array(data['a (m/s^2)'].values)  # External force from acceleration data
x = np.zeros_like(t)  # Array to store position values calculated by Euler's method
x_prime = np.zeros_like(t)  # Array to store velocity values calculated by Euler's method

# Implement Euler's method to solve the differential equation for motion
# This loop iterates over each time step, updating the position 'x' and velocity 'x_prime'
for i in range(len(t)-1):
    # Calculate the net acceleration from force balance: F = ma, considering damping and spring force
    x_2prime = (m*xf[i] - y*x_prime[i] - k * x[i]) / m
    # Update velocity using Euler's forward method: v = v0 + a*dt
    x_prime[i+1] = x_prime[i] + (t[i+1] - t[i]) * x_2prime
    # Update position using Euler's forward method: x = x0 + v*dt
    x[i+1] = x[i] + (t[i+1] - t[i]) * x_prime[i]    

dydx_x = np.array([0])

for i in range(len(x)-1):
    derivative = abs((x[i+1] - x[i]) /  (data['t'][i+1] - data['t'][i]))
    dydx_x = np.append(dydx_x, derivative)

dydx_max = np.argmax(dydx_x)

dydx_x = dydx_x / max(dydx_x)


def is_valid_point(idx, lst, depth):
    # Check if all points within 'depth' after 'idx' are below 'rest_treshold'
    return all(value < rest_treshold for value in lst[idx:idx+depth])

for i in range(len(dydx_x)):
    
    idx_offset = np.argmax(dydx_x) + i
    
    if is_valid_point(idx_offset, dydx_x, search_depth):
        balance_x = idx_offset
        break



for i in range(len(data['a (m/s^2)'])):
    if data['a (m/s^2)'][i] > peak_treshold:
        peak_x = i
        break


real_volt_response = np.array([])

for i in range(len(data['x'])):
    real_volt_response = np.append(real_volt_response, calculate_voltage(data['x'][i]))


integral_velocity = np.array([0])
for i  in range(len(x)-1):
    velocity_integral = (data['t'][i+1] - data['t'][i]) * ((x[i+1] + x[i]) / 2)
    integral_velocity = np.append(integral_velocity, integral_velocity[-1] + velocity_integral)


integral_position = np.array([0])
for i in range(len(integral_velocity)-1):
    position_integral = (data['t'][i+1] - data['t'][i]) * ((integral_velocity[i+1] + integral_velocity[i]) / 2)
    integral_position = np.append(integral_position, integral_position[-1] + position_integral)

accelerometer_volt_response = np.array([])
for i in range(len(integral_position)):
    accelerometer_volt_response = np.append(accelerometer_volt_response, calculate_voltage(integral_position[i]))


dt_peaks = data['t'][balance_x] - data['t'][peak_x]

print(dt_peaks)

# plot both lines in one graph
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()

# ax1.set_ylim(-1, 21.135)
# ax2.set_ylim(-0.233e-7, x[balance_x]+ 1.2e-7)

# do plots
ax1.plot(data['t'], real_volt_response, zorder=1)
ax2.plot(data['t'], accelerometer_volt_response, color='orange', zorder=2 , label='y: 0.0235')

plt.legend()

#plt.scatter(data['t'][balance_x], x[balance_x], color='gray', zorder=3)

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Real volt response (V)', color='C0')
ax2.set_ylabel('Accleramoter volt response (V)', color='orange')

plt.savefig('Volt_response.png', dpi=600)


