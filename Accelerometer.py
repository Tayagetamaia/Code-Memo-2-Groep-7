# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define physical parameters for the system
k = 84  # Spring constant in N/m   (42)
m = 1.9372e-6  # Mass of the object in kg
y = np.sqrt(4 * m * k) / 1.5  # Damping coefficient in kg/s     (0.045101)   np.sqrt(4 * m * k)

# Read in the position data from a text file, assuming whitespace delimiter remove or add the # to choose which position file to use
data = pd.read_csv('Posities\posities_1_Team_07.txt', sep='\s+')    
#data = pd.read_csv('Posities\posities_2_Team_07.txt', sep='\s+')   

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

dt_peak =  '{:.2e}'.format(data['t'][x.argmax()] - data['t'][data['a (m/s^2)'].idxmax()])


fig, ax1 = plt.subplots()

ax2 = ax1.twinx()


ax1.plot(data['t'], data['a (m/s^2)'])
ax2.plot(data['t'], x, color='orange')

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Real acceleration (m/s^2)', color='C0')
ax2.set_ylabel('Accleramoter acceleration (m/s^2)', color='orange')

plt.savefig('Acceleration_vs_time.png', dpi=600)



integral_real = 0
integral_accelerometer = 0
x_scale = data['a (m/s^2)'][data['a (m/s^2)'].idxmax()] / x[x.argmax()]

for i in range(len(data['t'])-1):
    integral_real += (data['t'][i+1] - data['t'][i]) * ((data['a (m/s^2)'][i] + data['a (m/s^2)'][i+1]) / 2)
    integral_accelerometer += (data['t'][i+1] - data['t'][i]) * ((x[i] + x[i+1]) / 2) * x_scale

accuracy = 100 - abs((integral_accelerometer - integral_real) / integral_real) * 100


print('delta time peaks:', dt_peak, 's')
print('accuracy:', round(accuracy, 2), "%")