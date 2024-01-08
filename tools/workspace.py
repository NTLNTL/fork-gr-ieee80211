import numpy as np
import matplotlib.pyplot as plt

# Parameters
w_0 = 2 * np.pi  # Angular frequency
theta = np.pi / 4  # Phase angle

# Frequency vector
w = np.linspace(-4 * np.pi, 4 * np.pi, 1000)

# Expression
signal = np.pi * (np.sinc((w - w_0) / np.pi) * np.exp(1j * theta) +
                 np.sinc((w + w_0) / np.pi) * np.exp(-1j * theta))

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(w, np.real(signal), label='Real Part')
plt.plot(w, np.imag(signal), label='Imaginary Part')
plt.title(r'$\pi[\delta(w-w_0)e^{j\theta} + \delta(w+w_0)e^{-j\theta}]$ in Frequency Domain (Approximation)')
plt.xlabel('Frequency (rad/s)')
plt.legend()
plt.grid(True)
plt.show()
