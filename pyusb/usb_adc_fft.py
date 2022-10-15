from ciaa import Ciaa, CiaaStatusCode
import numpy as np
import matplotlib.pyplot as plt

# ADC clock
f_adc = 4.5e6 
# Conversion time for a 10 bit sample
conversion_time = 4 / f_adc
# Real sample length
N = 32
# Get sampling frequency
fs = 1 / (N * conversion_time)

# Initialize CIAA usb handler
ciaa = Ciaa()
# Look up for a conection
if ciaa.find() != CiaaStatusCode.CIAA_FOUND:
    # Raise error
    raise ValueError("CIAA not found")


# If it was found, print message
print("CIAA found!")
print("-" * 80)
print()

while True:

    # Get the FFT values
    fft = ciaa.read_json(100, debug = True)["values"]
    # Get list length
    n = len(fft)
    # Go from 0 to sampling frequency with N steps
    fn = np.linspace(0, fs / 2, n)
    # Take error into account
    #for i, f in enumerate(fn):
    #    fn[i] = 1.46 * f + 8.91e-1

    # Normalize amplitudes
    for i, s in enumerate(fft):
        fft[i] = s / n

    # Get rid of DC level
    fft[0] = 0.0

    # Get peak value
    peak = max(fft)
    # Get frequency index of max value
    index = fft.index(peak)

    print()
    print(f"Max value of {peak:.2f} in {fn[index]:.2f} Hz")
    print()

    # Plot settings
    plt.grid(True)
    plt.title("Sample FFT")
    plt.ylabel("Amplitude")
    plt.xlabel("Frequency [Hz]")
    plt.ylim([0, 3.3])
    plt.xlim([0, 20e3])
    # Plot up to fs / 2 the FFT values
    plt.stem(fn, fft, linefmt="blue", markerfmt="none")
    # Pause plot
    plt.pause(0.001)
    # Clear plot
    plt.clf()
