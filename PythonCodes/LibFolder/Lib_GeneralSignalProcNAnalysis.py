import numpy as np
from scipy import fft
from scipy import signal

# Function to calculate the root-mean-square of a target list of data against a target list of data
def rmse(predictions, targets):
    return np.sqrt(((np.asarray(predictions)-np.asarray(targets))**2).mean()).tolist()

# Function to filter using a Butterworth filter
def Butterworth(Signal, Type = "low",CutoffFrequency = 7., SamplingFrequency = 200):
    NormFrequency = CutoffFrequency / (SamplingFrequency / 2.0)
    
    b,a = signal.butter(1, NormFrequency, Type)
    output =  signal.filtfilt(b, a, np.asarray(Signal)).tolist()
    return output

# Function to calculate the FFT of a dataset
def CalculateFFT(Data, sampling = 0.005, TMax = 4.0):
    N=int(TMax/sampling)
    yf = fft(Data)
    xf = np.linspace(0.0, 1.0/(2.0*sampling), N//2)
    yf = 2.0/N * np.abs(yf[0:N//2])
    
    return(xf,yf)