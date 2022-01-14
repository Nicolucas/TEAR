from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,mark_inset)
from obspy.imaging import spectrogram
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Lib_SigmoidProcessing import *

def SortReceiver(Receiver,AmplitudeList):
    timePre = np.asarray(Receiver.Time)
    AmplitudePre = np.asarray(AmplitudeList)

    TimeIdx = timePre.argsort()
    time = timePre[TimeIdx]
    Amplitude = AmplitudePre[TimeIdx]
    return time, Amplitude


def Receiver2Spectrogram(Receiver, AmplitudeList, yLabel="", **kwargs):
    timePre = np.asarray(Receiver.Time)
    AmplitudePre = np.asarray(AmplitudeList)

    TimeIdx = timePre.argsort()
    time = timePre[TimeIdx]
    Amplitude = AmplitudePre[TimeIdx]

    samplingrate = len(Amplitude)/time.max()

    fig = plt.figure(figsize = (10, 5),dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(1, 1)
    ax = fig.add_subplot(gs[:, :])
    ax2 = ax.twinx()

    spectrogram.spectrogram(Amplitude,samplingrate,log=True,axes=ax,show=False,**kwargs)
    ax.set_ylabel("Frequency [Hz]")
    ax.set_xlabel("time [s]")

    mappable = ax.collections[0]
    cbaxes = inset_axes(ax,width="40%",height="4%",loc=1, borderpad=2)
    plt.colorbar(mappable=mappable, cax=cbaxes, orientation="horizontal", label="Amplitude spectrum")
    cbaxes.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    cbaxes.xaxis.set_label_position('top')


    cbaxes.xaxis.label.set_color('white')
    cbaxes.tick_params(labelcolor='white')

    ax2.plot(time, Amplitude, lw=2, c='r')
    #ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax2.set_ylim([0,9])
    ax2.set_ylabel(yLabel)

    fig.suptitle('Spectrogram, st {}'.format(Receiver.Coord))

    return fig, ax, ax2



def Receiver2GradientSpectrogram(Receiver, AmplitudeList, yLabel="", **kwargs):
    timePre = np.asarray(Receiver.Time)
    AmplitudePre = np.asarray(AmplitudeList)

    TimeIdx = timePre.argsort()
    time = timePre[TimeIdx]
    Amplitude = AmplitudePre[TimeIdx]

    Amplitude = np.asarray(pd.Series(np.gradient(Amplitude, time).tolist()))

    samplingrate = len(Amplitude)/time.max()
    
    AxOut = kwargs.get("ax",None)
    if  AxOut == None:

        fig = plt.figure(figsize = (10, 5),dpi=300, constrained_layout=True)
        gs = fig.add_gridspec(1, 1)
        ax = fig.add_subplot(gs[:, :])
        
        fig.suptitle('Spectrogram, st {}'.format(Receiver.Coord))
    else:
        kwargs.pop("ax")
        ax = AxOut
        
        
    ax2 = ax.twinx()

    spectrogram.spectrogram(Amplitude,samplingrate,log=True,axes=ax,show=False,**kwargs)
    ax.set_ylabel("Frequency [Hz]")
    ax.set_xlabel("time [s]")

    mappable = ax.collections[0]
    cbaxes = inset_axes(ax,width="40%",height="4%",loc=1, borderpad=2)
    plt.colorbar(mappable=mappable, cax=cbaxes, orientation="horizontal", label="Amplitude spectrum")
    cbaxes.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    cbaxes.xaxis.set_label_position('top')


    cbaxes.xaxis.label.set_color('white')
    cbaxes.tick_params(labelcolor='white')

    ax2.plot(time, Amplitude, lw=1, c='k')
    #ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax2.set_ylim([0,9])
    ax2.set_ylabel(yLabel)

    if AxOut == None:
        return fig, ax, ax2
    else:
        return ax, ax2

