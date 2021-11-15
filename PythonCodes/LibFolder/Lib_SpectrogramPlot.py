from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,mark_inset)
from obspy.imaging import spectrogram
import numpy as np
import matplotlib.pyplot as plt
from Lib_SigmoidProcessing import *

def Receiver2Spectrogram(Receiver, Amplitude, yLabel=""):
    time = np.asarray(Receiver.Time)
    samplingrate = len(Amplitude)/Receiver.Time.max()

    fig = plt.figure(figsize = (10, 5),dpi=300, constrained_layout=True)
    gs = fig.add_gridspec(1, 1)
    ax = fig.add_subplot(gs[:, :])
    ax2 = ax.twinx()

    spectrogram.spectrogram(Amplitude,samplingrate,log=True,axes=ax,show=False)
    ax.set_xlim([0.15,3.8])
    ax.set_ylim([0.19,1e2])
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
    ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax2.set_ylim([0,10])
    ax2.set_ylabel(yLabel)

    ax.set_title('Spectrogram, st{}'.format(Receiver.Coord))