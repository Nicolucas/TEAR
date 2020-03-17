import numpy as np
import matplotlib.pyplot as plt 
import scipy.ndimage
import math


class MGridData:
    def __init__(self, ZGrid, NumPoints):
        self.coords = []
        self.ZGrid = ZGrid
        self.NumPoints = NumPoints
        self.fig = None
        self.cid = None


    def onclick(self, event):
        ix,iy = event.xdata, event.ydata
        
        print ('x = {}, y = {}'.format(ix,iy))

        self.coords.append((ix,iy))
        if len(self.coords) == 2:
            self.fig.canvas.mpl_disconnect(self.cid)
            plt.close(self.fig)


    def PlotNExtractProfile(self):
        fig, ax = plt.subplots(nrows=1)
        self.fig = fig
        ax.imshow(self.ZGrid, origin="lower",cmap ="RdBu_r")
        cid = self.fig.canvas.mpl_connect("button_press_event", self.onclick)
        self.cid = cid
        
        plt.show()

        x, y = np.linspace(self.coords[0][0], self.coords[1][0], self.NumPoints), \
               np.linspace(self.coords[0][1], self.coords[1][1], self.NumPoints)

        zi = scipy.ndimage.map_coordinates(Z, np.vstack((x,y)))
        Dist = math.sqrt((self.coords[1][0] - self.coords[0][0])**2.0 + (self.coords[1][1] - self.coords[0][1])**2.0)
        ZDist = np.linspace(0, Dist, self.NumPoints)

        fig, ax = plt.subplots(nrows=2)
        self.fig = fig
        self.cid = None

        ax[0].imshow(self.ZGrid,origin="lower",cmap="RdBu_r")
        ax[0].plot([self.coords[0][0], self.coords[1][0]], [self.coords[0][1], self.coords[1][1]], 'ro-')

        ax[1].plot(ZDist,zi)
        plt.show()


N = 100
X, Y = np.mgrid[-3:3:complex(0, N), -2:2:complex(0, N)]
Z1 = np.exp(-X**2 - Y**2)
Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
Z = (Z1 - Z2) * 2


MGrid = MGridData(NumPoints=1000,ZGrid=Z)
MGrid.PlotNExtractProfile()
print(MGrid.coords)
