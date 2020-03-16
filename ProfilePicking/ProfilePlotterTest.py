import numpy as np
import matplotlib.pyplot as plt 
import scipy.ndimage



def onclick(event):
    global ix, iy
    ix,iy = event.xdata, event.ydata
    
    print ('x = {}, y = {}'.format(ix,iy))
    global coords

    coords.append((ix,iy))

    if len(coords) == 2:
        fig.canvas.mpl_disconnect(cid)
        plt.close(fig)
    return coords

N = 100
X, Y = np.mgrid[-3:3:complex(0, N), -2:2:complex(0, N)]
Z1 = np.exp(-X**2 - Y**2)
Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
Z = (Z1 - Z2) * 2

coords = []


fig, ax = plt.subplots(nrows=1)

ax.imshow(Z, origin="lower",cmap="RdBu_r")
cid = fig.canvas.mpl_connect("button_press_event", onclick)
print(coords)
plt.show()

num = 1000
x, y = np.linspace(coords[0][0], coords[1][0], num), np.linspace(coords[0][1], coords[1][1], num)
zi = scipy.ndimage.map_coordinates(Z, np.vstack((x,y)))

fig, ax = plt.subplots(nrows=2)

ax[0].imshow(Z,origin="lower",cmap="RdBu_r")
ax[0].plot([coords[0][0], coords[1][0]], [coords[0][1], coords[1][1]], 'ro-')
ax[0].axis('image')

ax[1].plot(zi)
plt.show()

