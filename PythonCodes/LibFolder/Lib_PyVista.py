import pyvista as pv
import numpy as np

def pyvistaArraySorting(StressMesh, data, basis_degree):
    point= StressMesh.points
    Xpoint = point[:,0]
    Ypoint = point[:,1]
    
    nqp = basis_degree + 1
    nqp2d = nqp**2
    
    nx = len(Xpoint)
    ny = len(Ypoint)
    mx = int(np.sqrt(nx/nqp2d))
    my = int(np.sqrt(ny/nqp2d))
    
    # 2D Array initialization
    Xmatrix = np.zeros((mx*nqp,my*nqp))
    Ymatrix = np.zeros((mx*nqp,my*nqp))
    Zmatrix = np.zeros((mx*nqp,my*nqp))
    
    # separate per element
    XElements = np.hsplit(Xpoint, nx/nqp2d)
    YElements = np.hsplit(Ypoint, nx/nqp2d)
    ZElements = np.hsplit(data, nx/nqp2d)
    
    # Rearrange the point data so that it increases left to right and bottom up
    for j in range(1,my+1):
        Xmatrix[(j-1)*nqp:j*nqp,:mx*nqp] = np.stack([i.reshape(nqp,-1) for i in np.array(XElements[mx*(j-1):mx*j])],axis=1).reshape(-1,nqp*mx)
        Ymatrix[(j-1)*nqp:j*nqp,:mx*nqp] = np.stack([i.reshape(nqp,-1) for i in np.array(YElements[mx*(j-1):mx*j])],axis=1).reshape(-1,nqp*mx)
        Zmatrix[(j-1)*nqp:j*nqp,:mx*nqp] = np.stack([i.reshape(nqp,-1) for i in np.array(ZElements[mx*(j-1):mx*j])],axis=1).reshape(-1,nqp*mx)
        
    return Xmatrix,Ymatrix,Zmatrix