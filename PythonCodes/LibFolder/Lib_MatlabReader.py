from scipy.io import loadmat


def ReadSEM2DPACKFieldStructure(FolderPath,File,ColumnNames=["dt","nt","nsta","x","z","ux","uz"]):
    MatlabStruct = loadmat(FolderPath+File)
    MatlabStruct.pop('__header__')
    MatlabStruct.pop( '__version__')
    MatlabStruct.pop( '__globals__')
    
    MatlabStructData = MatlabStruct[list(MatlabStruct.keys())[0]]

    dt    = MatlabStructData[ColumnNames[0]][0][0][0][0]
    nt    = MatlabStructData[ColumnNames[1]][0][0][0][0]
    nsta  = MatlabStructData[ColumnNames[2]][0][0][0][0]

    x    = MatlabStructData[ColumnNames[3]][0][0]
    z    = MatlabStructData[ColumnNames[4]][0][0]
    ux   = MatlabStructData[ColumnNames[5]][0][0]
    uz   = MatlabStructData[ColumnNames[6]][0][0]


    ThisDictStruct = {"dt":dt,
                      "nt":nt,
                      "nsta":nsta,
                      "x": x,
                      "z": z,
                      "Field_x":ux,
                      "Field_z":uz}

    return ThisDictStruct

