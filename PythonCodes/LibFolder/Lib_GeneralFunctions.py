import os
from glob import glob
import pickle, datetime
import matplotlib.pyplot as plt 


def CreateFolder(PathFolder):
    try:
        os.makedirs(PathFolder)
    except:
        pass

def SavePickleFile(OutputFolder, Filename, Variable):
    CreateFolder(OutputFolder)
    with open(OutputFolder + Filename,"wb") as f:
        pickle.dump(Variable,f)

def LoadPickleFile(FolderPath, Filename):
    with open(FolderPath + Filename,"rb") as f:
        LoadedVariable = pickle.load(f)
    return LoadedVariable 

def GetTodayDate():
    mylist = []
    today = datetime.date.today()
    mylist.append(today)
    return str(mylist[0]).replace("-","")   

def FontSizeControlFreak(SMALL_SIZE,MEDIUM_SIZE,BIGGER_SIZE):
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title