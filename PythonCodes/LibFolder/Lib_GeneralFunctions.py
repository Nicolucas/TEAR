import os
from glob import glob
import pickle, datetime


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
