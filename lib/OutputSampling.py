from Libs.SmartOutput.outputfile import OutputFile
from Libs.SmartOutput.outputfile import TimeDomain

import numpy as np
from enum import Enum
from Libs.Scaling import ScalerListToArray

class Scenario(Enum):
    Control = 1
    Drought = 2

class StratifiedSample:

    def __init__(self, path):
        self.path = path
        file = OutputFile(self.path)
        self.lons = file.lons
        self.lats = file.lats
        self.nInsfiles = file.insFileDim.size
        self.nPfts = file.pftDim.size

    def PftNames(self):
        file = OutputFile(self.path)
        names = file.PFTVarNames
        #file.Close()
        return names;

    def GetUnit(self, group, var):
        file = OutputFile(self.path)
        name = file.GetUnit(group, var)
        #file.Close()
        return name;

    def PatchNames(self):
        file = OutputFile(self.path)
        names = file.PatchVarNames
        #file.Close()
        return names;

    def GetAllIns(self, var, year, gridcell):
        file = OutputFile(self.path)
        data = file._GetPFT_Ins(var, 0, file.insFileDim.size, yearBeginOfInterest=year, yearEndOfInterest= year,
                                gridcellBegin=gridcell, gridcellEnd=gridcell,
                                pftID_begin= file.pftDim.size-1, pftID_end= file.pftDim.size-1)
        file.Close()
        return data


    def GetSelection(self, var, yearBegin, yearEnd, gridcell, insfileBegin, insfileEnd, pftID):
        file = OutputFile(self.path)
        data = file._GetPFT_Ins(var, insfileBegin, insfileEnd,
                                pftID_begin= pftID,
                                pftID_end=pftID,
                                yearBeginOfInterest=yearBegin, yearEndOfInterest= yearEnd, gridcellBegin=gridcell, gridcellEnd=gridcell)
        file.Close()
        return data

    def GetSelectionPatch(self, var, yearBegin, yearEnd, gridcell, insfileBegin, insfileEnd):
        file = OutputFile(self.path)
        data = file._GetPatch(var, insfileBegin, insfileEnd,
                                yearBeginOfInterest=yearBegin, yearEndOfInterest= yearEnd, gridcellBegin=gridcell, gridcellEnd=gridcell)
        file.Close()
        return data


    def GetAllGridcells(self, var, yearBegin, yearEnd, insfileBegin, insfileEnd, pftID):
        file = OutputFile(self.path)
        data = file._GetPFT_Ins(var, insfileBegin, insfileEnd, pftID_begin=pftID,
                            pftID_end=pftID, yearBeginOfInterest=yearBegin,
                            yearEndOfInterest= yearEnd,
                            gridcellBegin=0, gridcellEnd=file.gridcellDim.size)
        return data

    def GetAllGridcellsPft(self, var, yearBegin, yearEnd, insfileBegin, insfileEnd, pftID_begin, pftID_end):
        file = OutputFile(self.path)
        data = file._GetPFT_Ins(var, insfileBegin, insfileEnd, pftID_begin=pftID_begin,
                            pftID_end=pftID_end, yearBeginOfInterest=yearBegin,
                            yearEndOfInterest= yearEnd,
                            gridcellBegin=0, gridcellEnd=file.gridcellDim.size)
        #file.Close()
        return data

    def GetAllGridcellsPatch(self, var, yearBegin,yearEnd, insfileBegin, insfileEnd):
        file = OutputFile(self.path)
        data = file._GetPatch(var, insfileBegin, insfileEnd, yearBeginOfInterest=yearBegin, yearEndOfInterest= yearEnd,gridcellBegin=0, gridcellEnd=file.gridcellDim.size)
        #file.Close()
        return data


    def GetAllGridcellsPatchNoIns(self, var, yearBegin,yearEnd):
        file = OutputFile(self.path)
        data = file.GetPatchMulti(varName = var,yearBeginOfInterest=yearBegin, yearEndOfInterest= yearEnd,gridcellBegin=0, gridcellEnd=file.gridcellDim.size)

        #file.Close()
        return data


    def CreateImage(self, res, data):

        lonScaler = ScalerListToArray(self.lons, res, False)
        latScaler = ScalerListToArray(self.lats, res, True)
        lonIndexes = lonScaler.Indexes
        latIndexes = latScaler.Indexes
        xlen = lonScaler.len
        ylen = latScaler.len
        self.IMG_extent = [lonScaler.min, lonScaler.max, latScaler.min, latScaler.max ]

        image = np.empty((ylen + 1, xlen + 1,)) * np.nan
        for i in range(len(self.lons)):
            image[latIndexes[i], lonIndexes[i]] = data[i]

        return image



