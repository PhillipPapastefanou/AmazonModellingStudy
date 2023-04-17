import numpy as np
from netCDF4 import Dataset
from enum import  Enum
from Libs.Scaling import ScalerListToArray
from datetime import datetime
from datetime import timedelta
import Libs.SmartOutput.julian

class TimeDomain(Enum):
    Daily = 3
    Monthly = 2
    Annually = 1

class Dimension:
    def __init__(self, name, id, size):
        self.id = id
        self.name = name
        self.size = size

class OutputFile:
    def __init__(self, path):
        self.path = path
        self.nc = Dataset(path, 'r')


        try:
            ncTimeDim = self.nc.dimensions["time"]
        except:
            ncTimeDim = self.nc.dimensions["Time"]

        self.timeDim = Dimension(ncTimeDim.name, ncTimeDim._dimid, ncTimeDim.size)



        try:
            ncGridcellDim = self.nc.dimensions["Gridcell"]
        except:
            ncGridcellDim = self.nc.dimensions["gridcell"]


        self.gridcellDim = Dimension(ncGridcellDim.name, ncGridcellDim._dimid, ncGridcellDim.size)


        try:
            ncPftDim = self.nc.dimensions["Pfts_and_total"]
        except:
            ncPftDim = self.nc.dimensions["pfts_and_total"]

        self.pftDim = Dimension(ncPftDim.name, ncPftDim._dimid, ncPftDim.size)

        self.nPfts = self.pftDim.size


        try:
            try:
                ncStandDim = self.nc.dimensions["Stands"]
            except:
                ncStandDim = self.nc.dimensions["stands"]

            self.standDim = Dimension(ncStandDim.name, ncStandDim._dimid, ncStandDim.size)
            self.hasStandDim = True
        except:
            self.hasStandDim = False


        try:

            try:
                ncInsfileDim = self.nc.dimensions["Insfiles"]
            except:
                ncInsfileDim = self.nc.dimensions["insfiles"]
            self.insFileDim = Dimension(ncInsfileDim.name, ncInsfileDim._dimid, ncInsfileDim.size)
            self.hasInsfileDim = True
        except Exception:
            self.hasInsfileDim = False



        try:
            base_group = self.nc.groups['Base']
            timeVar = base_group.variables['Time']
            self.lons = base_group.variables['Longitude'][:]
            self.lats = base_group.variables['Latitude'][:]
            self.pfts = base_group.variables['Pfts'][:]
        except Exception:

            try:
                timeVar = self.nc.variables['Time']
                self.lons = self.nc.variables['Longitude'][:]
                self.lats = self.nc.variables['Latitude'][:]
                self.pfts = self.nc.variables['Pfts'][:]

            except Exception:
                timeVar = self.nc.variables['Base_Time']
                self.lons = self.nc.variables['Base_Longitude'][:]
                self.lats = self.nc.variables['Base_Latitude'][:]
                self.pfts = self.nc.variables['Base_Pfts'][:]

        sinceString  = timeVar.units

        timeStringData = str.split(sinceString, ' since ')
        timeToken = timeStringData[0]
        beginDate = timeStringData[1]

        beginDateData = str.split(beginDate, '-')
        self.yearBegin = int(beginDateData[2])
        self.monthBegin = int(beginDateData[1])
        self.dayBegin = int(beginDateData[0])

        if (timeToken == "days"):
            self.multiplier = 365
            self.timeDomain = TimeDomain.Daily
        elif (timeToken == "months"):
            self.multiplier = 12
            self.timeDomain = TimeDomain.Monthly
        elif(timeToken == "years"):
            self.multiplier = 1
            self.timeDomain = TimeDomain.Annually

        else:
            print('Invalid TimeDomain')

        timeData = timeVar[:]
        self.timeSteps = timeData.size
        self.yearEnd = self.yearBegin + int(timeData[self.timeSteps-1] / self.multiplier)

        patch_group = self.nc.groups['Patch-Out']
        pft_group = self.nc.groups['Pft-Out']

        self.PFTVarNames = []
        self.PatchVarNames= []

        for pftVar in pft_group.variables:
            self.PFTVarNames.append(pftVar)

        for patchVar in patch_group.variables:
            self.PatchVarNames.append(patchVar)
        self.nc.close()


    def GetSinglePFT(self, varName, gridCellid, yearBeginOfInterest, yearEndOfInterest, pftID):
        self.nc = Dataset(self.path, 'r')
        data = self._GetPFT(varName, gridCellid, gridCellid, yearBeginOfInterest, yearEndOfInterest, pftID, pftID)
        return data
        self.nc.close()


    def GetUnit(self, groupName, varName):
        self.nc = Dataset(self.path, 'r')
        group = self.nc.groups[groupName]
        var = group[varName]
        return var.getncattr('unit')
        self.nc.close()

    def GetMultiPFT(self, varName, gridcellBegin, gridcellEnd, yearBeginOfInterest, yearEndOfInterest, pftID_begin, pftID_end):
        data = self._GetPFT(varName = varName, gridcellBegin = gridcellBegin, gridcellEnd = gridcellEnd,
                            yearBeginOfInterest= yearBeginOfInterest,
                            yearEndOfInterest = yearEndOfInterest,
                            pftID_begin= pftID_begin, pftID_end = pftID_end)
        return data

    def GetSinglePFTTotal(self, varName, gridCellid, yearBeginOfInterest, yearEndOfInterest):
        totalDimIDPFT = self.pftDim.size - 1
        return self._GetPFT(varName, gridCellid, gridCellid, yearBeginOfInterest, yearEndOfInterest, totalDimIDPFT, totalDimIDPFT)

    def GetAllPFTTotal(self, varName, yearBeginOfInterest, yearEndOfInterest):
        #self.nc = Dataset(self.path, 'r')
        totalDimIDPFT = self.pftDim.size - 1
        data = self._GetPFT(varName, 0, self.gridcellDim.size, yearBeginOfInterest, yearEndOfInterest, totalDimIDPFT, totalDimIDPFT)
        return data

    def GetPft_allgridcells(self, varName, yearBeginOfInterest, yearEndOfInterest):
        #self.nc = Dataset(self.path, 'r')
        totalDimIDPFT = self.pftDim.size - 1
        data = self._GetPFT(varName, 0, self.gridcellDim.size, yearBeginOfInterest, yearEndOfInterest, 0, totalDimIDPFT)
        return data


    def GetTime(self, yearBegin, yearEnd):
        return self.split(yearBegin, yearEnd + 1, (yearEnd-yearBegin + 1)*self.multiplier)

    def _GetPFT(self, varName, gridcellBegin, gridcellEnd, yearBeginOfInterest, yearEndOfInterest, pftID_begin, pftID_end):

        self.nc = Dataset(self.path, 'r')
        timeOffset= (yearBeginOfInterest - self.yearBegin)*self.multiplier
        timeEnd = (yearEndOfInterest - self.yearBegin + 1)*self.multiplier

        ncVariable = self.nc.groups["Pft-Out"].variables[varName]

        offset = np.zeros(ncVariable.ndim,int)
        end = np.zeros(ncVariable.ndim, int)

        id = 0
        for varDim in ncVariable.dimensions:
            if (varDim == self.gridcellDim.name):
                offset[id] = gridcellBegin
                end[id] = gridcellEnd + 1
            if (varDim == self.pftDim.name):
                offset[id] = pftID_begin
                end[id] = pftID_end + 1
            if (varDim== self.standDim.name):
                offset[id] = 0
                end[id] = 1
            if (varDim == self.timeDim.name):
                offset[id] = timeOffset
                end[id] = timeEnd
            id+=1

        #Todo can this be done in a nicer way with generic offsets?
        data = ncVariable[offset[0]:end[0], offset[1]: end[1], offset[2]:end[2],offset[3]:end[3]]
        self.nc.close()
        return np.squeeze(data)

    def _GetPFT_Ins(self, varName, insfileBegin, insfileEnd, gridcellBegin, gridcellEnd, yearBeginOfInterest, yearEndOfInterest, pftID_begin, pftID_end):

        timeOffset= (yearBeginOfInterest - self.yearBegin)*self.multiplier
        timeEnd = (yearEndOfInterest - self.yearBegin + 1)*self.multiplier


        self.nc = Dataset(self.path, 'r')
        ncVariable = self.nc.groups["Pft-Out"].variables[varName]

        offset = np.zeros(ncVariable.ndim,int)
        end = np.zeros(ncVariable.ndim, int)

        id = 0
        for varDim in ncVariable.dimensions:
            if (varDim == self.gridcellDim.name):
                offset[id] = gridcellBegin
                end[id] = gridcellEnd + 1
            if (varDim == self.pftDim.name):
                offset[id] = pftID_begin
                end[id] = pftID_end + 1
            if (varDim== self.standDim.name):
                offset[id] = 0
                end[id] = 1
            if (varDim == self.timeDim.name):
                offset[id] = timeOffset
                end[id] = timeEnd
            if (varDim == self.insFileDim.name):
                offset[id] = insfileBegin
                end[id] = insfileEnd + 1
            id+=1
        #Todo can this be done in a nicer way with generic offsets?
        data = ncVariable[offset[0]:end[0], offset[1]: end[1], offset[2]:end[2] ,offset[3]:end[3], offset[4]:end[4]]

        self.nc.close()
        return data

    def _GetPatch(self, varName,insfileBegin, insfileEnd, gridcellBegin, gridcellEnd, yearBeginOfInterest, yearEndOfInterest):


        self.nc = Dataset(self.path, 'r')
        timeOffset= (yearBeginOfInterest - self.yearBegin)*self.multiplier
        timeEnd = (yearEndOfInterest - self.yearBegin + 1)*self.multiplier

        ncVariable = self.nc.groups["Patch-Out"].variables[varName]

        offset = np.zeros(ncVariable.ndim,int)
        end = np.zeros(ncVariable.ndim, int)

        id = 0
        for varDim in ncVariable.dimensions:
            if (varDim == self.gridcellDim.name):
                offset[id] = gridcellBegin
                end[id] = gridcellEnd + 1
            if (varDim== self.standDim.name):
                offset[id] = 0
                end[id] = 1
            if (varDim == self.timeDim.name):
                offset[id] = timeOffset
                end[id] = timeEnd
            if (varDim == self.insFileDim.name):
                offset[id] = insfileBegin
                end[id] = insfileEnd + 1
            id+=1

        data = ncVariable[offset[0]:end[0], offset[1]: end[1], offset[2]:end[2] ,offset[3]:end[3]]
        self.nc.close()
        return data



    def GetPatchSingle(self, varName, gridcell, yearBeginOfInterest, yearEndOfInterest):
        self.nc = Dataset(self.path, 'r')
        data = self.GetPatchMulti(varName, gridcell, gridcell, yearBeginOfInterest, yearEndOfInterest)
        self.nc.close()
        return data

    def GetPatchSingleAllGridcells(self, varName, yearBeginOfInterest, yearEndOfInterest):
        self.nc = Dataset(self.path, 'r')
        data = self.GetPatchMulti(varName, 0, self.gridcellDim.size, yearBeginOfInterest, yearEndOfInterest,0)
        #self.nc.close()
        return data


    def GetPatchMulti(self, varName, gridcellBegin, gridcellEnd, yearBeginOfInterest, yearEndOfInterest, insfileID):
        self.nc = Dataset(self.path, 'r')
        timeOffset= (yearBeginOfInterest - self.yearBegin)*self.multiplier
        timeEnd = (yearEndOfInterest - self.yearBegin + 1)*self.multiplier

        if not varName in self.PatchVarNames:
            print("Variable not found in output")
            exit(1)


        ncVariable = self.nc.groups["Patch-Out"].variables[varName]

        offset = np.zeros(ncVariable.ndim,int)
        end = np.zeros(ncVariable.ndim, int)

        id = 0
        for varDim in ncVariable.dimensions:
            if (varDim == self.gridcellDim.name):
                offset[id] = gridcellBegin
                end[id] = gridcellEnd + 1
            if (varDim== self.standDim.name):
                offset[id] = 0
                end[id] = 1
            if (varDim == self.timeDim.name):
                offset[id] = timeOffset
                end[id] = timeEnd
            #if (varDim == self.insFileDim.name):
           #     offset[id] = insfileID
            #    end[id] = insfileID+1
            id+=1

        data = ncVariable[offset[0]:end[0], offset[1]: end[1], offset[2]:end[2]]
        self.nc.close()
        return np.squeeze(data)






    def Close(self):
        self.nc.close()


    def split(self, x, y, n):
        vals = np.zeros(n)
        for i in range(0,n):
            vals[i] = x + i*(y-x)/n
        return vals

    def CreateImage(self, res, data, lat_flip):

        lonScaler = ScalerListToArray(self.lons, res, False)
        latScaler = ScalerListToArray(self.lats, res, lat_flip)
        lonIndexes = lonScaler.Indexes
        latIndexes = latScaler.Indexes
        xlen = lonScaler.len
        ylen = latScaler.len
        self.IMG_extent = [lonScaler.min, lonScaler.max, latScaler.min, latScaler.max]

        image = np.empty((ylen + 1, xlen + 1,)) * np.nan
        for i in range(len(data)):
            image[latIndexes[i], lonIndexes[i]] = data[i]

        return image













