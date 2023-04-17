from Libs.OutputSampling import StratifiedSample
from Libs.SmartOutput.outputfile import TimeDomain
from Libs.Standard import MCWDCalc
import pandas
import numpy as np

class AmazonBasinStratified:
    def __init__(self, filePath, configPath, configheader = None):
        self.stratifiedOutput = StratifiedSample(filePath)


        self.inputDataFrame = pandas.read_csv(configPath,
                     sep = '\t',
                     header=configheader)

        self.inputData = self.inputDataFrame.values

        self.filePath = filePath

        self.npft = self.stratifiedOutput.nPfts
        self.nIns = self.stratifiedOutput.nInsfiles


    def CreateImage(self, data):
        return  self.stratifiedOutput.CreateImage(0.5, data)

    def GetDataSlices(self, varName, yearBefore, yearAfter, pftId = -1):
        if pftId == -1:
            pftId = self.npft - 1

        dataBefore = self.stratifiedOutput.GetAllGridcells(varName, yearBefore, yearBefore, 0, self.nIns, pftId)
        dataAfter = self.stratifiedOutput.GetAllGridcells(varName, yearAfter, yearAfter, 0, self.nIns, pftId)

        return np.squeeze(dataBefore),  np.squeeze(dataAfter)

    def GetDataSlicesRange(self, varName, yearBefore, yearAfter, pftId = -1):
        if pftId == -1:
            pftId = self.npft - 1

        data = self.stratifiedOutput.GetAllGridcells(varName, yearBefore, yearAfter, 0, self.nIns, pftId)
        return np.squeeze(data)

    def GetDataSlicesYear(self, varName, year, pftId = -1):
        if pftId == -1:
            pftId = self.npft - 1

        data = self.stratifiedOutput.GetAllGridcells(varName, year, year, 0, self.nIns, pftId)

        return np.squeeze(data)

    def GetDataSlicesPatchYear(self, varName, year, pftId = -1):
        if pftId == -1:
            pftId = self.npft - 1

        data =  self.stratifiedOutput.GetAllGridcellsPatch(varName, year, year, 0, self.nIns)


        return np.squeeze(data)

    def GetDataSlicesYearPft(self, varName, year, pftID_begin, pftID_end):

        data = self.stratifiedOutput.GetAllGridcellsPft(varName, year, year, 0, self.nIns, pftID_begin, pftID_end)

        return np.squeeze(data)

    def GetDataSlicesCond(self, varName, yearBefore, yearAfter, less, pftId = -1):

        if pftId == -1:
            pftId = self.npft - 1

        dataBefore = self.stratifiedOutput.GetAllGridcells(varName, yearBefore, yearBefore, 0, self.nIns, pftId)
        dataAfter = self.stratifiedOutput.GetAllGridcells(varName, yearAfter, yearAfter, 0, self.nIns, pftId)

        ind = dataBefore < less

        dataBefore[ind] = np.NAN
        dataAfter[ind] = np.NAN

        return np.squeeze(dataBefore),  np.squeeze(dataAfter)


    def GetDataSlicesCondContinous(self, varName, yearBefore, yearAfter, less, pftId = -1):

        if pftId == -1:
            pftId = self.npft - 1

        data = self.stratifiedOutput.GetAllGridcells(varName, yearBefore, yearAfter, 0, self.nIns, pftId)

        ind = data < less

        data[ind] = np.NAN

        return np.squeeze(data)


    def GetDataSlicesPatch(self, varName, yearBefore, yearAfter):

        dataBefore = self.stratifiedOutput.GetAllGridcellsPatch(varName, yearBefore, yearBefore, 0, self.nIns)
        dataAfter = self.stratifiedOutput.GetAllGridcellsPatch(varName, yearAfter, yearAfter, 0, self.nIns)

        return np.squeeze(dataBefore),  np.squeeze(dataAfter)

    def GetDataPatchRange(self, varName, yearBefore, yearAfter):

        data = self.stratifiedOutput.GetAllGridcellsPatch(varName, yearBefore, yearAfter, 0, self.nIns)

        return np.squeeze(data)

    def GetDataPatch(self, varName, yearBegin, yearEnd, insfile):

        dataBefore = self.stratifiedOutput.GetAllGridcellsPatch(varName, yearBegin, yearEnd, insfile, insfile)

        return np.squeeze(dataBefore)


    def GetDataPatchAllIns(self, varName, yearBegin, yearEnd):

        dataBefore = self.stratifiedOutput.GetAllGridcellsPatch(varName, yearBegin, yearEnd, 0, self.nIns)

        return np.squeeze(dataBefore)



    def GetMCWD(self, year):

        return MCWDCalc(self.filePath, year).DeltaMCWD