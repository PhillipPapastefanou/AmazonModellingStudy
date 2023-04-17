from SmartOutput.outputfile import OutputFile
from SmartOutput.outputfile import TimeDomain
import numpy as np
from enum import Enum


class Scenario(Enum):
    Control = 1
    Drought = 2


class TFE_Extractor:
    def __init__(self, path, nins):
        self.path = path
        self.nins = nins

        self.controlFiles = [[0 for x in range(nins)] for y in range(3)]
        self.droughtedFiles = [[0 for x in range(nins)] for y in range(3)]


        print(self.controlFiles[1][1])
        for i in range(0,nins):
            s = path + "\\" + str(i) + "run\\AnnuallyOut.nc"
            self.controlFiles[0][i] = s
            s = path + "\\" + str(i) + "run\\MonthlyOut.nc"
            self.controlFiles[1][i] = s
            s = path + "\\" + str(i) + "run\\DailyOut.nc"
            self.controlFiles[2][i] = s


        for i in range(nins,2*nins):
            s = path + "\\" + str(i) + "run\\AnnuallyOut.nc"
            self.droughtedFiles[0][i-nins] = s
            s = path + "\\" + str(i) + "run\\MonthlyOut.nc"
            self.droughtedFiles[1][i-nins] = s
            s = path + "\\" + str(i) + "run\\DailyOut.nc"
            self.droughtedFiles[2][i-nins] = s




    def PftNames(self, timeDomain):
        timeID = timeDomain.value - 1
        files = self.controlFiles[timeID]

        file = OutputFile(files[0])
        names = file.PFTVarNames
        file.Close()
        return names;

    def PatchNames(self, timeDomain):
        timeID = timeDomain.value - 1
        files = self.controlFiles[timeID]

        file = OutputFile(files[0])
        names = file.PatchVarNames
        file.Close()
        return names;


    def GetData(self, scenario, timeDomain, gridCellID, pftID, varName, yearBegin, yearEnd):


        timeID = timeDomain.value - 1
        if (scenario == Scenario.Control):
            files = self.controlFiles[timeID]
        elif(scenario == Scenario.Drought):
            files = self.droughtedFiles[timeID]
        else:
            print("Error, Invalid scenario")

        file = OutputFile(files[0])
        self.timeLength = len(file.GetTime(yearBegin,yearEnd))
        file.Close()

        data = np.zeros((self.nins, self.timeLength))

        for i in range(0, self.nins):
            file = OutputFile(files[i])
            data[i,:] = file.GetSinglePFT(varName, gridCellID, yearBegin, yearEnd, pftID)
            file.Close()
        return data

    def GetPatchData(self, scenario, timeDomain, gridCellID, varName, yearBegin, yearEnd):

        timeID = timeDomain.value - 1
        if (scenario == Scenario.Control):
            files = self.controlFiles[timeID]
        elif(scenario == Scenario.Drought):
            files = self.droughtedFiles[timeID]
        else:
            print("Error, Invalid scenario")

        file = OutputFile(files[0])
        self.timeLength = len(file.GetTime(yearBegin,yearEnd))
        file.Close()

        data = np.zeros((self.nins, self.timeLength))

        for i in range(0, self.nins):
            file = OutputFile(files[i])
            data[i,:] = file.GetPatchSingle(varName, gridCellID, yearBegin, yearEnd)
            file.Close()
        return data




    def GetTimes(self, scenario, timeDomain, yearBegin, yearEnd):

        timeID = timeDomain.value - 1
        if (scenario == Scenario.Control):
            filename = self.controlFiles[timeID][0]
        elif(scenario == Scenario.Drought):
            filename = self.droughtedFiles[timeID][0]
        else:
            print("Error, Invalid scenario")
            exit(1)

        file = OutputFile(filename)
        times = file.GetTime(yearBegin,yearEnd)
        file.Close()
        return times