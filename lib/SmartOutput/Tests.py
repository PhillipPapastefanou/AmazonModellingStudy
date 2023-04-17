from Libs.SmartOutput.outputfile import OutputFile
import numpy as np
import matplotlib.pyplot as plt

testFile = OutputFile(r'F:\SourceTreeRepos\guess4_hydraulics2021\build\Release\9run\AnnuallyOut.nc')

print(testFile.PFTVarNames)

data = testFile.GetMultiPFT("cmasstotal", 0,1946, 1990, 1990, 0, 0)
print(np.mean(data))


testFile = OutputFile('F:\\Simulations\\V1.5.2.0sd_AGU_PftList_npatch=2\\0run\\MonthlyOut.nc')
data = testFile.GetMultiPFT("nppp", 1990, 1990, 0)
print(np.sum(data))
testFile.Close()

testFile = OutputFile('F:\\Simulations\\V1.5.2.0sd_AGU_PftList_npatch=2\\0run\\AnnuallyOut.nc')
data = testFile.GetSinglePFT("gpp", 1, 1990, 2010, 0)
print(np.sum(data))
testFile.Close()

testFile = OutputFile('F:\\Simulations\\V1.5.2.0sd_AGU_PftList_npatch=2\\0run\\AnnuallyOut.nc')
data = testFile.GetMultiPFT("gpp", 0, 5, 1990, 2010, 3)
print(data[1])

print(np.mean(data, axis=1))
testFile.Close()


testFile = OutputFile('F:\\Simulations\\V1.5.2.0sd_AGU_PftList_npatch=2\\0run\\AnnuallyOut.nc')
data = testFile.GetSinglePFTTotal("gpp", 5, 2000, 2005)
timeValues = testFile.GetTime(2000, 2005)



fig, ax = plt.subplots()
ax.plot(timeValues, data)
plt.show()

testFile.Close()
