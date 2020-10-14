import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

path = "/Users/a27inchMac/Desktop/SARS2 qPCR Data/Analytical_Precision_20200909/"

files = []
files.append(
    "SARS2_PLATE_ONE_27893-27938_FAM_HEX_2020-08-28 11-31-29_CT042716 -  Quantification Amplification Results.xlsx")
files.append(
    "SARS2_PLATE_TWO_27939-28000_FAM_HEX_2020-08-28 11-35-45_BR003296 -  Quantification Amplification Results.xlsx")
files.append(
    "SARS2_PLATE_THREE_27957-27973_28001-28028_FAM_HEX_2020-08-28 14-22-01_CT042716 -  Quantification Amplification Results.xlsx")
files.append(
    "SARS2_PLATE_FOUR_28029-28074_FAM_HEX_2020-08-28 16-46-55_BR003296 -  Quantification Amplification Results.xlsx")
files.append(
    "SARS2_PLATE_FIVE_28075-28138_FAM_HEX_2020-08-28 17-30-22_CT042716 -  Quantification Amplification Results.xlsx")
# print(files)

colNames = ['Plate 1', 'Plate 2', 'Plate 3', 'Plate 4', 'Plate 5']
rownames = range(1, 46)

wellEplusN1 = "G12"
wellEplusN3 = "H12"

tmpN1 = []
tmpN3 = []

for file in files:
    dataFam = pd.read_excel(path + file, sheet_name='FAM')
    tmpN1.append(pd.DataFrame(dataFam, columns=[wellEplusN1]).values)
    tmpN3.append(pd.DataFrame(dataFam, columns=[wellEplusN3]).values)

tmpN1 = np.array(tmpN1).T[0]
cn = ['N1 (' + colname + ')' for colname in colNames]
dataN1 = pd.DataFrame(tmpN1, columns=cn, index=rownames)

tmpN3 = np.array(tmpN3).T[0]
cn = ['N3 (' + colname + ')' for colname in colNames]
dataN3 = pd.DataFrame(tmpN3, columns=cn, index=rownames)

BT = 180  # manually set baseline threshold value
dataBT = pd.DataFrame([BT, BT], columns=['BT'], index=[1, 45])

# print(dataBT)

CqLabel = 'Cq @ BT = ' + str(BT)

tmp = []
for (colName, colData) in dataN1.iteritems():
    diffBT = colData.values - BT
    idxBeforeCrossing = np.where(diffBT < 0)[0][-1]
    Cq = dataN1.index[idxBeforeCrossing] + (BT - colData.values[idxBeforeCrossing]) / (
            colData.values[idxBeforeCrossing + 1] - colData.values[idxBeforeCrossing])
    tmp.append(Cq)

for (colName, colData) in dataN3.iteritems():
    diffBT = colData.values - BT
    idxBeforeCrossing = np.where(diffBT < 0)[0][-1]
    Cq = dataN3.index[idxBeforeCrossing] + (BT - colData.values[idxBeforeCrossing]) / (
            colData.values[idxBeforeCrossing + 1] - colData.values[idxBeforeCrossing])
    tmp.append(Cq)

tmp = np.array(tmp)
tmp2 = np.ones_like(tmp) * BT
dataCq = pd.DataFrame(tmp2, columns=[CqLabel], index=tmp)

# print(dataCq)

plt.rcParams['figure.figsize'] = [10, 9]
ax = plt.gca()

dataN1.plot(kind='line', marker='.', ax=ax)
dataN3.plot(kind='line', marker='+', ax=ax)
dataBT.plot(kind='line', ax=ax)
dataCq.plot(kind='line', marker='*', markerfacecolor='k', markeredgecolor='k', linestyle='none', ax=ax)
plt.xlabel('Cycle')
plt.ylabel('RFU')
plt.title('Positive Extraction Controls for Five Plates on 8/28/20')
plt.ylim([0, 400])
plt.xlim([30, 35])

plt.show()
