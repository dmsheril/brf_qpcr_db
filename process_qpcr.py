import pandas as pd
import numpy as np
import glob
import time

print("This is the top of the process_qpcr module")


def read_data_from_curve_file(file):
    print("Reading data from curve file " + file)
    dataFam = pd.read_excel(file, sheet_name='FAM')
    dataHex = pd.read_excel(file, sheet_name='HEX')
    dataInfo = pd.read_excel(file, sheet_name='Run Information')
    return dataFam, dataHex, dataInfo


def read_data_from_cq_file(file):
    print("Reading data from Cq file " + file)
    dataCq = pd.read_excel(file, sheet_name='0')
    dataInfo = pd.read_excel(file, sheet_name='Run Information')
    return dataCq, dataInfo


def process_data(data):
    print("Beginning data processing...")
    modified_data = data + " that has been modified"
    time.sleep(2)
    print("Data processing finished.")
    return modified_data


def write_data_to_file(outfile, dataCq, dataCurve):
    print("Writing Cq & curve data to excel file: " + outfile)
    writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
    dataCq.to_excel(writer, sheet_name="Cq", index=False)
    dataCurve.to_excel(writer, sheet_name="Curve", index=False)
    writer.save()


def standardizeWellName(nameIn):
    nameAlpha = nameIn[0]
    nameNumber = int(nameIn[1:])
    nameOut = nameAlpha + "{:02d}".format(nameNumber)
    return nameOut


def formatCurveData(curveDataIn, filename, fluorName):
    curveData = curveDataIn.copy()
    curveData = curveData.transpose()
    curveData.columns = ["Cycle{:02d}".format(n) for n in range(1, curveData.shape[1] + 1)]
    curveData = curveData.drop(curveData.index[0:2])
    curveData.index.name = "Well"
    curveData.reset_index(inplace=True)
    curveData["Well"] = [standardizeWellName(name) for name in curveData['Well']]
    colFluor = [fluorName] * curveData.shape[0]
    curveData.insert(loc=1, column="Fluor", value=colFluor)
    colFilename = [filename] * curveData.shape[0]
    curveData.insert(loc=2, column="File name", value=colFilename)

    return curveData


def formatCqData(cqDataIn, filename, notes, runStart, baseSerial):
    cqData = cqDataIn.copy()
    colsToKeep = ["Well", "Fluor", "Target", "Content", "Sample", "Biological Set Name", \
                  "Cq", "Starting Quantity (SQ)", "Cq Std. Dev", "SQ Std. Dev"]
    cqData = cqData[colsToKeep]
    colFilename = [filename] * cqData.shape[0]
    colNotes = [notes] * cqData.shape[0]
    colRunStart = [runStart] * cqData.shape[0]
    colBaseSerial = [baseSerial] * cqData.shape[0]
    cqData.insert(cqData.shape[1], column="Date (Run start)", value=colRunStart)
    cqData.insert(cqData.shape[1], column="File name", value=colFilename)
    cqData.insert(cqData.shape[1], column="Base Serial #", value=colBaseSerial)
    cqData.insert(cqData.shape[1], column="Notes", value=colNotes)
    return cqData


def formatInfoData(infoDataIn):
    infoData = infoDataIn.copy()
    infoData = infoData.T.reset_index().T
    infoData.reset_index(inplace=True, drop=True)
    return infoData


def findBTforCq(curvedata, Cq):
    # find the baseline threshold given a curve and a Cq value
    x1 = np.floor(Cq)
    x2 = np.ceil(Cq)
    y1 = curvedata.loc[[int(x1)], ["RFU"]].values[0]
    y2 = curvedata.loc[[int(x2)], ["RFU"]].values[0]
    m = (y2 - y1) / (x2 - x1)
    b = y2 - m * x2
    BT = m * Cq + b
    BT = BT[0]
    return BT


def findCqforBT(curvedata, BT):
    # find Cq given a curve and a baseline threshold value
    diffBT = curvedata - BT
    if np.any(diffBT > 0):
        idxBeforeCrossing = np.where(diffBT < 0)[0][-1]
        Cq = curvedata.index[idxBeforeCrossing] + (BT - curvedata.values[idxBeforeCrossing]) \
             / (curvedata.values[idxBeforeCrossing + 1] - curvedata.values[idxBeforeCrossing])
        Cq = Cq[0]
    else:
        Cq = np.nan
    return Cq


def main():
    path = "/Users/a27inchMac/Desktop/SARS2 qPCR Data/data_to_process/"
    # path = "/Users/a27inchMac/Desktop/SARS2 qPCR Data/Pipeline_test_20201005/"

    allCurve = []
    allCq = []

    files = glob.glob(path + '[!~]*Amplification Results.xlsx')  # ignores temp files
    files.sort()
    for file in files:
        (currFam, currHex, currInfoCurve) = read_data_from_curve_file(file)
        currInfoCurve2 = formatInfoData(currInfoCurve)
        filename = currInfoCurve2.at[0, 1]
        currFam2 = formatCurveData(currFam, filename, fluorName="FAM")
        currHex2 = formatCurveData(currHex, filename, fluorName="HEX")
        allCurve.append(pd.concat([currFam2, currHex2]))

    files = glob.glob(path + '[!~]*Cq Results.xlsx')  # ignores temp files
    files.sort()
    for file in files:
        (currCq, currInfoCq) = read_data_from_cq_file(file)
        currInfoCq2 = formatInfoData(currInfoCq)
        filename = currInfoCq2.at[0, 1]
        notes = currInfoCq2.at[2, 1]
        runStart = currInfoCq2.at[4, 1]
        baseSerial = currInfoCq2.at[10, 1]
        currCq2 = formatCqData(currCq, filename, notes, runStart, baseSerial)
        allCq.append(currCq2)

    outfile = path + "output_" + time.strftime("%Y%m%d-%H%M%S", time.localtime()) + ".xlsx"

    cqFinal = allCq.pop(0)
    curveFinal = allCurve.pop(0)

    # for cq in allCq:
    #     cqFinal = pd.concat([cqFinal, cq])
    #
    # for curve in allCurve:
    #     curveFinal = pd.concat([curveFinal, curve])

    for (cq, curve) in zip(allCq, allCurve):
        cqFinal = pd.concat([cqFinal, cq])
        curveFinal = pd.concat([curveFinal, curve])

    # write_data_to_file(outfile, allCq[0], allCurve[0])
    write_data_to_file(outfile, cqFinal, curveFinal)

    print("done!")


if __name__ == "__main__":
    main()
