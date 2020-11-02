import pandas as pd
import glob
import time
import os
import sys
from xlrd import XLRDError

print("This is the top of the process_qpcr module")


def read_data_from_curve_file(file):
    print("Reading data from curve file " + file)
    try:
        dataFam = pd.read_excel(file, sheet_name='FAM')
    except XLRDError:
        print("Sheet FAM not found in file: " + file)
        dataFam = pd.DataFrame()

    try:
        dataHex = pd.read_excel(file, sheet_name='HEX')
    except XLRDError:
        print("Sheet HEX not found in file: " + file)
        dataHex = pd.DataFrame()

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
    if not curveData.empty:
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


def doFileCleanup(sourcepath, destpath):
    filesToKeep = glob.glob(os.path.join(sourcepath, '*Amplification Results.xlsx')) \
                  + glob.glob(os.path.join(sourcepath, '*Cq Results.xlsx'))

    if len(filesToKeep) > 0:
        try:
            os.makedirs(destpath)
        except OSError as error:
            print(error)

        for file in filesToKeep:
            head, tail = os.path.split(file)
            os.rename(file, os.path.join(destpath, tail))
            print("Moved " + tail + " to " + destpath)

    filesToDelete = glob.glob(os.path.join(sourcepath, '*ANOVA Results.xlsx')) \
                    + glob.glob(os.path.join(sourcepath, '*End Point Results.xlsx')) \
                    + glob.glob(os.path.join(sourcepath, '*Bar Chart.xlsx')) \
                    + glob.glob(os.path.join(sourcepath, '*Melt Curve Plate View Results.xlsx')) \
                    + glob.glob(os.path.join(sourcepath, '*Quantification Plate View Results.xlsx')) \
                    + glob.glob(os.path.join(sourcepath, '*Quantification Summary.xlsx')) \
                    + glob.glob(os.path.join(sourcepath, '*Allelic Discrimination Results.xlsx')) \
                    + glob.glob(os.path.join(sourcepath, '*Standard Curve Results.xlsx'))

    for file in filesToDelete:
        os.remove(file)
        head, tail = os.path.split(file)
        print("Deleted " + tail + " from " + head)

    return filesToKeep, filesToDelete


def main():
    # path = "/Users/a27inchMac/Desktop/SARS2 qPCR Data/data_to_process/"
    path = "/Users/a27inchMac/Desktop/Delsey/SARS2/data_proc_test"

    filesAmp = glob.glob(os.path.join(path, '[!~]*Amplification Results.xlsx'))  # ignores temp files
    filesCq = glob.glob(os.path.join(path, '[!~]*Cq Results.xlsx'))  # ignores temp files

    if (len(filesAmp) < 1) | (len(filesCq) < 1):
        print("No files to process in " + path + " - done!")
        sys.exit(0)

    assert (len(filesAmp) == len(filesCq)), "Number of files for Amp and Cq must match - check for missing files"

    filesAmp.sort()
    filesCq.sort()
    allCurve = []
    allCq = []

    for file in filesAmp:
        (currFam, currHex, currInfoCurve) = read_data_from_curve_file(file)
        currInfoCurve2 = formatInfoData(currInfoCurve)
        filename = currInfoCurve2.at[0, 1]
        currFam2 = formatCurveData(currFam, filename, fluorName="FAM")
        currHex2 = formatCurveData(currHex, filename, fluorName="HEX")
        allCurve.append(pd.concat([currFam2, currHex2]))

    for file in filesCq:
        (currCq, currInfoCq) = read_data_from_cq_file(file)
        currInfoCq2 = formatInfoData(currInfoCq)
        filename = currInfoCq2.at[0, 1]
        notes = currInfoCq2.at[2, 1]
        runStart = currInfoCq2.at[4, 1]
        baseSerial = currInfoCq2.at[10, 1]
        currCq2 = formatCqData(currCq, filename, notes, runStart, baseSerial)
        allCq.append(currCq2)

    timestr = time.strftime("%Y%m%d-%H%M%S", time.localtime())
    outfile = os.path.join(path, "output_" + timestr + ".xlsx")
    procdir = os.path.join(path, "processed_" + timestr)

    cqFinal = allCq.pop(0)
    curveFinal = allCurve.pop(0)
    for (cq, curve) in zip(allCq, allCurve):
        cqFinal = pd.concat([cqFinal, cq])
        curveFinal = pd.concat([curveFinal, curve])

    write_data_to_file(outfile, cqFinal, curveFinal)
    doFileCleanup(path, procdir)
    print("done!")


if __name__ == "__main__":
    main()
