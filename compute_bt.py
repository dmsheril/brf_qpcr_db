import sys
import pandas as pd
import numpy as np
import glob
import os
import time


def findBTforCq(curvedata, Cq):
    # find the baseline threshold given a curve and a Cq value
    x1 = np.floor(Cq)
    x2 = np.ceil(Cq)
    y1 = curvedata.iloc[0,int(x1)-1]
    y2 = curvedata.iloc[0,int(x2)-1]
    m = (y2 - y1) / (x2 - x1)
    b = y2 - m * x2
    BT = m * Cq + b
    return BT


def findCqforBT(curvedata, BT):
    # find Cq given a curve and a baseline threshold value
    diffBT = curvedata - BT
    positiveZeroCrossings = np.where(np.diff(np.sign(diffBT)) > 0)[1]

    if len(positiveZeroCrossings) > 0:
        # need to get the index of the curve before the rightmost positive zero crossing
        idxBeforeCrossing = positiveZeroCrossings[-1]
        x1 = idxBeforeCrossing + 1
        x2 = x1 + 1
        y1 = curvedata.iloc[0, int(x1) - 1]
        y2 = curvedata.iloc[0, int(x2) - 1]
        m = (y2 - y1) / (x2 - x1)
        b = y1 - m * x1
        Cq = (BT - b) / m
    else:
        Cq = np.nan
    return Cq


def read_from_processing_output(file):
    print("Reading data from file " + file)
    dataCq = pd.read_excel(file, sheet_name='Cq')
    dataCurve = pd.read_excel(file, sheet_name='Curve')
    return dataCq, dataCurve

def write_data_to_file(outfile, dataCq, dataCurve):
    print("Writing Cq & curve data with BT values to excel file: " + outfile)
    writer = pd.ExcelWriter(outfile, engine='xlsxwriter')
    dataCq.to_excel(writer, sheet_name="Cq", index=False)
    dataCurve.to_excel(writer, sheet_name="Curve", index=False)
    writer.save()

def main():
    path = "/Users/a27inchMac/Desktop/SARS2 qPCR Data/data_to_process/"

    files = glob.glob(path + 'output_*.xlsx')
    files = [file for file in files if "_BT_" not in file] # ignore results of previous processing by this script

    if len(files) < 1:
        print("No files to process in " + path + " - done!")
        sys.exit(0)

    files.sort()

    for j, file in enumerate(files):
        (currCq, currCurve) = read_from_processing_output(file)

        BTvals = np.empty(currCq.shape[0])
        BTvals[:] = np.nan
        CqAtBT120 = np.empty(currCq.shape[0])
        CqAtBT120[:] = np.nan
        CqAtBT200 = np.empty(currCq.shape[0])
        CqAtBT200[:] = np.nan

        for i, rowCq in currCq.iterrows():
            if pd.notna(rowCq['Cq']):
                rowCurve = currCurve.loc[(currCurve["File name"] == rowCq["File name"])
                                     & (currCurve["Well"] == rowCq["Well"])
                                     & (currCurve["Fluor"] == rowCq["Fluor"])].filter(regex='Cycle')
                BTvals[i] = findBTforCq(rowCurve, rowCq["Cq"])
                CqAtBT120[i] = findCqforBT(rowCurve, BT=120)
                CqAtBT200[i] = findCqforBT(rowCurve, BT=200)

                # Cq2 = findCqforBT(rowCurve, BT)
                # print("i = {:d}, Well = {:s}, Fluor = {:s}, Sample = {:s}, Cq = {:f}, #Curves = {:d}, BT = {:f}, "
                #       "Cq2 = {:f} ".format(i, rowCq["Well"], rowCq["Fluor"], rowCq["Sample"], rowCq["Cq"], len(rowCurve), BT, Cq2))

        currCq.drop(labels=["Starting Quantity (SQ)", "Cq Std. Dev", "SQ Std. Dev"], axis="columns", inplace=True)
        currCq.insert(loc=7, column="BT", value=BTvals)
        currCq.insert(loc=8, column="CqAtBT120", value=CqAtBT120)
        currCq.insert(loc=9, column="CqAtBT200", value=CqAtBT200)

        parentdir, fname = os.path.split(file)
        fprefix, ext = os.path.splitext(fname)
        timestr = time.strftime("%Y%m%d-%H%M%S", time.localtime())

        outfile = os.path.join(parentdir, fprefix + "_BT_" + timestr + ext)
        write_data_to_file(outfile, currCq, currCurve)

    print("done!")
    return


if __name__ == "__main__":
    main()
