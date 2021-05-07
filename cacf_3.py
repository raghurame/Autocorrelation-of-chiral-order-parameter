import statistics 
import numpy as np
import os
from time import sleep
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl
import csv
import decimal
import math as m

def acf(series):
	n = len(series)
	data = np.asarray(series)
	mean = np.mean(data)
	c0 = np.sum((data - mean) ** 2) / float(n)

	def r(h):
		acf_lag = ((data[:n - h] - mean) * (data[h:] - mean)).sum() / float(n) / c0
		return round(acf_lag, 3)

	x = np.arange(n)
	acf_coeffs = list (map(r, x))
	return acf_coeffs

def computeCorrelation (timeSeries, outputFilename, dirs):

	arraySize = len (timeSeries)

	timeSeriesNumpy = np.array (timeSeries)
	time_axisvalues = [x for x in range(0, arraySize)]
	result = acf (timeSeriesNumpy)

	result_trimmed = result[:arraySize]

	rows = zip(time_axisvalues, result_trimmed)
	with open(dirs + "/acf_" + outputFilename + ".csv", "w") as csvfile:
		writer = csv.writer(csvfile)
		for row in rows:
			writer.writerow(row)

	binSize = 0.01
	nBins = m.ceil (2 / binSize)
	binPopulation = {}
	summation = 0

	for bin in range (nBins):
		binMin = -1 + bin * binSize
		binMax = binMin + binSize
		binPopulation[bin] = 0

		for x in timeSeries:
			if (x < binMax and x >= binMin):
				binPopulation[bin] += 1

	for x in range (60, 140, 1):
		summation += binPopulation[x]

	with open (dirs + "/COP_" + outputFilename + ".csv", "w") as output:
		for x, y in binPopulation.items ():
			output.write ("{}, {}\n".format (x, y))


	with open (dirs + "/fixedCOP_" + outputFilename + ".csv", "w") as output:
		for x, y in binPopulation.items ():
			try:
				y *= (4999/summation) # Change the normalization equation based on the input sample
			except:
				pass
			output.write ("{}, {}\n".format (x, y))

	with open (dirs + "/normCOP_" + outputFilename + ".csv", "w") as output:
		for x, y in binPopulation.items ():
			try:
				y /= summation
			except:
				pass
			output.write ("{}, {}\n".format (x, y))

def readfile(filename, dirs):
	chiral_overall = []
	chiral_inside = []
	chiral_outside = []

	with open(filename, "r") as file:
		next(file)
		for line in file:
			separated_values = line.split(",")
			chiral_overall.append(float(separated_values[0]))
			chiral_inside.append(float(separated_values[1]))
			chiral_outside.append(float(separated_values[2]))

	computeCorrelation (chiral_overall[:5000], "overall1", dirs)
	computeCorrelation (chiral_inside[:5000], "inside1", dirs)
	computeCorrelation (chiral_outside[:5000], "outside1", dirs)

	computeCorrelation (chiral_overall[5001:], "overall2", dirs)
	computeCorrelation (chiral_inside[5001:], "inside2", dirs)
	computeCorrelation (chiral_outside[5001:], "outside2", dirs)

def browseDirectories (filetemplate):
	# Goes through all directories and checks for the file "chiral_order_parameter.log"
	for dirs, dirname, files in os.walk("."):
		for file in files:
			if filetemplate in file:
				filePath = dirs + "/" + file
				print ("reading ", dirs + "/" + file)
				readfile (filePath, dirs)

if __name__ == '__main__':
	browseDirectories ("chiral_order_parameter.log")