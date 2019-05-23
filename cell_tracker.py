"""
Filename: cell_tracker.py
Author: Jack Lonnborn
Date: updated May 2019

This program processes raw data from a proteome profiling experiment in cancer research. In particular, it calculates the mean-squared-displacement of migrating cells in micronenvironments characterised by either cancer-associated (CAF) or non-malignant (NPF) fibroblasts. The mean-squared-displacement is (roughly) a measure of the area explored by the cells. The key question this calculation assists in addressing is: do cells in a CAF microenvironment migrate more aggressively than those in an NPF?

The results of this research have been published in:
	Nguyen, E. V. and B. A. Pereria et al, "Proteomic profiling of human prostate cancer-associated fibroblasts (CAF) reveals LOXL2-dependent regulation of the tumor microenvironment" (2019) Molecular & Cellular Proteomics
eprint available at: https://www.mcponline.org/content/early/2019/05/06/mcp.RA119.001496
"""
import pandas as pd
import numpy as np
import itertools
import re
from os import path

class Experiment(object):

	def __init__(self, filename, sheet):
		self.filename = filename
		self.sheet = sheet
		self.data = pd.read_excel(self.filename, sheet_name=self.sheet, header=1, skiprows=0)
		self.times = self.get_unique_vals('Time')
		self.cell_ids = self.get_unique_vals('TrackID')
		self.num_cells = len(self.cell_ids)
		
	def get_unique_vals(self, column):
		"""Returns all unique values in `column`"""
		all_vals = self.data[column].values
		return np.unique(all_vals)

	@property
	def microenvironment(self):
		"""Reads filename to determine if experiment was in a CAF or NPF microenvironment"""
		_, tail = path.split(self.filename)
		if re.search("CAF", tail):
			microenvironment = "CAF"
		elif re.search("NPF", tail):
			microenvironment = "NPF"
		else:
			microenvironment = "Unknown"
		return microenvironment

	def get_sheet_units(self):
		"""Gets the units in which the quantity on `sheet` is reported"""
		unit_vec = np.unique(self.data['Unit'])
		if len(unit_vec)==1:
			return np.unique(unit_vec)[0]
		else:
			print("Warning: units multiply defined in sheet '{}'".format(self.sheet))
			return "Unknown"

	def map_timesTodata(self, column='Value'):
		"""Creates a dictionary with key: time, value: numpy array containing data from `column` which corresponds to that time"""
		data_of_t = {}
		for time in self.times:
			mask = self.data['Time'].values==time
			data_of_t[time] = self.data[column].values[mask]
		return data_of_t

	def calculate_MSD(self):
		ds_data = self.map_timesTodata(column='Value')
		msd_data = []
		for time in ds_data.keys():
			msd_data.append(np.mean(ds_data[time]))
		return msd_data

class Cell(object):
	
	def __init__(self, experiment, cell_id):
		self.cell_id = cell_id
		self.experiment = experiment

	@property
	def times_w_data(self):
		"""Returns the times for which we have data for this cell"""
		times_w_data = [x[0] for x in self.build_data_vec('Time')]
		return times_w_data

	def build_data_vec(self, column):
		"""Returns a list of (time, value) tuples"""
		all_data = self.experiment.data
		mask = all_data['TrackID'].values==self.cell_id
		cell_data = list((t, val) for t, val in zip(all_data['Time'].values[mask], all_data[column].values[mask]))
		return cell_data

# Usage example
######################################
# exp = Experiment("test_data.xls", 'Displacement^2')
# print(exp.microenvironment, exp.num_cells)
# print(exp.times)
# print(exp.cell_ids)
# print(exp.sheet)
# cell = Cell(exp, exp.cell_ids[16])
# print(cell.times_w_data)
# print(cell.build_data_vec('Value'))

# MSD = exp.calculate_MSD()
# units = exp.get_sheet_units()
# import matplotlib.pyplot as plt
# plt.plot(exp.times, MSD)
# plt.xlabel("Time", size=15)
# plt.ylabel(r"Mean-squared-displacement, $\langle x^2 \rangle \quad$ (${}$)".format(units), size=15)
# plt.show()
######################################