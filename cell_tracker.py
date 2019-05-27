"""
Filename: cell_tracker-2.py
Version: 2.0
Author: Jack Lonnborn
Date: updated May 2019

This program processes raw data from a proteome profiling experiment in cancer research. In particular, it calculates the mean-squared-displacement of migrating cells in micronenvironments characterised by either cancer-associated (CAF) or non-malignant (NPF) fibroblasts. The mean-squared-displacement is (roughly) a measure of the area explored by the cells. The key question this calculation assists in addressing is: do cells in a CAF microenvironment migrate more aggressively than those in an NPF?

The results of this research have been published in:
	Nguyen, E. V. and B. A. Pereria et al, "Proteomic profiling of human prostate cancer-associated fibroblasts (CAF) reveals LOXL2-dependent regulation of the tumor microenvironment" (2019) Molecular & Cellular Proteomics
eprint available at: https://www.mcponline.org/content/early/2019/05/06/mcp.RA119.001496

The experimental data comes from time-lapse videos of samples which are analysed using propeitary object tracking software by Imaris. More information: https://imaris.oxinst.com/products/imaris-for-tracking
"""

import pandas as pd
import numpy as np
import re
from os import path
import sys

class Sheet(object):

	def __init__(self, filename, sheet):
		self.filename = filename
		self.sheet = sheet
		if re.search("Track", self.sheet):
			# cells are labelled by `ID` or `TrackID` depending on the sheet
			self.id_label = 'ID'
		else:
			self.id_label = 'TrackID'
		self.data = pd.read_excel(self.filename, sheet_name=self.sheet, header=1, skiprows=0)
		self.times = self.get_unique_vals('Time')
		self.cell_ids = self.get_unique_vals(self.id_label)
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

	@property
	def num_data_cols(self):
		return self.data.columns.get_loc('Unit')
	
	def get_sheet_units(self):
		"""Gets the units in which the quantity on `sheet` is reported"""
		unit_vec = np.unique(self.data['Unit'])
		if len(unit_vec)==1:
			return np.unique(unit_vec)[0]
		else:
			print("Warning: units multiply defined in sheet '{}'".format(self.sheet))
			return "Unknown"

	def data_array(self):
		"""Returns a len(times) numpy array, each element of which is a tuple containing
		(time, x) where x is a numpy array containing values from the data columns of `sheet` which correspond to `time`."""
		data_array = np.empty(len(self.times), dtype=object)
		for i, time in enumerate(self.times):
			mask = self.data['Time'].values==time
			data_array[i] = (time, self.data.iloc[:,:self.num_data_cols].values[mask])
		return data_array
		
	def calculate_MSD(self):
		if self.sheet != 'Displacement^2':
			sys.exit("Initialise Sheet object with sheet=`Displacement^2` to calculate mean-squared-displacement!")
		ds_data = self.data_array()
		msd_data = []
		for row in ds_data:
			msd_data.append(np.mean(row[1]))
		return msd_data

class Cell(object):
	
	def __init__(self, experiment, cell_id):
		self.cell_id = cell_id
		self.experiment = experiment

	@property
	def times_w_data(self):
		"""Returns the times for which we have data for this cell"""
		mask = self.experiment.data[self.experiment.id_label].values==self.cell_id
		times_w_data = np.sort(self.experiment.data['Time'][mask].values)
		return times_w_data

	def data_array(self):
		"""Returns a len(times) numpy array, each element of which is a tuple containing
		(time, x) where x is a numpy array containing values from the data columns of `sheet` which correspond to `time` for this cell."""
		data_array = np.empty(len(self.times_w_data), dtype=object)
		exp = self.experiment
		cell_mask = exp.data[exp.id_label].values==self.cell_id
		for i, time in enumerate(self.times_w_data):
			time_mask = exp.data['Time'].values==time
			full_mask = [all(mask) for mask in zip(time_mask, cell_mask)]
			data_columns = exp.data.iloc[:,:exp.num_data_cols].values
			data_array[i] = (time, data_columns[full_mask][0])
		return data_array