# cell-tracker

Author: Jack Lonnborn  
Date: updated May 2019

This program processes raw data from a proteome profiling experiment in cancer research. In particular, it calculates the mean-squared-displacement of migrating cells in micronenvironments characterised by either cancer-associated (CAF) or non-malignant (NPF) fibroblasts. The mean-squared-displacement is (roughly) a measure of the area explored by the cells. The key question this calculation assists in addressing is: do cells in a CAF microenvironment migrate more aggressively than those in an NPF?

The results of this research have been published in:  
---
	Nguyen, E. V. and B. A. Pereira et al, "Proteomic profiling of human prostate cancer-associated fibroblasts (CAF) reveals LOXL2-dependent regulation of the tumor microenvironment" (2019) Molecular & Cellular Proteomics  
---  
eprint available at: https://www.mcponline.org/content/early/2019/05/06/mcp.RA119.001496

The program is intended to be extensible and could be modified perform additional calculations on the experimental data.

Usage example:  
Open a terminal and run `python3 -i cell_tracker.py` to start an interactive session.  
```Python
# Create an object for an experiment:
exp = Experiment("test_data.xls")
exp.microenvironment, exp.num_cells
exp.times
exp.cell_ids

# Create an object for one of the tracked cells:
cell = Cell(exp, exp.cell_ids[16])
cell.times_w_data
cell.build_data_vec('Value')

# Calculate and plot the mean-squared-displacement
# (of all cells tracked in the experiment):
msd = exp.calculate_MSD()
import matplotlib.pyplot as plt
plt.plot(exp.times, msd)
plt.xlabel("Time", size=15)
plt.ylabel(r"Mean-squared-displacement, $\langle x^2 \rangle \quad$ ($\mu m^2$)", size=15)
plt.show()
```