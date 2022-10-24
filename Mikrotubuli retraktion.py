import os
import sys
from ij import IJ, ImagePlus, WindowManager, ImageStack, Macro
from ij.gui import GenericDialog
from ij.plugin.frame import RoiManager
from ij.gui import Roi


def main():
	# User input directory
	image_dir = IJ.getDirectory("Input_directory")
	
	
	# Make dialog box with channel input
	gd = GenericDialog("Channels")
	gd.addStringField("DAPI in channel: ",'1')
	gd.addStringField("p25a in channel: ",'2')
	gd.addStringField("Tubulin in channel: ",'3')
	gd.showDialog()
	
	dapi_ch, p25a_ch, tubulin_ch = gd.getNextString(), gd.getNextString(), gd.getNextString()

	bog = filesort(image_dir)
	mtr_prog(bog, dapi_ch, p25a_ch, tubulin_ch)
	

# Program to sort files into library with structure - sub directory path : "ORG" tif files in that folder
def filesort(image_dir):
	dir_list = []
	bog = {}
	
	for root, dirs, files in os.walk(image_dir, topdown=True):
		for d in dirs:
			dir_list.append(os.path.join(root,d))
	
	
	for i in dir_list:
		for root, dirs, files in os.walk(i, topdown=True):
			org_files = [x for x in files if x.endswith("ORG.tif")]
			bog[i] = org_files

	return(bog)

def mtr_prog(bog, dapi_ch, p25a_ch, tubulin_ch):


	# Clear ROI manager
	rm = RoiManager.getInstance();
	if (rm==None):
   		rm = RoiManager();
   	rm.reset()


	# Locate all p25a positive cells and add to roi manager
	
	for i in bog:
		for x in bog[i]:
			if "_c{channel}_".format(channel=p25a_ch) in x:
				imp = IJ.openImage(i+'//'+x)
				IJ.setAutoThreshold(imp, "Li dark")
				IJ.run(imp, "Analyze Particles...", "size=1000-100000 pixel show=Masks exclude")
				IJ.run("Invert")
				IJ.run("Create Selection")
				rm.addRoi(imp.getRoi())
				IJ.run("Close")


	# Locate all nuclei and add to roi manager
	# Assume size (mainly to exclude small fragments) and circularity (excludes clumped nuclei)
	
		for x in bog[i]:
			if "_c{channel}_".format(channel=dapi_ch) in x:
				imp = IJ.openImage(i+'//'+x)
				IJ.setAutoThreshold(imp, "Default dark")
				#IJ.run(imp, "Convert to Mask", "")
				#IJ.run(imp, "Watershed", "");
				IJ.run(imp, "Analyze Particles...", "size=300-8000 pixel circularity=0.50-1.00 show=Masks exclude")
				IJ.run("Invert")
				IJ.run("Create Selection")
				rm.addRoi(imp.getRoi())
				IJ.run("Close")

	# Select only rois that overlap (nuclei in p25a positive cells), 
	# clear outside and then select only circular particles (exclude multi nulcei complexes and any nuclei segmented by the roi overlap)

		for x in bog[i]:
			if "_c{channel}_".format(channel=dapi_ch) in x:
				imp = IJ.openImage(i+'//'+x)
				IJ.run(imp, "Convert to Mask","")
				
				if p25a_ch == "":
					rm.select(imp,0)


				else:					
					rm.setSelectedIndexes([0,1])
					rm.runCommand(imp,"AND")

				
				IJ.run(imp, "Clear Outside", "")
				IJ.run(imp, "Select None", "")
				rm.reset()
				IJ.run(imp, "Analyze Particles...", "size=300-8000 pixel circularity=0.50-1.00 show=Masks exclude add")
				IJ.run("Close")


	# Apply nulcie selection of p25 positive cells to tubulin channel and measure median intensity close to and far away nuclei


		n = rm.getCount()
		IJ.run("Set Measurements...", "area mean integrated median display redirect=None decimal=4")
		
		for x in bog[i]:
			if "_c{channel}_".format(channel=tubulin_ch) in x:
				imp = IJ.openImage(i+'//'+x)
				for i in range(n):
					rm.select(imp, i)
					IJ.run(imp, "Enlarge...", "enlarge=15 pixel")
					rm.addRoi(imp.getRoi())
					IJ.run(imp, "Measure","")
					rm.select(imp, n+i*2)
					IJ.run(imp, "Enlarge...", "enlarge=10 pixel")
					rm.addRoi(imp.getRoi())
					rm.setSelectedIndexes([n+i*2,n+i*2+1])
					rm.runCommand(imp,"XOR")
					IJ.run(imp, "Measure","")
					#imp.show()
					rm.deselect()
				rm.reset()

main()
