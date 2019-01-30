import string
import csv
import os
import mesh
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from load import *
from mesh import *
from visualization import *
import numpy as np
import re

def importEBSD(inpFileName,CSVfileName):
	fileName, file_extension = os.path.splitext(inpFileName)
	if not CSVfileName or CSVfileName=='<leave empty if the CSV file is in the same folder>':
		CSVfileName=fileName+'.csv'
	while True:
		# Load grain properties
		try:
			file = open(CSVfileName, "r")
			reader = csv.DictReader(file,delimiter='\t',lineterminator='\n',quoting = csv.QUOTE_NONNUMERIC)
			phase=[];phi1=[];Phi=[];phi2=[]
			for row in reader:
				phase.append(row['Phase'],)
				phi1.append(row['phi1'],)
				Phi.append(row['Phi'],)
				phi2.append(row['phi2'],)
			file.close()
		except IOError:
			print 'Error:',CSVfileName,'not found.'
			break

		mdbName=os.path.basename(fileName)
		# Import INP file
		try:
			mdb.ModelFromInputFile(name=mdbName,inputFileName=inpFileName)
			pk=mdb.models[mdbName].parts.keys()
			partName=pk[0]
		except IndexError:
			print 'Error:',fileName+'.inp','not found.'
			break

		a=mdb.models[mdbName].rootAssembly
		sets=a.sets
		sets_list=sets.keys()
		p = mdb.models[mdbName].parts[partName]

		# Copy sets from assembly to part
		ng=0	# Number of grains
		elemSurf = mesh.ElemType(elemCode=S4, elemLibrary=STANDARD, secondOrderAccuracy=OFF)
		elemCohe = mesh.ElemType(elemCode=COH3D8, elemLibrary=STANDARD)	# Cohesive elements at interfaces
		for setName in sets_list:
			elset=sets[setName].elements
			if setName.startswith('GRAIN_') | setName.startswith('GB_'):
				if setName.startswith('GRAIN_'):	# Grain volume
					Labels=[j.label for j in elset]
					p.SetFromElementLabels(elementLabels=Labels,name=setName)
					ng+=1
				elif setName.startswith('GB_'):		# Interface surface
					mask=elset.getMask()
					pickedRegions =(p.elements.getSequenceFromMask(mask), )
					p.setElementType(regions=pickedRegions, elemTypes=(elemSurf, ))
					p.Surface(side1Elements=pickedRegions, name=setName)	# Generate surface before offset
					nelem=len(p.elements)
					p.generateMeshByOffset(region=p.surfaces[setName], meshType=SOLID, totalThickness=0.0, numLayers=1, offsetDirection=OUTWARD, shareNodes=True, deleteBaseElements=False)
					newelems=p.elements[nelem:]		# Indices of new elements
					region=(newelems,)
					p.Set(elements=region, name=setName)
					p.setElementType(regions=region, elemTypes=(elemCohe, ))
					m=re.search(r"^(GB_.+)_\d$",setName)		# Remove trailing number
					if m is not None:
						matName=m.group(1)
						if not(matName in mdb.models['small'].materials.keys()):	# Create new material
							mdb.models[mdbName].Material(name=matName)
							mdb.models[mdbName].CohesiveSection(name=matName,material=matName, response=TRACTION_SEPARATION,outOfPlaneThickness=None)
					p.SectionAssignment(region=region, sectionName=matName, offset=0.0,offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)
					#p.deleteElement(elements=p.elements.getSequenceFromMask(mask), deleteUnreferencedNodes=ON)



		# Set the new part visible
		session.viewports['Viewport: 1'].setValues(displayedObject=p)
		
		# Create materials
		phaseList=set(phase)
		for i in list(phaseList):
			mdb.models[mdbName].Material(name=i)

		if len(phase)!=ng:
			print 'Error: the number of grains in INP file and CSV file are inconsistent.'
			break

		# Create grain sections and assign properties
		for i in range(0,len(phase)):
			# Compute the rotation matrix
			R1=np.matrix([[np.cos(phi1[i]),np.sin(phi1[i]),0],[-np.sin(phi1[i]),np.cos(phi1[i]),0],[0,0,1]])
			R2=np.matrix([[1,0,0],[0,np.cos(Phi[i]),np.sin(Phi[i])],[0,-np.sin(Phi[i]),np.cos(Phi[i])]])
			R3=np.matrix([[np.cos(phi2[i]),np.sin(phi2[i]),0],[-np.sin(phi2[i]),np.cos(phi2[i]),0],[0,0,1]])
			M=R3*R2*R1

			# Extract the coordinates of (rotated) X and Y vectors
			X=M[0,:]
			Y=M[1,:]

			# Assign the local coordinates to each section
			sectionID='GRAIN_{:d}'.format(i+1)
			mdb.models[mdbName].HomogeneousSolidSection(name=sectionID, material=phase[i],thickness=None)
			region = p.sets[sectionID]
			p.SectionAssignment(region=region, sectionName=sectionID,offset=0.0,offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)
			datumName='ORIENT_{:d}'.format(i+1)
			p.DatumCsysByThreePoints(name=datumName, coordSysType=CARTESIAN, origin=(.0,.0,.0), point1=(X[0,0], X[0,1], X[0,2]), point2=(Y[0,0], Y[0,1], Y[0,2]))
			id=p.features[datumName].id
			orientation = p.datums[id]
			p.MaterialOrientation(region=region,orientationType=SYSTEM,axis=AXIS_3,localCsys=orientation,fieldName='',additionalRotationType=ROTATION_NONE, angle=0.0,additionalRotationField='', stackDirection=STACK_3)
		
				
		break