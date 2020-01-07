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
import regionToolset

def importEBSD(inpFileName,CSVfileName,cohesive):
	fileName, file_extension = os.path.splitext(inpFileName)
	if not CSVfileName or CSVfileName=='<leave empty if the CSV file is in the same folder>':
		CSVfileName=fileName+'.csv'
	while True:
		# Load grain properties
		try:
			file = open(CSVfileName, "r")
			reader = csv.DictReader(file,delimiter='\t',lineterminator='\n')
			grainId=[];phase=[];phi1=[];Phi=[];phi2=[]
			for row in reader:
				grainId.append(int(row['GrainID']),)
				phase.append(row['Phase'],)
				phi1.append(float(row['phi1']),)
				Phi.append(float(row['Phi']),)
				phi2.append(float(row['phi2']),)
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
		elemCohe = mesh.ElemType(elemCode=COH3D8, elemLibrary=STANDARD, elemDeletion=ON, viscosity=0.001)	# Cohesive elements at interfaces
		
		for setName in sets_list:
			elset=sets[setName].elements
			if setName.startswith('GRAIN_'):
				Labels=[j.label for j in elset]
				p.SetFromElementLabels(elementLabels=Labels,name=setName)
				ng+=1
			if setName=='MEDIUM':
				Labels=[j.label for j in elset]
				p.SetFromElementLabels(elementLabels=Labels,name=setName)
		
		# Create materials
		phaseList=set(phase)
		GB_format='Grain Boundary ({})'
		for i in list(phaseList):
			mdb.models[mdbName].Material(name=i)
			mdb.models[mdbName].HomogeneousSolidSection(name=i, material=i,thickness=None)
		if cohesive:
			GB_name='Cohesive boundary'
			mdb.models[mdbName].Material(name=GB_name)
			mdb.models[mdbName].CohesiveSection(name=GB_name,material=GB_name,response=TRACTION_SEPARATION,outOfPlaneThickness=1.0e-3)

		if len(grainId)!=ng:
			print 'Error: the number of grains in INP file and CSV file are inconsistent.'
			break

		# Create grain sections and assign properties
		for i in range(0,ng):
			# Compute the rotation matrix
			R1=np.matrix([[np.cos(phi1[i]),np.sin(phi1[i]),0],[-np.sin(phi1[i]),np.cos(phi1[i]),0],[0,0,1]])
			R2=np.matrix([[1,0,0],[0,np.cos(Phi[i]),np.sin(Phi[i])],[0,-np.sin(Phi[i]),np.cos(Phi[i])]])
			R3=np.matrix([[np.cos(phi2[i]),np.sin(phi2[i]),0],[-np.sin(phi2[i]),np.cos(phi2[i]),0],[0,0,1]])
			M=R3*R2*R1

			# Extract the coordinates of (rotated) X and Y vectors
			X=M[0,:]
			Y=M[1,:]

			# Assign the local coordinates to each section
			sectionID='GRAIN_{:d}'.format(grainId[i])
			region = p.sets[sectionID]
			p.SectionAssignment(region=region, sectionName=phase[i],offset=0.0,offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)
			datumName='ORIENT_{:d}'.format(grainId[i])
			p.DatumCsysByThreePoints(name=datumName, coordSysType=CARTESIAN, origin=(.0,.0,.0), point1=(X[0,0], X[0,1], X[0,2]), point2=(Y[0,0], Y[0,1], Y[0,2]))
			id=p.features[datumName].id
			orientation = p.datums[id]
			p.MaterialOrientation(region=region,orientationType=SYSTEM,axis=AXIS_3,localCsys=orientation,fieldName='',additionalRotationType=ROTATION_NONE, angle=0.0,additionalRotationField='', stackDirection=STACK_3)
		
			if cohesive:
				elset=region.elements
				elemLabels=[k.label for k in elset]
				labels=list()
				faceN=list()
				for p1 in elset:
					f1=p1.getElemFaces()
					for idf in range(0,len(f1)):
						elems=f1[idf].getElements()
						for k in range(0,len(elems)):
							if not(elems[k].label in elemLabels):# and (elems[k].type.name!='COH3D8'):
								corres=elems[k]
								labels.append(corres.label)
								f2=corres.getElemFaces()
								for k2 in range(0,len(f2)):
									if p1 in f2[k2].getElements():
										faceN.append(f2[k2].face)
										break
								break
				strcmd=''
				face=list()
				for j in range(1,7):
					facej=[p.elements.getFromLabel(labels[ii]) for ii in range(len(labels)) if faceN[ii].name == 'FACE{}'.format(j)]
					face.append(mesh.MeshElementArray(facej))
					if facej!=[]:
						strcmd+='face{}Elements=face[{}],'.format(j,j-1)
						
				if strcmd!='':
					newregion=eval('regionToolset.Region(' + strcmd +')')
					nelem=len(p.elements)
					p.generateMeshByOffset(region=newregion, meshType=SOLID, totalThickness=10,numLayers=1, shareNodes=True, offsetDirection=OUTWARD)
					newelems=p.elements[nelem:]
					newregion=(newelems,)
					#setName='GRAIN BOUNDARIES'
					#p.Set(elements=newregion, name=setName)
					p.setElementType(regions=newregion, elemTypes=(elemCohe, ))
					p.SectionAssignment(region=newregion, sectionName=GB_name,offset=0.0,offsetType=MIDDLE_SURFACE, offsetField='',thicknessAssignment=FROM_SECTION)

		# Set the new part visible
		session.viewports['Viewport: 1'].setValues(displayedObject=p)
		
		break