#import dependancies 
import pickle
import pandas as pd
from rdkit import Chem
import torch
from os import path
import matplotlib.pyplot as plt

#import pkasolver
import pkasolver
from pkasolver.query import QueryModel
from pkasolver.ml_architecture import GINPairV1
from pkasolver.query import draw_pka_map 
from pkasolver.query import calculate_microstate_pka_values, draw_pka_reactions
from IPython.display import display
from IPython.display import HTML

#load trained model
base_path=path.dirname(pkasolver.__file__)

###Functions
def RunPkasolver(x): #takes input of smile 
	mol=Chem.MolFromSmiles(x)
	protonation_states = calculate_microstate_pka_values(mol) #performs internal calculations and stores as object
	sites=len(protonation_states) #get the number of ionization sites 
	
	lst=[]
	depSmi=[]
	proSmi=[]
	
	for j in range(len(protonation_states)):
		state=protonation_states[j]
		depSmi.append(Chem.MolToSmiles(state.deprotonated_mol))
		proSmi.append(Chem.MolToSmiles(state.protonated_mol))
		lst.append(round(state.pka,2)) #get pka values for all sites for a given molecule store in a list
	yield sites, lst,proSmi,depSmi

def GetMonoPlot(pka_list,minph,maxph,step,proSmi,depSmi):
	pka=pka_list[0]
	x=[]
	a0=[]
	a1=[]

	ph=minph

	while ph<=maxph:
		ph+=step
		x.append(ph)

		a0.append(round((1/(1+10**(ph-pka))),pka_dec))
		a1.append(1-(round((1/(1+10**(ph-pka))),pka_dec))) # or round((10**(ph-pka))/(1+10**(ph-pka)),pka_dec)

	plt.plot(x,a0,color='red',label='a0')
	plt.plot(x,a1, color='green',label='a1')
	plt.legend()
	#make a list of smiles in order 
	smiles=[pro[0]]+dep

	print("GetMonoPlot smiles: {}".format(smiles))

	#make a dictionary where the keys are the a_index (ex.0= a0, 1=a1, etc.) and values are the microspecies smiles strings
	microspecies=dict(list(enumerate(smiles)))

	print("GetMonoPlot microspecies: {}".format(microspecies))

	print("GetMonoPlot a0: {}".format(a0))
	print("GetMonoPlot a1: {}".format(a1))

	return plt,microspecies


def GetDiPlot(pka_list,minph,maxph,step,proSmi,depSmi):
	pka1=pka_list[0]
	pka2=pka_list[1]
	
	x=[]
	a0=[]
	a1=[]
	a2=[]
	
	ph=minph
	
	while ph<=maxph:
		ph+=step
		x.append(ph)
	
		ka1=10**(ph-pka1)
		ka2=10**(ph-pka2)
		E=((1+ka1)+(ka1*ka2))
	
		a0.append(round(((1**2)/E),pka_dec))
		a1.append(round(((1*ka1)/E),pka_dec))
		a2.append(round(((ka1*ka2)/E),pka_dec))
	
	plt.plot(x,a0,color='red',label='a0')
	plt.plot(x,a1, color='green',label='a1')
	plt.plot(x,a2,color='blue',label='a2')
	plt.legend()
	
	#make a list of smiles in order 
	smiles=[pro[0]]+dep
	#make a dictionary where the keys are the a_index (ex.0= a0, 1=a1, etc.) and values are the microspecies smiles strings
	microspecies=dict(list(enumerate(smiles)))

	return plt,microspecies

def GetTriPlot(pka_list,minph,maxph,step,proSmi,depSmi):
	pka1=pka_list[0]
	pka2=pka_list[1]
	pka3=pka_list[2]
	
	x=[]
	a0=[]
	a1=[]
	a2=[]
	a3=[]
	
	ph=minph
	
	while ph<=maxph:
		ph+=step
		x.append(ph)
	
		ka1=10**(ph-pka1)
		ka2=10**(ph-pka2)
		ka3=10**(ph-pka3)
		D=((1+ka1)+(ka1*ka2)+(ka1*ka2*ka3))
	
		a0.append(round(((1**2)/D),pka_dec))
		a1.append(round(((1*ka1)/D),pka_dec))
		a2.append(round(((ka1*ka2)/D),pka_dec))
		a3.append(round(((ka1*ka2*ka3)/D),pka_dec))
	
	plt.plot(x,a0,color='red',label='a0')
	plt.plot(x,a1, color='green',label='a1')
	plt.plot(x,a2,color='blue',label='a2')
	plt.plot(x,a3,color='yellow',label='a3')

	#make a list of smiles in order 
	smiles=[pro[0]]+dep
	#make a dictionary where the keys are the a_index (ex.0= a0, 1=a1, etc.) and values are the microspecies smiles strings
	microspecies=dict(list(enumerate(smiles)))

	return plt,microspecies

def GetMultiPlot(pka_sites,pka_list,minph,maxph,step,proSmi,depSmi):
	df=pd.DataFrame()
	x=[]
	a0=[]
	ax=[]
	
	ph=minph
	
	while ph<=maxph:
		ph+=step
		x.append(ph)
		count=0
		D=1
		numTerms=[]
		for i in range(0,(pka_sites-1)):
			n1=10**(ph-pka_lst[i])
			n2=10**(ph-pka_lst[1+i])
			#if there is only one term, return D+n1 (should not happend)
			if pka_sites ==1:
				D+=n1
			else:
				while count < pka_sites:
					#get the numerator term and save to list
					N=n1
					numTerms.append(N)
					#caluclate denominator
					D+=n1
					#calculate next term
					nth=n1 *n2
					#update values
					n1=nth
					n2=10**(ph-pka_lst[i+2])
					count+=1

		#calculate ionization fraction for each numTerm
		a0.append(round(((1**2)/D),pka_dec))
		a=[]
		for t in numTerms:
			a.append(round((t/D),pka_dec))
		ax.append(a)
	   
	df['pH']=x
	df['ax']=ax
	df['a0']=a0
	
	
	#separate out ionization fractions into their own columns (based on ka)
	points=pd.DataFrame(df.ax.tolist()).add_prefix('a')
	
	#combine pka columns with rest of data
	data=pd.concat([df,points],axis=1)

	#make a list of smiles in order 
	smiles=[pro[0]]+dep
	#make a dictionary where the keys are the a_index (ex.0= a0, 1=a1, etc.) and values are the microspecies smiles strings
	microspecies=dict(list(enumerate(smiles)))

	
	#plot
	for i in data.columns[2:]:
		x=data['pH']
		y=data[i]
		plt.plot(x,y,lable=i)
		plt.legend()
	return plt, microspecies

### DEMO
## user input
pka_dec= 2 
step=0.2 
minph=int(0)
maxph=int(14)

# parent='CC(=O)OC1=CC=CC=C1C(O)=O'
parent = 'OC(=O)CC(O)(CC(O)=O)C(O)=O'
#calculating pka for input chemical
data=RunPkasolver(parent)

#function return number of pka sites(n), list of pka values(p), and lists of protonated (pro) and deprotonated (dep) microspecies smiles
for n,p,ps,ds in data:
	pkaSites=n
	pka_list=p
	pro=ps
	dep=ds

print("pkaSites: {}\npka_list: {}\nprotonated sites: {}\ndeprotonated: {}".format(pkaSites, pka_list, pro, dep))

#determine which set of equations to use based on number of pka sites
if pkaSites == 1:
	#user monoprotic function
	graph,species=GetMonoPlot(p,minph,maxph,step,pro,dep)
elif pkaSites == 2:
	#use diprotic funcion
	graph,species=GetDiPlot(p,minph,maxph,step,pro,dep)
elif pkaSites == 3:
	#use triprotic function
	graph,species=GetTriPlot(p,minph,maxph,step,pro,dep)
elif pkaSites > 3:
	#use generalized function
	graph,species=GetMultiPlot(n,p,minph,maxph,step,pro,dep)