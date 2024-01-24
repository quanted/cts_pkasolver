import logging
import pandas as pd
from rdkit import Chem
from os import path
import pkasolver
from pkasolver.query import calculate_microstate_pka_values


# Loads trained model:
base_path=path.dirname(pkasolver.__file__)


class CTSPkasolver:
	"""
	Function calls to pkasolver used by CTS.
	"""

	def __init__(self):
		self.pka_dec= 2
		self.step=0.2  # step size for charts
		self.minph=0
		self.maxph=14

	def format_chart_data(self, chart_data):
		"""
		Formats cts-pkasolver chart_data into something the
		speciation output page is expecting.

		Reformats this:
		{
			"x": x,
			"a0": a0,
			"a1": a1,
			"a2": a2,
			"a3": a3
		}
		Into this:
		{
			"microspecies1": [[x[0], a0[0]], [x[1], a0[1]], ...],
			"microspecies2": [[x[0], a1[0]], ...] 
		}

		"""
		results_obj = {}
		# chart_data = results_obj.get("chart_data", {})

		for key in chart_data:
			if key != "x":
				result = [[chart_data["x"][i], (100.0 * chart_data[key][i])] for i in range(len(chart_data[key]))]
				new_key = "microspecies{}".format(int(key[1]) + 1)
				results_obj[new_key] = result


		return results_obj

	def run_pka_solver(self, smiles):
		"""
		Generator function that returns pKa values and more.
		"""
		mol=Chem.MolFromSmiles(smiles)  # creates mol object from smiles
		protonation_states = calculate_microstate_pka_values(mol) # performs internal calculations and stores as object
		sites=len(protonation_states) # gets the number of ionization sites 
		
		lst=[]
		depSmi=[]
		proSmi=[]
		
		for j in range(len(protonation_states)):
			state=protonation_states[j]
			depSmi.append(Chem.MolToSmiles(state.deprotonated_mol))
			proSmi.append(Chem.MolToSmiles(state.protonated_mol))
			lst.append(round(state.pka,2)) # get pka values for all sites for a given molecule store in a list

		yield sites, lst, proSmi, depSmi

	def get_mono_plot(self, pka_list, proSmi, depSmi):
		pka = pka_list[0]
		x = []
		a0 = []
		a1 = []
		
		ph = self.minph
		
		while ph <= self.maxph:
			ph += self.step
			x.append(ph)
		
			a0.append(round((1/(1+10**(ph-pka))),self.pka_dec))
			a1.append(1-(round((1/(1+10**(ph-pka))),self.pka_dec))) # or round((10**(ph-pka))/(1+10**(ph-pka)),self.pka_dec)

		chart_data = {
			"x": x,
			"a0": a0,
			"a1": a1
		}

		
		# Makes a list of smiles in order 
		smiles=[proSmi[0]]+depSmi
		
		# Makes a dictionary where the keys are the a_index (ex.0= a0, 1=a1, etc.) and values are the microspecies smiles strings
		microspecies=dict(list(enumerate(smiles)))

		return chart_data, microspecies

	def get_di_plot(self, pka_list, proSmi, depSmi):
		pka1 = pka_list[0]
		pka2 = pka_list[1]
		
		x = []
		a0 = []
		a1 = []
		a2 = []
		
		ph = self.minph
		
		while ph <= self.maxph:
			ph += self.step
			x.append(ph)
		
			ka1 = 10**(ph-pka1)
			ka2 = 10**(ph-pka2)
			E = ((1+ka1) + (ka1*ka2))
		
			a0.append(round(((1**2)/E),self.pka_dec))
			a1.append(round(((1*ka1)/E),self.pka_dec))
			a2.append(round(((ka1*ka2)/E),self.pka_dec))


		chart_data = {
			"x": x,
			"a0": a0,
			"a1": a1,
			"a2": a2
		}
		
		# Makes a list of smiles in order 
		smiles=[proSmi[0]] + depSmi
		
		# Makes a dictionary where the keys are the a_index (ex.0= a0, 1=a1, etc.) and values are the microspecies smiles strings
		microspecies=dict(list(enumerate(smiles)))

		return chart_data, microspecies

	def get_tri_plot(self, pka_list, proSmi, depSmi):
		pka1 = pka_list[0]
		pka2 = pka_list[1]
		pka3 = pka_list[2]
		
		x = []
		a0 = []
		a1 = []
		a2 = []
		a3 = []
		
		ph = self.minph
		
		while ph <= self.maxph:
			ph += self.step
			x.append(ph)
		
			ka1 = 10**(ph-pka1)
			ka2 = 10**(ph-pka2)
			ka3 = 10**(ph-pka3)
			D = ((1+ka1)+(ka1*ka2)+(ka1*ka2*ka3))
		
			a0.append(round(((1**2)/D),self.pka_dec))
			a1.append(round(((1*ka1)/D),self.pka_dec))
			a2.append(round(((ka1*ka2)/D),self.pka_dec))
			a3.append(round(((ka1*ka2*ka3)/D),self.pka_dec))

		chart_data = {
			"x": x,
			"a0": a0,
			"a1": a1,
			"a2": a2,
			"a3": a3
		}

		# Makes a list of smiles in order 
		smiles = [proSmi[0]] + depSmi

		# Makes a dictionary where the keys are the a_index (ex.0= a0, 1=a1, etc.) and values are the microspecies smiles strings
		microspecies=dict(list(enumerate(smiles)))

		return chart_data, microspecies

	def get_multi_plot_orig(self, pka_sites, pka_list, proSmi, depSmi):
		df = pd.DataFrame()
		x = []
		a0 = []
		ax = []
		
		ph = self.minph
		
		while ph <= self.maxph:
			ph += self.step
			x.append(ph)
			count = 0
			D = 1
			numTerms=[]
			for i in range(0,(pka_sites-1)):
				n1 = 10**(ph-pka_list[i])
				n2 = 10**(ph-pka_list[1+i])
				#if there is only one term, return D+n1 (should not happend)
				if pka_sites == 1:
					D += n1
				else:
					while count < pka_sites:
						#get the numerator term and save to list
						N = n1
						numTerms.append(N)
						#caluclate denominator
						D += n1
						#calculate next term
						nth = n1 *n2
						#update values
						n1 = nth
						n2 = 10**(ph-pka_list[i+2])
						count += 1

			# Calculates ionization fraction for each numTerm
			# a0.append(round(((1**2)/D), self.pka_dec))
			a = []
			for t in numTerms:
				a.append(round((t/D), self.pka_dec))
			ax.append(a)
		   
		df['pH'] = x
		df['ax'] = ax
		# df['a0'] = a0

		# Separate out ionization fractions into their own columns (based on ka)
		points = pd.DataFrame(df.ax.tolist()).add_prefix('a')
		
		# Combine pka columns with rest of data
		data = pd.concat([df,points],axis=1)

		# Makes a list of smiles in order 
		smiles = [proSmi[0]]+depSmi
		
		# Makes a dictionary where the keys are the a_index (ex.0= a0, 1=a1, etc.) and values are the microspecies smiles strings
		microspecies = dict(list(enumerate(smiles)))

		chart_data = {}

		# Creates chart data for speciation workflow:
		for i in data.columns[2:]:
			chart_data[i] = data[i].tolist()

		chart_data['x'] = data['pH'].tolist()
			
		return chart_data, microspecies

	def get_multi_plot(self, pka_sites,pka_list,proSmi,depSmi):
		df=pd.DataFrame()
		x=[]
		a0=[]
		ax=[]
		
		ph=self.minph
		
		while ph<=self.maxph:
			ph+=self.step
			x.append(ph)
			count=0
			D=1
			numTerms=[]
			for i in range(0,(pka_sites-1)):
				n1=10**(ph-pka_list[i])
				n2=10**(ph-pka_list[1+i])
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
						n2=10**(ph-pka_list[i+2])
						count+=1

			#calculate ionization fraction for each numTerm
			a0.append(round(((1**2)/D),self.pka_dec))
			a=[]
			for t in numTerms:
				a.append(round((t/D),self.pka_dec))
			ax.append(a)
		   
		df['pH']=x
		df['ax']=ax
		df['a0']=a0
		
		#make list of column names 
		col_names=[]
		for i in range(1,len(df.ax[0])+1):
			name='a'+str(i)
			col_names.append(name)
				
		#separate out ionization fractions into their own columns (based on ka)
		points=pd.DataFrame(df.ax.tolist())

		#assign columns appropriate name    
		points.columns=col_names
		
		#combine pka columns with rest of data
		data=pd.concat([df,points],axis=1)
	 

		#make a list of smiles in order 
		smiles=[proSmi[0]]+depSmi
		#make a dictionary where the keys are the a_index (ex.0= a0, 1=a1, etc.) and values are the microspecies smiles strings
		microspecies=dict(list(enumerate(smiles)))

		
		# Creates chart data for speciation workflow:
		chart_data = {}
		for i in data.columns[2:]:
			chart_data[i] = data[i].tolist()
		chart_data['x'] = data['pH'].tolist()
		
		return chart_data, microspecies





	def main(self, parent, data_type=None):
		"""
		Main function for returning pkas and/or microspecies chart data.
		Examples: ['CC(O)=O','CC(C)C(N)C(O)=O','C(O)1=CC=C(N)C=C1','NC(CCS)C(O)=O','NC(CC1=CN=CN1)C(O)=O']
		"""

		# Calculates pka for input chemical:
		data = self.run_pka_solver(parent)

		# Returns number of pka sites(n), list of pka values(p),
		# and lists of protonated (pro) and deprotonated (dep) microspecies smiles:
		for n, p, ps, ds in data:
			pkaSites = n  # num of sites
			pka_list = p  # list of pka values
			pro = ps  # deprotonated_mol
			dep = ds  # protonated_mol

		# Option to just return pKa list:
		if data_type == "pka":
			return pka_list

		chart_data, species = None, None

		if pkaSites == 1:
			# print("Calling get_mono_plot.")
			chart_data, species = self.get_mono_plot(p, pro, dep)
		elif pkaSites == 2:
			# print("Calling get_di_plot.")
			chart_data, species = self.get_di_plot(p, pro, dep)
		elif pkaSites == 3:
			# print("Calling get_tri_plot.")
			chart_data, species = self.get_tri_plot(p, pro, dep)
		elif pkaSites > 3:
			# print("Calling get_multi_plot.")
			chart_data, species = self.get_multi_plot(n, p, pro, dep)

		reformatted_chart_data = self.format_chart_data(chart_data)

		return reformatted_chart_data, species, pka_list


if __name__ == "__main__":


	# test_smiles = "CC(=O)OC1=CC=CC=C1C(O)=O"
	test_smiles = "OC(=O)CC(O)(CC(O)=O)C(O)=O"

	cts_pkasolver = CTSPkasolver()
	chart_data, species, pka_list = cts_pkasolver.main(test_smiles)
