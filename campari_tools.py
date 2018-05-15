import sys
import MDAnalysis as mda
import MDAnalysis.analysis.distances as mdad
import subprocess
import numpy as np

pdb = sys.argv[1]
campari_home = '/home/dillion/pkg/campari/'

class CAMPARI:

	def __init__(self, pdb, pbc):
		self.pdb = pdb
		self.pbc = pbc
		self.name = pdb.split('.')[0]  

	def check_caps(self):
		"""
		check to see if protein is capped with ACE and NME

		"""

		uni = mda.Universe(pdb)
		sel = uni.select_atoms('protein')
		if sel.atoms.resnames[0] != 'ACE':
			return 1
		else:
			return 0
	
	def merge_structures(self,ace,ref,nme):
		"""
		MDAnalysis is very slopy with assigning resids. The Merge function
		by default will not reassign resids, so this function takes care of it
	
		"""
		
		a = ace.select_atoms('resname ACE') ; al = len(a.atoms.residues)
		r = ref.select_atoms('all')	    ; rl = len(r.atoms.residues)
		n = nme.select_atoms('resname NME') ; nl = len(n.atoms.residues)
		
		a.atoms.residues[0].resid = 1
		
		curr_res = 2
		for i in range(rl):
			r.atoms.residues[i].resid = curr_res
			curr_res += 1
		
		n.atoms.residues[0].resid = curr_res
		
		return mda.Merge(a.atoms,r.atoms,n.atoms)

	def cap_protein(self):
		"""
		CAMPARI requires all proteins to be capped with ACE and NME. Most don't have that. Add it if needed.

		"""

		from MDAnalysis.analysis.align import alignto
	
		ref = mda.Universe(pdb)
		ace = mda.Universe("ace.pdb")
		nme = mda.Universe("nme.pdb")
	
		resids = ref.select_atoms("all").resids
		resid_min, resid_max = min(resids), max(resids)
		
		alignto(ace, ref, select={"mobile": "resname ALA and backbone", 
					  "reference": "resid {0} and backbone".format(resid_min)})
		alignto(nme, ref, select={"mobile": "resname ALA and backbone",
					  "reference": "resid {0} and backbone".format(resid_max)})

		A = ace.select_atoms('resname ACE')
		N = nme.select_atoms('resname NME')

		uni = self.merge_structures(ace,ref,nme)
		self.name += '_capped'
		uni.atoms.write(self.name+'.pdb')

		return uni

	def rename_HIS(self):
		"""
		Rosetta names all Histidines 'HIS', but really there are 2 kinds of histidines, and Rosetta 
		makes which ever one it thinks is best. This function detects the type and renames the
		HIS residues with the appropriate names

		"""

		uni = mda.Universe(self.name+'.pdb')
		sel = uni.select_atoms('resname HIS and name CA')
		for r in sel.resids:
			HD1 = uni.select_atoms('resid '+str(r)+' and name HD1')

			ND1 = uni.select_atoms('resid '+str(r)+' and name ND1')
			d1 = mdad.distance_array(HD1.positions,ND1.positions)[0][0]

			NE2 = uni.select_atoms('resid '+str(r)+' and name NE2')
			d2 = mdad.distance_array(HD1.positions,NE2.positions)[0][0]
			
			his = uni.select_atoms('resid ' + str(r))
			if d1 < d2:
				his.residues.resnames = 'HID'
				print "changed residue"
			else:
				his.residues.resnames = 'HIE'
				print "changed residue"

		self.name += '_HIScheck'
		uni.atoms.write(self.name+'.pdb')

		return None

	def write_pdb(self,coors,atom):
		"""
		MDAnalysis doesnt allow users to add atoms, so you have to make a pdb and load it

		"""

		charge = {'NA':'+1', 'CL':'-1'}
		ion_type = {'NA':'NA+', 'CL':'CL-'}  

		outfile = open(self.name+'_charge.pdb',"w")
		it = 0
		# loop contains extra residue at beginning and extra residue at end. Don't put those in PDB
		for i in coors:
			t1 = "ATOM"					# ATOM
			t2 = it						# INDEX
			t3 = atom					# ATOM NAME
			t4 = ""						# ALTERNATE LOCATION INDICATOR
			t5 = ion_type[atom]				# RESIDUE NAME
			t6 = "A"					# CHAIN
			t7 = it						# RESIDUE NUMBER
			t8 = ""						# INSERTION CODE
			t9 =  float(i[0])				# X
			t10 = float(i[1])				# Y
			t11 = float(i[2])				# Z
			t12 = 0.0					# OCCUPANCY
			t13 = 0.0					# TEMPERATURE FACTOR
			t14 = ""					# ELEMENT SYMBOL
			t15 = str(charge[atom])				# CHARGE ON ATOM
			outfile.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15))
			it+=1
		outfile.close()
		return 0

	def get_charge(self, seq):
		"""
		Sum positive and negative charges
	
		"""
	
		positive_charge = ['ARG', 'HID', 'HIS', 'HIE', 'LYS']
		negative_charge = ['ASP', 'GLU'] 
		plus = 0 ; minus = 0
		for aa in positive_charge:
			plus+=len(np.where(seq==aa)[0])
			print aa, len(np.where(seq==aa)[0])
		for aa in negative_charge:
			minus+=len(np.where(seq==aa)[0])
			print aa, -1*len(np.where(seq==aa)[0])
		return plus-minus-1

	def add_charge(self, sel, ion, charge):
		"""
		Generate coordinates for ions and write them to a pdb

		"""

		pos = []
		for c in range(charge):
			pos.append(np.random.rand(3)*self.pbc)
		self.write_pdb(pos,ion)

		return None

	def ionize(self):
		"""
		This is a very simple function that adds counter ions to neutralize the charge of
		the system.

		"""

		uni = mda.Universe(self.name+'.pdb')
		sel = uni.select_atoms('name CA')
		net_charge = self.get_charge(sel.resnames)

		print net_charge

		if net_charge == 0:
			return None
		elif net_charge > 0:
			self.add_charge(sel, 'CL', net_charge)
		elif net_charge < 0:
			self.add_charge(sel, 'NA', net_charge)
		else:
			print "something is wrong. net charge is", net_charge; exit()

		chg = mda.Universe(self.name+'_charge.pdb')
		ions = chg.select_atoms('all')
		protein = uni.select_atoms('all')

		dist = mdad.distance_array(protein.positions,ions.positions)

		if dist.min() < 2:
			print "OVERLAP!!", dist.min()
			exit()

		c = chg.select_atoms('all') ; cl = len(c.atoms.residues)
		p = uni.select_atoms('all') ; pl = len(p.atoms.residues)
		
		for i in range(cl):
			c.atoms.residues[i].resid = i+pl
		
		merged = mda.Merge(p.atoms,c.atoms)
		print merged
		self.name += '_ionized'
		merged.atoms.write(self.name+'.pdb')
		return None
	
	def pdb2seq(self):
		"""
		Generate the CAMPARI specific .seq file

		"""

		uni = mda.Universe(self.name+'.pdb')
		sel = uni.select_atoms('all')
		seqfile = open(self.name+'.seq','w')
		old_r = -1
		for c in sel:
			r = c.resnum
			if old_r == r:
				continue
			else:
				seqfile.write(c.resname.lower()+"\n")
			old_r = r
	
		seqfile.write('END')
		seqfile.close()
	
		return None
	
	def make_input(self):
		"""
		Generate an input script for REMC 

		"""

		input_file = open(self.name+'.key','w' )
		input_file.write( "FMCSC_SEQFILE "+ self.name + ".seq # NAME OF CAMPARI SPECIFIC SEQUENCE FILE" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_PDBANALYZE 0 # NOT SURE" + "\n" )
		input_file.write( "FMCSC_PDB_FORMAT 1 # SINGLE PDB FILE CONTAINING TRAJECTORY" + "\n" )
		input_file.write( "FMCSC_PDB_READMODE 1 # READ HEAVY ATOMS" + "\n" )
		input_file.write( "FMCSC_PDB_HMODE 1 # READ HYDROGENS, REGENERATE ONES THAT CLASH" + "\n" )
		input_file.write( "FMCSC_PDBFILE " + self.name + ".pdb # NAME OF PDB" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_BOUNDARY 4 # PBC" + "\n" )
		input_file.write( "FMCSC_SHAPE 2 # RECTANGULAR BOX" + "\n" )
		input_file.write( "FMCSC_SIZE" + str(self.pbc) + " # LENGTH OF EDGE OF BOX" + "\n" )
		input_file.write( "FMCSC_RANDOMIZE 1 # MINIMAL RANDOMIZATION WILL BE DONE" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_NRSTEPS 40000000 # TOTAL NUMBER OF STEPS INCLUDING EQUILIBRATION" + "\n" )
		input_file.write( "FMCSC_EQUIL 10000000 # NUMBER OF EQUILIBRATION STEPS" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_SC_ATTLJ 1.0 # ENABLES DISPERSIVE INTERACTIONS" + "\n" )
		input_file.write( "FMCSC_SC_POLAR 1.0 # ENABLES PARTIAL CHARGE INTERACTIONS" + "\n" )
		input_file.write( "FMCSC_SC_IMPSOLV 1.0 # ENABLES DIRECT MEAN FIELD INTERACTIONS AND CHARGE SCREENING" + "\n" )
		input_file.write( "\n" )
		input_file.write( "PARAMETERS " + campari_home + "params/abs3.2_charmm36.prm # ABSINTH FORCE FIELD WITH CHARMM PARTIAL CHARGES" + "\n" )
		input_file.write( "FMCSC_UAMODEL 0 # USE ALL ATOM MODEL, NOT UNITED ATOM" + "\n" )
		input_file.write( "FMCSC_INTERMODEL 1 # EXCLUDE FROZEN-IN INTERACTIONS, i.e. ATOMS IN AROMATIC RINGS" + "\n" )
		input_file.write( "FMCSC_ELECMODEL 2 # USE THE ABSINTH EXCLUSION MODEL FOR SHORT RANGE ELEECTROSTATICS" + "\n" )
		input_file.write( "FMCSC_MODE_14 1 # EXCLUDE 1-4 INTERACTIONS" + "\n" )
		input_file.write( "FMCSC_FUDGE_EL_14 1.0" + "\n" )
		input_file.write( "FMCSC_FUDGE_ST_14 1.0" + "\n" )
		input_file.write( "FMCSC_EPSRULE 2" + "\n" )
		input_file.write( "FMCSC_SIGRULE 1" + "\n" )
		input_file.write( "FMCSC_SC_BONDED_B 0.0" + "\n" )
		input_file.write( "FMCSC_SC_BONDED_A 0.0" + "\n" )
		input_file.write( "FMCSC_SC_BONDED_T 1.0" + "\n" )
		input_file.write( "FMCSC_SC_BONDED_I 0.0" + "\n" )
		input_file.write( "FMCSC_SC_EXTRA 0.0" + "\n" )
		input_file.write( "FMCSC_FOSMID 0.1" + "\n" )
		input_file.write( "FMCSC_FOSTAU 0.25" + "\n" )
		input_file.write( "FMCSC_SCRMID 0.9" + "\n" )
		input_file.write( "FMCSC_SCRTAU 0.5" + "\n" )
		input_file.write( "FMCSC_SAVPROBE 2.5" + "\n" )
		input_file.write( "FMCSC_IMPDIEL 78.2" + "\n" )
		input_file.write( "FMCSC_SCRMODEL 2" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_CUTOFFMODE 4 # topology-assisted cutoffs"+"\n" )
		input_file.write( "FMCSC_MCCUTMODE 2"+"\n" )
		input_file.write( "FMCSC_NBCUTOFF 10.0 # cutoff value in Angstrom for IPP/ATTLJ/DMFI"+"\n" )
		input_file.write( "FMCSC_ELCUTOFF 14.0 # cutoff value in Angstrom for POLAR"+"\n" )
		input_file.write( "FMCSC_LREL_MC 3"+"\n" )
		input_file.write( ""+"\n" )
		input_file.write( "FMCSC_RIGIDFREQ 0.2"+"\n")
		input_file.write( "FMCSC_RIGIDRDFREQ 0.4"+"\n" )
		input_file.write( "FMCSC_RIGIDRDBUF 1.1"+"\n" )
		input_file.write( "FMCSC_CHIFREQ 0.05"+"\n" )
		input_file.write( "FMCSC_CHIRDFREQ 0.4"+"\n" )
		input_file.write( "FMCSC_CHISTEPSZ 20.0"+"\n" )
		input_file.write( "FMCSC_OMEGAFREQ 0.1"+"\n" )
		input_file.write( "FMCSC_OMEGARDFREQ 0.1"+"\n" )
		input_file.write( "FMCSC_OMEGARDFREQ 0.1"+"\n" )
		input_file.write( "FMCSC_PIVOTRDFREQ 0.2"+"\n" )
		input_file.write( "FMCSC_PIVOTSTEPSZ 10.0"+"\n" )
		input_file.write( ""+"\n" )
		input_file.write( "FMCSC_DISABLE_ANALYSIS 1"+"\n" )
		input_file.write( "FMCSC_FLUSHTIME 2.0 # in minutes"+"\n" )
		input_file.write( "FMCSC_SEGCALC 50"+"\n" )
		input_file.write( "FMCSC_BBSEGFILE " + campari_home + "data/bbseg2.dat"+"\n" )
		input_file.write( "FMCSC_DSSPCALC 50"+"\n" )
		input_file.write( "FMCSC_INSTDSSP 0"+"\n" )
		input_file.write( "FMCSC_ANGCALC 50"+"\n" )
		input_file.write( "FMCSC_ENOUT 10000"+"\n" )
		input_file.write( "FMCSC_XYZOUT 2000" + " # FREQUENCY OF SNAPSHOT WRITE OUT " + "\n" )
		input_file.write( "FMCSC_XYZMODE 2" + " # WRITE AS TRAJECTORY (NOT PDBs) " + "\n" )
		input_file.write( "FMCSC_XYZPDB 3"+" # CHARMM STYLE .DCD OUTPUT" + "\n" )
		input_file.write( "FMCSC_RSTOUT 10000" + " # FREQUENCY TO OUTPUT RESTART FILES " + "\n" )
		input_file.write( ""+"\n" )
		input_file.write( "FMCSC_REFILE " + campari_home + "examples/tutorial3/tutorial3.rex"+"\n" )
		input_file.write( "FMCSC_REPLICAS 16"+"\n" )
		input_file.write( "FMCSC_REDIM 1 # the number of exchange dimensions"+"\n" )
		input_file.write( "FMCSC_REMC 1 # enables the replica exchange method in an MPI-parallel simulation run"+"\n" )
		input_file.write( "FMCSC_RESWAPS 15"+"\n" )
		input_file.write( "FMCSC_RENBMODE 2 # only swap with neighboring replicas"+"\n" )
		input_file.write( "FMCSC_REFREQ 2000"+"\n" )
	
	def run_MC(self):
		"""
		Run the CAMPARI subprocess

		"""

		key = self.name + '.key'
		output  = self.name + '.out'
		subprocess.call([campari_home+'bin/Gnu/campari', '-k', key, '&>', output, '&'])
		
		return None

if __name__ == "__main__":
	pbc = 40 # side length of box
	camp = CAMPARI(pdb,pbc)
	check = camp.check_caps()
	if check == 1:
		camp.cap_protein()
	camp.rename_HIS()
	camp.ionize()
	camp.pdb2seq()
	camp.make_input()
	camp.run_MC()
