
##################################
###      CAMPARI TOOLS         ### 
###                            ###
###       Dillion Fox          ###
###          5/2018            ###
### University of Pennsylvania ###
###   Oak Ridge National Lab   ###
##################################

import sys
import os
import MDAnalysis as mda
import MDAnalysis.analysis.distances as mdad
import subprocess
import numpy as np

if sys.argv[1] == '-h' or sys.argv[1] == '--h' or sys.argv[1] == '--help':
	print
	print
	print "                                ~~~ CAMPARI TOOLS ~~~" 
	print
	print
	print "This code makes it easy to generate input files to run CAMPARI (http://campari.sourceforge.net/),"
	print "a Monte Carlo package developed by the Pappu lab. This code is only set up generate input files"
	print "for Replica Exchange Monte Carlo simulations with a pre-compiled, MPI version of CAMPARI"
	print "This code can either 1. Take a CAMPARI-style .seq file and generate and run a REMC simulation,"
	print "or it can take a pdb of a protein, add ACE and NME residues to the ends, determine the protonation"
	print "states of histidine and rename the residues accordingly, and neutralize the system by adding ions"
	print 
	h = 'y'
h = 'n'

if len(sys.argv) != 2 or h == 'y':
	print "USEAGE: python campari_tools.py input"
	print "USEAGE: input is EITHER a PDB or a SEQ file"
	print
	print "example:"
	print
	print "python campari_tools.py 1ova.pdb"
	print "or"
	print "python campari_tools.py 1ova.seq"
	exit()

# Run CAMPARI, or just generate inputs?
run = 'n'

# The extension of the input argument determines the functions that need to be run
NAME = sys.argv[1].split('.')[0]
MODE = sys.argv[1].split('.')[1]

# Provide the path to the force field. The code will also check the current directory for the force field.
campari_home = '/home/dillion/pkg/campari/'
ff_path = campari_home+'params/abs3.2_charmm36.prm'
bbseg_path = campari_home+'data/bbseg2.dat'

#
# This class contains all of the main functions. It can generate the ACE and NME pdb's if they don't exist, along
# with all of the input files EXCEPT the force field.
#
class CAMPARI:

	def __init__(self, pbc, pdb='NULL'):
		self.pdb = pdb
		self.pbc = pbc
		self.name = NAME

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
		This function requires "ace.pdb" and "nme.pdb". I can provide them upon request.

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
		HIS residues with the appropriate names. This function relies on ROSETTA's naming scheme.

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
			else:
				his.residues.resnames = 'HIE'

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
		for aa in negative_charge:
			minus+=len(np.where(seq==aa)[0])
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
		input_file.write( "FMCSC_SEQFILE "+ self.name + ".seq # name of campari specific sequence file" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_PDBANALYZE 0 # PDB- keywords are for inputs. This one is a bit ambiguous" + "\n" )
		input_file.write( "FMCSC_PDB_FORMAT 1 # input format: 1. single pdb containing trajectory" + "\n" )
		input_file.write( "FMCSC_PDB_READMODE 1 # read heavy atoms" + "\n" )
		input_file.write( "FMCSC_PDB_HMODE 1 # read hydrogens, regenerate ones that clash" + "\n" )
		if MODE == 'pdb':
			input_file.write( "FMCSC_PDBFILE " + self.name + ".pdb # input file (pdb)" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_BOUNDARY 4 # 1. pbc, 2. hard-wall boundary, 3. residue-based soft wall, 4. atom-based soft wall" + "\n" )
		input_file.write( "FMCSC_SHAPE 2 # 1. rectangular box, 2. sphere, 3. cylinder" + "\n" )
		input_file.write( "FMCSC_SIZE" + str(self.pbc) + " # depends on 'SHAPE'. 1. 3-vector, 2. scalar, 3. two floats" + "\n" )
		input_file.write( "FMCSC_RANDOMIZE 1 # minimal randomization will be done" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_NRSTEPS 10000000 # total number of steps including equilibration" + "\n" )
		input_file.write( "FMCSC_EQUIL 1000000 # number of equilibration steps" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_SC_ATTLJ 1.0 # enables dispersive interactions" + "\n" )
		input_file.write( "FMCSC_SC_WCA 0.0 # scales Weeks-Chandler-Andersen potential" + "\n" )
		input_file.write( "FMCSC_SC_POLAR 1.0 # enables partial charge interactions" + "\n" )
		input_file.write( "FMCSC_SC_IMPSOLV 1.0 # enables direct mean field interactions and charge screening" + "\n" )
		input_file.write( "\n" )
		input_file.write( "PARAMETERS " + ff_path + " # ABSINTH force field with CHARMM partial charges" + "\n" )
		input_file.write( "FMCSC_UAMODEL 0 # use all-atom model, not united atom" + "\n" )
		input_file.write( "FMCSC_INTERMODEL 1 # exclude frozen-in interactions, i.e. atoms in aromatic rings" + "\n" )
		input_file.write( "FMCSC_ELECMODEL 2 # use the ABSINTH exclusion model for short range electrostatics" + "\n" )
		input_file.write( "FMCSC_MODE_14 1 # exclude 1-4 interactions" + "\n" )
		input_file.write( "FMCSC_FUDGE_EL_14 1.0 # scale factor for elec interactions" + "\n" )
		input_file.write( "FMCSC_FUDGE_ST_14 1.0 # scale factor for steric interactions" + "\n" )
		input_file.write( "FMCSC_EPSRULE 2 # NOT IN DOCUMENTATION!!!" + "\n" )
		input_file.write( "FMCSC_SIGRULE 1 # NOT IN DOCUMENTATION!!!" + "\n" )
		input_file.write( "FMCSC_SC_BONDED_B 0.0 # linear scaling for all bond lengths" + "\n" )
		input_file.write( "FMCSC_SC_BONDED_A 1.0 # linear scaling for all bond angles" + "\n" )
		input_file.write( "FMCSC_SC_BONDED_T 1.0 # linear scaling for all torsional potentials" + "\n" )
		input_file.write( "FMCSC_SC_BONDED_I 1.0 # linear scaling for all improper dihedrals" + "\n" )
		input_file.write( "FMCSC_SC_EXTRA 0.0 # scaling factor applied to rotatable bonds with electronic effects not captured by atomic pairwise interactions" + "\n" )
		input_file.write( "FMCSC_SC_IPP 1.0 # scales inverse power potential" + "\n" )
		input_file.write( "FMCSC_FOSTAU 0.25 # parameter determining the steepness of the sigmoidal interpolation for mapping SA volumes to SA states" + "\n" )
		input_file.write( "FMCSC_FOSMID 0.1 # related to FOSTAU, but this is the 'shift' applied to the sigmoid" + "\n" )
		input_file.write( "FMCSC_SCRMID 0.9 # just like FOSMID, but for charge screening" + "\n" )
		input_file.write( "FMCSC_SCRTAU 0.5 # just like FOSTAU, but for charge screening" + "\n" )
		input_file.write( "FMCSC_SAVPROBE 2.5 # IMPORTANT. size of solvation shell around atoms (Ang)." + "\n" )
		input_file.write( "FMCSC_IMPDIEL 78.2 # 'implicit dielectric': in ABSINTH this is the parameter that controls electrostatic screening" + "\n" )
		input_file.write( "FMCSC_SCRMODEL 2 # IMPORTANT. This determines how dielectric screening is done. Option 2 goes with IMPDIEL and localizes and strengthens specific interactions" + "\n" )
		input_file.write( "\n" )
		input_file.write( "FMCSC_CUTOFFMODE 4 # topology-assisted cutoffs"+"\n" )
		input_file.write( "FMCSC_MCCUTMODE 2 # monte carlo cutoff approach. Option 2 applies a residue level exclusion approach similar to MD"+"\n" )
		input_file.write( "FMCSC_NBCUTOFF 10.0 # cutoff value in Angstrom for IPP/ATTLJ/DMFI"+"\n" )
		input_file.write( "FMCSC_ELCUTOFF 14.0 # cutoff value in Angstrom for POLAR"+"\n" )
		input_file.write( "FMCSC_LREL_MC 3 # monopole-monopole int. computed explicitly. polyatomic monopole groups are represented by total charge at center of geo. dipole int skipped"+"\n" )
		input_file.write( ""+"\n" )
		input_file.write( "FMCSC_RIGIDFREQ 0.2 # very poor description in docs. I think its the fraction of MC moves that are 'rigid'."+"\n")
		input_file.write( "FMCSC_RIGIDRDFREQ 0.1 # randomizes a fraction of the underlying tree of possible MC moves"+"\n" )
		input_file.write( "FMCSC_RIGIDRDBUF 1.1 # avoids a systematic bias from too many interactions with boundary potentials? it scales distances to reject moves that get close to wall"+"\n" )
		input_file.write( "FMCSC_CHIFREQ 0.3 # fraction of all sidechain moves including a specialized move type used for analysis only"+"\n" )
		input_file.write( "FMCSC_CHIRDFREQ 0.4 # docs not clear. I think it randomizes the tree to promote the acceptance of bigger moves"+"\n" )
		input_file.write( "FMCSC_CHISTEPSZ 20.0 # step size for chi twists"+"\n" )
		input_file.write( "FMCSC_CRFREQ 0.3 # frequency of concerted rotation moves"+"\n" )
		input_file.write( "FMCSC_OMEGAFREQ 0.3 # Omega only ever takes 2 values for proline, and one for others. But small sampling is important"+"\n" )
		input_file.write( "FMCSC_OMEGARDFREQ 0.1 # also not explained clearly, but something to do with randomizing probability tree"+"\n" )
		input_file.write( "FMCSC_PIVOTRDFREQ 0.3 # docs not clear. I think it randomizes the tree to promote the acceptance of bigger moves (phi/psi)"+"\n" )
		input_file.write( "FMCSC_PIVOTSTEPSZ 10.0 # maximum step size for perturbing phi/psi"+"\n" )
		input_file.write( "FMCSC_TRANSSTEPSZ 10.0 # maximum step size for perturbing translations"+"\n" )
		input_file.write( "FMCSC_ROTSTEPSZ 20.0 # maximum step size for perturbing rotations"+"\n" )
		input_file.write( "FMCSC_ROTFREQ 0.1 # frequency of purely rotational moves, requiring COUPLRERIGID be false"+"\n" )
		input_file.write( "FMCSC_PKRFREQ 0.2 # frequency of proline pucker moves"+"\n" )
		input_file.write( "FMCSC_PKRRDFREQ 0.02 # frequency of reflection moves in proline pucke"+"\n" )
		input_file.write( "FMCSC_PUCKERSTEP_DI 4.0 # max step size for torsions in pucker moves"+"\n" )
		input_file.write( "FMCSC_PUCKERSTEP_AN 2.0 # max step size for angles in pucker moves"+"\n" )
		input_file.write( "FMCSC_ALIGN 4 # how the molecule swivels relative to the C and N termini"+"\n" )
		input_file.write( ""+"\n" )
		input_file.write( "FMCSC_DISABLE_ANALYSIS 2 # 1. all analysis and inst. output disabled. 2. all analysis disabled. 3. all inst. output disabled."+"\n" )
		input_file.write( "FMCSC_FLUSHTIME 2.0 # how often performance estimates are given (in minutes)"+"\n" )
		input_file.write( "FMCSC_SEGCALC 50 # how often secondary structure is determined based on phi/psi"+"\n" )
		input_file.write( "FMCSC_BBSEGFILE " + bbseg_path + " # Ramachandran plot. used for SS determination"+"\n" )
		input_file.write( "FMCSC_DSSPCALC 50 # how often secondary structure is estimated based on h-bond patterns"+"\n" )
		input_file.write( "FMCSC_DSSP 0.0 # scaling factor for DSSP aligning potential"+"\n" )
		input_file.write( "FMCSC_INSTDSSP 0 # NOT IN DOCUMENTATION!!!"+"\n" )
		input_file.write( "FMCSC_ZSEC 0.0 # scaling factor for global secondary structure bias"+"\n" )
		input_file.write( "FMCSC_TOR 0.0 # scaling factor for controlling external torsional bias terms"+"\n" )
		input_file.write( "FMCSC_ANGCALC 50 # NOT IN DOCUMENTATION!!!"+"\n" )
		input_file.write( "FMCSC_DREST 0.0 # Scaling factor for externally defined harmonic distance restraints"+"\n" )
		input_file.write( "FMCSC_TABUL 0.0 # Scaling factor for externally defined tabulated potentials"+"\n" )
		input_file.write( "FMCSC_POLY 0.0 # Scaling factor for restraint potentials on polymeric properties"+"\n" )
		input_file.write( "FMCSC_ENOUT 1000 # NOT IN DOCUMENTATION!!! But it has to be how often energies are written to file"+"\n" )
		input_file.write( "FMCSC_XYZOUT 1000" + " # how often to write coordinates to trajectory " + "\n" )
		input_file.write( "FMCSC_XYZMODE 2" + " # (2) write as trajctory (not individual PDBs) " + "\n" )
		input_file.write( "FMCSC_XYZPDB 3"+" # (3) is CHARMM style .dcd output" + "\n" )
		input_file.write( "FMCSC_RSTOUT 10000" + " # how often to write restart files " + "\n" )
		input_file.write( ""+"\n" )
		input_file.write( "FMCSC_REFILE " + self.name + ".rex" + "# this file defines the replica exchange method"+"\n" )
		input_file.write( "FMCSC_REPLICAS 16 # desired number of replicas. NOTE: this can be overridden if MPI is used"+"\n" )
		input_file.write( "FMCSC_REDIM 1 # the number of exchange dimensions"+"\n" )
		input_file.write( "FMCSC_REMC 1 # enables the replica exchange method in an MPI-parallel simulation run"+"\n" )
		input_file.write( "FMCSC_RESWAPS 15 # number of replicas to randomly swap every REFREQ steps"+"\n" )
		input_file.write( "FMCSC_RENBMODE 2 # only swap with neighboring replicas (similar to hamiltonian RE)"+"\n" )
		input_file.write( "FMCSC_REFREQ 2000 # how often to swap replicas"+"\n" )
		input_file.close()

	def make_HREX(self, N, Ti = 250, Tf = 400):
		"""
		This function writes the temperatures to a file for Hamiltonian Replica Exchange

		"""
		HREX_file = open(self.name+'.rex','w' )
		HREX_file.write( "1"+"\n" )
		for i in np.linspace(Ti, Tf, N):
			HREX_file.write( str(i) + "\n" )
		return None
	
	def make_nme(self):
		"""
		If you don't already have a pdb for this, it will make one

		"""
		nme_pdb = open("nme.pdb", 'w')
		nme_pdb.write("ATOM      1  N   ALA     1      25.072  26.282   1.453  1.00  0.00           N \n") 
		nme_pdb.write("ATOM      2  CA  ALA     1      25.768  27.185   2.366  1.00  0.00           C \n")
		nme_pdb.write("ATOM      3  C   ALA     1      27.248  26.890   2.395  1.00  0.00           C \n")
		nme_pdb.write("ATOM      4  O   ALA     1      27.749  25.990   1.714  1.00  0.00           O \n")
		nme_pdb.write("ATOM      5  N   NME     1A     28.047  27.661   3.214  1.00  0.00           N \n")
		nme_pdb.write("ATOM      6  CA  NME     1A     29.480  27.414   3.276  1.00  0.00           C \n")
		nme_pdb.write("ATOM      7  H   NME     1A     27.637  28.397   3.771  1.00  0.00           H \n")
		nme_pdb.write("ATOM      8 1HA  NME     1A     29.728  26.570   2.594  1.00  0.00           H \n")
		nme_pdb.write("ATOM      9 2HA  NME     1A     30.033  28.327   2.960  1.00  0.00           H \n")
		nme_pdb.write("ATOM     10 3HA  NME     1A     29.773  27.155   4.318  1.00  0.00           H \n")
		nme_pdb.close()
		return None

	def make_ace(self):
		"""
		If you don't already have a pdb for this, it will make one

		"""
		ace_pdb = open("ace.pdb", 'w')
		ace_pdb.write("ATOM      1  CH3 ACE A   1      14.909  41.084  28.582  1.00  0.00           C \n") 
		ace_pdb.write("ATOM      2 HH32 ACE A   1      14.736  41.787  29.399  1.00  0.00           H \n")
		ace_pdb.write("ATOM      3 HH33 ACE A   1      13.978  40.912  28.038  1.00  0.00           H \n")
		ace_pdb.write("ATOM      4 HH31 ACE A   1      15.253  40.130  28.989  1.00  0.00           H \n")
		ace_pdb.write("ATOM      5  C   ACE A   1      15.962  41.641  27.676  1.00  0.00           C \n")
		ace_pdb.write("ATOM      6  O   ACE A   1      17.150  41.342  27.794  1.00  0.00           O \n")
		ace_pdb.write("ATOM      7  N   ALA A   2      15.437  42.533  26.748  1.00  0.00           N \n")
		ace_pdb.write("ATOM      8  CA  ALA A   2      16.343  43.170  25.795  1.00  0.00           C \n")
		ace_pdb.write("ATOM      9  CB  ALA A   2      15.736  43.188  24.393  1.00  0.00           C \n")
		ace_pdb.write("ATOM     10  C   ALA A   2      16.717  44.562  26.282  1.00  0.00           C \n")
                ace_pdb.write("ATOM     11  O   ALA A   2      15.867  45.302  26.799  1.00  0.00           O \n")
		ace_pdb.close()
		return None

	def locate_ff(self):
		"""
		This code tries to generate all missing files, but the force field is too long to do this. This
		function tries to find the force field if it's not in the directory

		"""
		print "Force field not found. We'll try to find it for you..."
		print "locate abs3.2_charmm36.prm"
		subprocess.call(['locate', 'abs3.2_charmm36.prm'])
		print "if the force field was not found, download it from: http://campari.sourceforge.net/"
		exit()

	def make_bbseg(self):
		"""
		If you don't already have this file, it will make it

		"""

		bbseg_file = open("bbseg2.dat", 'w')
		bbseg_file.write(" 1 1 1 1 1 1 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n") 
		bbseg_file.write(" 1 1 1 1 1 1 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 \n")
		bbseg_file.write(" 1 1 1 1 1 1 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 1 1 1 1 1 1 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 1 1 1 1 1 1 2 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 \n")
		bbseg_file.write(" 1 1 1 1 1 1 2 2 2 2 2 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 1 1 1 1 1 0 2 2 2 5 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 1 1 1 0 0 0 0 5 5 5 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 5 5 5 5 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 3 3 3 0 0 5 5 5 5 5 5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 3 3 3 3 3 0 5 5 5 5 5 0 0 0 0 0 0 0 0 0 0 8 8 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 3 3 3 3 3 0 5 5 5 5 0 0 0 0 0 0 0 0 0 0 0 8 8 8 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 3 3 3 3 3 0 5 5 0 0 0 0 0 0 0 0 0 0 0 8 8 8 8 8 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 3 3 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 8 8 8 8 8 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 3 3 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 8 8 8 8 8 8 8 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 3 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8 8 8 8 8 8 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 3 3 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 8 8 8 8 8 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 8 8 8 8 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 8 8 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 7 7 7 0 0 0 0 0 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 \n")
		bbseg_file.write(" 7 7 7 0 0 0 0 0 0 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 \n")
		bbseg_file.write(" 7 7 7 7 0 0 0 0 0 0 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7 \n")
		bbseg_file.write(" 7 7 7 7 0 0 0 0 0 0 0 4 4 4 0 0 0 0 0 0 0 0 0 0 0 6 6 0 0 0 0 0 0 0 0 7 \n")
		bbseg_file.write(" 7 7 7 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 6 0 0 0 0 0 0 0 0 7 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 6 6 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 6 6 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 6 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 6 6 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 1 1 1 1 1 1 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.write(" 1 1 1 1 1 1 2 2 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \n")
		bbseg_file.close()

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

	if MODE == 'seq':
		# In this case, CAMPARI will generate a random structure
		camp = CAMPARI(pbc)
	elif MODE == 'pdb':
		# If we give CAMPARI a structure, it has to meet specific criteria
		pdb = sys.argv[1]
		camp = CAMPARI(pbc,pdb)
		# all proteins must be 'capped' with ACE and NME. If you don't have pdbs for these, write them
		check = camp.check_caps()
		if check == 1:
			if os.path.isfile("nme.pdb") != True:
				camp.make_nme()
			if os.path.isfile("ace.pdb") != True:
				camp.make_ace()
			camp.cap_protein()
		# Check to see if you already have bbseg2.dat file
		if os.path.isfile(bbseg_path) != True:
			bbseg_path = os.getcwd() + '/bbseg2.dat'
			if os.path.isfile(bbseg_path) != True:
				# if it's not in CAMPARI_HOME or in the current directory, just make a new one
				camp.make_bbseg()
		# All Histidine residues must be named according to protonation states. This function takes care of that.
		camp.rename_HIS()
		# Add ions to neutralize sytem
		camp.ionize()
		# Make the CAMPARI specific '.seq' file
		camp.pdb2seq()
	# the force field is the only non-input requirement that the code won't generate. Make sure it exists.
	if os.path.isfile(ff_path) != True:
		ff_path = os.getcwd() + '/abs3.2_charmm36.prm'
		if os.path.isfile(ff_path) != True:
			print "force field is not in expected location. We will try to find it for you"
			camp.locate_ff()

	# Write the Hamiltonian Replica Exchange (HREX) temperatures into a file
	# Note: HREX will only execute if the number of temperatures is equal to the number of cores. Serial executions do not run HREX even if specified.
	camp.make_HREX(16)
	# Write a CAMPARI input file and include all of the above information
	camp.make_input()
	# Run CAMPARI
	if run == 'y' or run == 'yes' or run == 'on':
		camp.run_MC()
