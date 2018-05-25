import sys
import subprocess
import os
import numpy as np
import mdtraj as md
import structure_analysis as sa

"""
Requires "structure_analysis.py", which can be found in the reflectin/structure_analysis folder

"""


def read_hrex(fil):
	"""
	Read in temperatures used in replica exchange

	"""
	t = []
	for line in open(fil):
		l = line.split('\n')[0]
		t.append(float(l))
	return t[1:]

def rextrace(fil):
	"""
	CAMPARI prints out trajectories for each replica, not for each temperature.
	We use the conversion table they give us to construct trajectories with 
	constant temperature. Here we take the inverse of the conversion table.

	"""
	skip = 1000000
	dcdfreq = 1000
	hrexfreq = 2000
	c = 0 
	# pre-allocate memory
	for line in open(fil): 
		if int(line.split()[0]) < skip:
			continue
		if c == 0:
			l = len(line.split())-1
		c+=1
	repmat = np.zeros((c,l))

	# populate relica matrix
	t = np.zeros(c)
	c = 0
	for line in open(fil):
		if int(line.split()[0]) < skip:
			continue
		line = line.split()
		t[c] = (float(line[0]) - skip)/ dcdfreq
		reps = np.array(line[1:])
		reps = [int(i) for i in reps]
		repmat[c] = reps
		c+=1
	return [repmat, t]

def reconstruct(repmat,t,T,dir):
	"""
	Use the inverted conversion table to construct new trajectories with
	constant temperature. Center each frame, too.

	"""
	trajs = []
	psf = os.getcwd() + '/' + dir + '/' + dir + '_autogen.psf'
	pdb = os.getcwd() + '/' + dir + '/' + 'N_000_yourcalc_START.pdb'
	# load each trajectory and append it to a tuple. Also make a pdb with same exact coordinates.
	for i in range(nreps):
		if i < 10:
			ii = '0' + str(i)
		else:
			ii = str(i)
		dcd = os.getcwd() + '/' + dir + '/' + 'N_0' + ii + '_yourcalc_traj.dcd'
		traj = md.load_dcd(dcd,top=pdb)
		trajs.append(traj)
		peptide = traj.atom_slice(traj.topology.select('(not resname ACE) and (not resname NME)'))[0]
		peptide.save(dir + '/T.pdb')

	# For each replica, go through each step in the inverted conversion table and write the new trajectory
	for replica in range(nreps):
		print "Writing new trajectory. This may take a while... Replica:", replica
		with md.formats.DCDTrajectoryFile(dir + '/T_'+str(T[replica]) + '.dcd', 'w') as f: 
			for interval in range(0,t.size*2-2,2):
				n = int(repmat[:,replica][interval/2])-1
				structure = trajs[n][interval] 
				peptide = structure.atom_slice(structure.topology.select('(not resname ACE) and (not resname NME)'))[0]
				peptide.center_coordinates()
				f.write(peptide.xyz*10)
	return None

def analyze(dir,T):
	pdb = dir + '/T.pdb'
	for replica in range(nreps):
		print "working on T =", T[replica], "..."
		dcd = dir + '/T_'+str(T[replica]) + '.dcd'
		analysis = sa.SA(dcd, pdb, '_'+str(T[replica]))
		analysis.run(dir)	

nreps = 16
dirlist = sys.argv[1]
REXfile = 'N_000_REXTRACE.dat'
pdbname = 'N_000_yourcalc_START.pdb'

for dir in open(dirlist):
	dir = dir.split('\n')[0] 
	make_dcds = False
	T = read_hrex(os.getcwd() + '/' + 'HREX.rex')
	for replica in range(nreps):
		if os.path.isfile(dir + '/T_'+str(T[replica]) + '.dcd') != True:
			make_dcds = True
	if make_dcds == True:
		[repmat,t] = rextrace(dir+'/'+REXfile)		# Read replica exchange file
		reconstruct(repmat,t,T,dir)
	analyze(dir,T)
