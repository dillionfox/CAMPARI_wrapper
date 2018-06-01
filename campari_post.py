import numpy as np
import sys
import os
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import campari_analysis

T = campari_analysis.read_hrex('HREX.rex')
#T = [310.0]
nreps = 16
dirlist = sys.argv[1]
FONT = 10

def plot_SS(T,H_av,H_err,E_av,E_err,C_av,C_err):
	plt.clf()
	plt.figure(dpi=600)
	nrows = 5 ; ncols = 4
	for i in range(nrows*ncols):
		ax = plt.subplot(nrows,ncols,i+1)
		ax.plot(H_av[i],label='Helix',c='r') ;      ax.errorbar(range(len(H_av[i])), H_av[i],yerr=H_err[i], fmt='o',color='r',elinewidth=0.5,markersize=0.01)
		ax.plot(E_av[i],label='Beta Sheet',c='b') ; ax.errorbar(range(len(H_av[i])), E_av[i],yerr=E_err[i], fmt='o',color='b',elinewidth=0.5,markersize=0.01)
		ax.plot(C_av[i],label='Coil',c='k') ;       ax.errorbar(range(len(H_av[i])), C_av[i],yerr=C_err[i], fmt='o',color='k',elinewidth=0.5,markersize=0.01)

		if i == 3:
			ax.legend(bbox_to_anchor=(1, 1))

		if i < 16:
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on',labelsize=FONT)
			ax.set_xlim([0,40])
		if i%4 != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on',labelsize=FONT)
		i+=1

	for i in range(nrows*ncols):
		plt.text(-130 + (i%4)*48, 8.7 - (np.floor(i/4)%5) * 1.85 ,'(S'+str(i)+')',fontsize=FONT-3)

	plt.subplots_adjust(top=0.93,bottom=0.1,left=0.08,right=0.78,hspace=0.2,wspace=0.2)
	plt.text(-73, -1.2,'Residue Number',fontsize=FONT)
	plt.text(-160, 4.85,'Frequency',rotation=90,fontsize=FONT)
	plt.savefig('SS_'+str(T)+'.png')
	plt.close()
	#plt.show()

	return None
def plot_cmaps(T,T_cmaps):
	plt.clf()
	plt.figure(dpi=600)
	fig, axes = plt.subplots(nrows=5, ncols=4)
	i = 0
	for ax in axes.flat:
		im = ax.imshow(T_cmaps[i], vmin=0, vmax=1)

		xs = i%4
		ys = np.floor(i/4)%5 * 41.5
		plt.text(-124+xs*45.4, -167+ys,'(S'+str(i)+')',fontsize=FONT-3)

		if i < 16:
			ax.tick_params(axis='x', which='both', bottom='off',top='off',labelbottom='off')
		else:
			ax.tick_params(axis='x', which='both', bottom='on',top='off',labelbottom='on',labelsize=FONT)
		if i%4 != 0:
			ax.tick_params(axis='y', which='both', left='off',right='off',labelleft='off')
		else:
			ax.tick_params(axis='y', which='both', left='on',right='off',labelleft='on',labelsize=FONT)
		i+=1
	plt.text(-75, -175,'T = ' + str(T) + ' K',fontsize=FONT)
	plt.text(-75, 55,'Residue Number',fontsize=FONT)
	plt.text(-160, -90,'Residue Number',rotation=90,fontsize=FONT)
	fig.subplots_adjust(top=0.9,bottom=0.09,left=0.25,right=0.8,hspace=0.15,wspace=0.0)
	cbar_ax = fig.add_axes([0.8, 0.25, 0.01, 0.5])
	cbar = fig.colorbar(im, cax=cbar_ax)
	cbar.ax.tick_params(labelsize=FONT)	

	plt.savefig('CMAPS_'+str(T)+'.png')
	plt.close()

	return None

for replica in range(nreps):
	T_cmaps = []
	H = []; He = []; E = []; Ee = []; C = []; Ce = []
	print T[replica]
	for dir in open(dirlist):
		dir = dir.split('\n')[0]
		name_mod = '_'+str(T[replica]) 
		if os.path.isfile(dir+'/'+"CMAPS"+ name_mod + ".npy") == True:
			cmaps = np.loadtxt(dir+'/'+"CMAPS"+ name_mod + ".npy")
			T_cmaps.append(cmaps)
			[H_av,H_err] = np.loadtxt(dir+'/'+"SS_H" + name_mod + ".npy")
			[E_av,E_err] = np.loadtxt(dir+'/'+"SS_E" + name_mod + ".npy")
			[C_av,C_err] = np.loadtxt(dir+'/'+"SS_C" + name_mod + ".npy")
			H.append(H_av)
			He.append(H_err)
			E.append(E_av)
			Ee.append(E_err)
			C.append(C_av)
			Ce.append(C_err)
	plot_cmaps(T[replica], T_cmaps)
	plot_SS(T[replica],H,He,E,Ee,C,Ce)
