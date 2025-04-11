from input import *

myhf = scf.RHF(mol).run()
#myhf.kernel()

fout = open('mo_energies.txt','w')
orb = myhf.mo_coeff
orb_energy = myhf.mo_energy
nmo = len(orb_energy)
for i, e in enumerate(orb_energy):
 print(i,e,file=fout)
fout.close()

# localized orbitals

#orb = lo.cholesky.cholesky_mos(myhf.mo_coeff)

#mo = Boys(mol, myhf.mo_coeff[:,myhf.mo_occ>0])
#mo.kernel()
#loc_orb = mo.mo_coeff

#mo = Boys(mol, myhf.mo_coeff[:,myhf.mo_occ==0])
#mo.kernel()
#vloc_orb = mo.mo_coeff

#orb = np.concatenate((loc_orb, vloc_orb), axis=1)

#mo = Boys(mol, myhf.mo_coeff[:,0:27])
#mo.kernel()
#loc_orb = mo.mo_coeff

#mo = Boys(mol, myhf.mo_coeff[:,27:nmo-1])
#mo.kernel()
#vloc_orb = mo.mo_coeff

#torb = np.concatenate((loc_orb, vloc_orb), axis=1)
#orb = np.concatenate((torb, myhf.mo_coeff[:,nmo-1:nmo]), axis=1)

#mo = lo.PM(mol,myhf.mo_coeff)
#mo.kernel()
#orb = mo.mo_coeff

#orb = lo.orth_ao(myhf, 'nao')

# lets start with the Wigner-Weisskopf approach # following notation from "On the computations of decay widths of Fano resonances"
argon = [ i  for i in  range(12,32)] # 3p orbitals of Ar

# saves the two-electron integrals in the file ftmp.name
ao2mo.kernel(mol, orb, erifile = 'hf.h5', dataname = 'test')
# load 2e integrals by filename and dataname
with ao2mo.load('hf.h5', 'test') as eri:
    eri1 = ao2mo.restore(1, np.asarray(eri), orb.shape[1])

# Define the path to the compressed .npz file
#compressed_file_path = 'data_compressed.npz'

# Save the large NumPy array to a compressed .npz file
#np.savez_compressed(compressed_file_path, data_array=eri1)

# Load the NumPy array from the compressed .npz file
#loaded_data = np.load(compressed_file_path)
#eri = loaded_data['data_array']

for k in argon:
   print()
   for l in argon:
       #if(k>=l):
          print(k,l,eri1[k,k,l,l],eri1[k,l,k,l])

