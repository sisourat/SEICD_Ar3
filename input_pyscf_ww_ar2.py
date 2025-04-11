from input import *
#print(mol.ao_labels())

myhf = scf.RHF(mol)
myhf.kernel()

orb = myhf.mo_coeff
orb_energy = myhf.mo_energy
nmo = len(orb_energy)

# lets start with the Wigner-Weisskopf approach # following notation from "On the computations of decay widths of Fano resonances"
a4s1  = 18 # 4s orbital of fake-Argon
e4s1 = orb_energy[a4s1]

a3p1 = [12,13,14] # 3p orbitals of fake-Argon
#orb_energy[a3p1] = -0.6575151921924877

#a4s2  = 20 # 4s orbital of Argon
#e4s2 = orb_energy[a4s2]

a3p2 = [15,16,17] # 3p orbitals of Argon

virt = [ i  for i in  range(22,125)] # virtual (continuum) orbitals

for i, e in enumerate(orb_energy):
 print(i,e)
print()
# load 2e integrals by filename and dataname
with ao2mo.load('hf.h5', 'test') as eri:
    eri1 = ao2mo.restore(1, np.asarray(eri), orb.shape[1])

i1 =  a3p1[0]
ei1 = orb_energy[i1]
i2 =  a3p1[1]
ei2 = orb_energy[i2]

b = a4s1
e0 = e4s1 -ei1 -ei2 + eri1[i1,i1,i2,i2] + eri1[i1,i2,i1,i2] - eri1[i1,i1,b,b] + 0.5*eri1[i1,b,b,i1] - eri1[i2,i2,b,b] + 0.5*eri1[i2,b,b,i2]

fout = open('fanotot0.txt','w')
for k in  a3p1:
    for l in  a3p2:
        for a in  virt:
            ea = orb_energy[a]
            ek = orb_energy[k]
            el = orb_energy[l]
            # singlet
            ef = ea -ek -el + eri1[k,k,l,l] + eri1[k,l,l,k] - eri1[k,k,a,a] + 0.5*eri1[k,a,a,k] - eri1[l,l,a,a] + 0.5*eri1[l,a,a,l]
            coup = 0.0
            if(k==i1):
             coup += eri1[b,i2,a,l]
             coup -= 0.5*eri1[b,l,a,i2]
            if(k==i2):
             coup += eri1[b,i1,a,l]
             coup -= 0.5*eri1[b,l,a,i1]
            print(ef-e0,coup**2,file=fout)
            # triplet
            ef = ea -ek -el + eri1[k,k,l,l] - eri1[k,l,l,k] - eri1[k,k,a,a] + 1.5*eri1[k,a,a,k] - eri1[l,l,a,a] + 1.5*eri1[l,a,a,l]
            if(k==i1):
             coup += eri1[b,i2,a,l]
             coup += eri1[b,l,a,i2]
            if(k==i2):
             coup += eri1[b,i1,a,l]
             coup += eri1[b,l,a,i1]
            print(ef-e0,coup**2/3.0,file=fout)

sys.exit()

fout = open('fanotot0.txt','w')

for a in virt:
    for k in argon:
        for l in argon:
# S1
            if(k==l):
             ef = (orb_energy[a]-orb_energy[k]-orb_energy[l]) + eri1[k,k,l,l]  - eri1[k,k,a,a] + 0.5*eri1[k,a,a,k] - eri1[l,l,a,a] + 0.5*eri1[l,a,a,l]
            else:
             ef = (orb_energy[a]-orb_energy[k]-orb_energy[l]) + eri1[k,k,l,l] + eri1[k,l,l,k]  - eri1[k,k,a,a] + 0.5*eri1[k,a,a,k] - eri1[l,l,a,a] + 0.5*eri1[l,a,a,l]
            if(eri1[k,k,l,l]<0.1):
             coup = np.sqrt(2.0)*(eri1[k,i,l,i]-eri1[a,a,k,l])+eri1[a,l,k,a]/np.sqrt(2.0)
             print(ef-e0,0.5*coup**2+1e-10,a,k,l,eri1[k,k,l,l],eri1[k,l,l,k], file=fout)
# T
            #ef = (orb_energy[a]-orb_energy[k]-orb_energy[l]) + eri1[k,k,l,l] - eri1[k,l,l,k] #- eri1[k,a,k,a] + 0.5*eri1[k,a,a,k] - eri1[l,a,l,a] + 0.5*eri1[l,a,a,l]
            #coup = (eri1[a,l,i,k]-eri1[a,k,i,l])
            #print(ef-ei,1.5*coup**2+1e-10)


fout.close()
