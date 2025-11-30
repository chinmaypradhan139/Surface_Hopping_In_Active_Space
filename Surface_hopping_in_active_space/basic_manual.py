import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

au2J=4.35974417e-18
kb=1.3806503e-23
mol=1
benchmark=0
Ne=10
Hi=20
band_width=1.0e-2
Gamma=0.005
spacing=2



KT=0.001
#------------------
# input to fortran

mass=2000
lamb=0.008
omg=0.0002
exo=-0.003


inarray=[str(mass),str(lamb),str(omg),str(exo)]
#----------------------

coup = np.sqrt(Gamma/(2*np.pi))
E_corr=0.00
mol=list(range(1,mol+1))
print(mol,'mol')
metal=list(range(len(mol)+1,Hi+1))
print(metal,'metal',int(len(metal)/2))
orbitals=mol+list(metal)
print(orbitals,'orbitals')
filled_met=list(range(len(mol)+1,Ne+len(mol)))
print(filled_met,'filled met')
print(mol[-1],'type mol')
#gs=[]
#gs.append(mol[-1])
#gs=gs+filled_met
#gs=filled_met
gs=list(range(2,Ne+len(mol)))
print(gs,'gs')
vacant_met=list(set(metal).symmetric_difference(set(filled_met)))
print(vacant_met,'vacant')

equal_energies=list(np.linspace(-band_width/2,band_width/2,len(metal)))


def gaussian_quadrature_mapping_with_weights(n: int, bw: float = band_width):
    # Step 1: Generate the list of integers from 3 to n
    integer_list = list(range(len(mol)+1, n + 1))
    num_integers = len(integer_list)

    # Step 2: Calculate Gaussian quadrature nodes and weights for the given number of integers
    nodes, weights = np.polynomial.legendre.leggauss(int(num_integers/2))
 

#    exit()
    # Step 3: Scale the nodes from [-1, 1] to [-bw/2, bw/2]
    scaled_nodes = [-(0.5+0.5*i) * (bw / 2) for i in nodes]
    scaled_nodes.reverse()
    scaled_nodes1=[(0.5+0.5*i) * (bw / 2) for i in nodes]

    scaled_nodes=scaled_nodes+scaled_nodes1
    weights=[np.sqrt(coup/(2*np.pi))*np.sqrt(band_width*i)/2.0 for i in weights]+[np.sqrt(coup/(2*np.pi))*np.sqrt(band_width*i)/2.0 for i in weights]

#    weights=weights+weights
    # Step 4: Map integers to the corresponding scaled nodes and weights
    mapping = {integer: (node, weight) for integer, node, weight in zip(integer_list, scaled_nodes, weights)}
    
    return mapping

#print(int(len(metal)/2),'new_metal')
gauss_orbs = gaussian_quadrature_mapping_with_weights(len(metal)+2, band_width)
if (spacing==1):
    met_map={i:gauss_orbs[i][0] for i in range(len(mol)+1,Hi+1)}
if (spacing==2):
    met_map={metal[i]:equal_energies[i] for i in range(len(metal))}

#print(met_map,'met_map')
#plt.plot([i*0 for i in range(3,Hi)],[met_map[i] for i in range(3,Hi)],'o')
#plt.show()
frac=3.0
if (benchmark==0):
    active=[]
    for orb in filled_met:
        diff=met_map[vacant_met[0]]-met_map[orb]
        if (abs(diff)<frac*KT):
            active.append(orb)
    for orb in vacant_met:
        diff=met_map[filled_met[len(filled_met)-1]]-met_map[orb]
        if (abs(diff)<frac*KT):
            active.append(orb)
else:
    active=metal

active=mol+active


common=list(set(gs) & set(active))
print(common,'common')
active_electrons=len(common)

print(active,active_electrons,'active','active_electrons')
print(common)
comb = list(combinations(active, active_electrons))

frozen=list(set(orbitals).symmetric_difference(set(active)))
print(frozen,'frozen')
for i in range(len(comb)):
    print(i+1,comb[i])

nquant=len(comb)

print(nquant,'nquant')

filled_frozen=list(set(frozen) & set(filled_met))


frozen_energy=sum([met_map[i] for i in filled_frozen])

Hamil_site=np.zeros((len(comb),len(comb)))

for i in range(len(comb)):
    Hamil_site[i][i]=frozen_energy+sum([met_map[k] for k in comb[i] if k>len(mol)])
    if (len(mol)==2):
        if (1 in comb[i] and 2 in comb[i]):
            Hamil_site[i][i]=E_corr+frozen_energy


rho=len(metal)/band_width
equi_coup=np.sqrt(Gamma/(2*np.pi*rho))

print(equi_coup,'coup')

print(len(mol),'lenmol')

if (len(mol)==1):
    for com in comb:
        for bom in comb:
            diff=list(set(com).symmetric_difference(set(bom)))
            #if (len(diff)==len(mol)):

            if (1 in diff):
                    cp=[l for l in diff if l!=1]
                    if (spacing==1):
                        Hamil_site[comb.index(com)][comb.index(bom)]=gauss_orbs[cp[0]][1]
                    elif (spacing==2):
                        Hamil_site[comb.index(com)][comb.index(bom)]=equi_coup




if (len(mol)==2):
    for com in comb:
        for bom in comb:
            diff=list(set(com).symmetric_difference(set(bom)))
            if (len(diff)==2):
                if (1 in diff) ^ (2 in diff):
           #   if not (1 in diff and 2 in diff):
#                    print(diff)
                    cp=[l for l in diff if l!=1 and l!=2]
                    if (spacing==1):
                        Hamil_site[comb.index(com)][comb.index(bom)]=gauss_orbs[cp[0]][1]
                    elif (spacing==2):
                        Hamil_site[comb.index(com)][comb.index(bom)]=equi_coup#                    
#                        print(com,bom) 


print(Hamil_site)
print(Hamil_site.size,"size",np.count_nonzero(Hamil_site),"non-zero")
print(np.array_equal(Hamil_site,Hamil_site.conjugate().T))

pot_a=[]


# think on this
if (len(mol)==2):
    for i in comb:
        if 1 in i:
            if 2 in i:
                pot_a.append(3)
            else:
                pot_a.append(1)
        elif 2 in i:
            pot_a.append(2)
        else:
            pot_a.append(0)
else:
    for i in comb:
        if 1 in i:
            pot_a.append(1)
        else:
            pot_a.append(0)

print(pot_a)
init=[]
for i in range(len(comb)):
    #print(list(set(gs) & set(comb[i])),"inits")
    init.append(len(list(set(gs) & set(comb[i]))))

init_state=init.index(max(init))+1

print(init_state,'init_state')




np.savetxt('Hamil_site.txt', Hamil_site,fmt='%.3e')
with open('Hamil_site.txt', 'r') as file:
    content = file.read().replace('e', 'd')

with open('Hamil_site.txt', 'w') as file:
    file.write(content)


with open('inpot.txt', 'w') as file:
    for item in inarray:
        file.write(f"{item}\n")

with open('pot_in.txt', 'w') as file:
    for item in pot_a:
        file.write(f"{item}\n")


def overwrite(file_name,line_to_overwrite,new_line_content):


    with open(file_name, 'r') as file:
        lines = file.readlines()

    # Check if the line to overwrite is within the bounds of the file
    if len(lines) >= line_to_overwrite:
        lines[line_to_overwrite - 1] = new_line_content  

    with open(file_name, 'w') as file:
        file.writelines(lines)

KT=0.002
KT_in_K=KT*au2J/kb
KT_input=f"{KT_in_K:.1f}"+"                         !! KT (selected in python).\n"
nquant_input = str(nquant)+"                           !! nquant (Auto-selected).\n"
mol_input=str(len(mol))+"                              !! mol (Auto-selected).\n"
init_state_input=str(init_state)+"                           !! init_state.\n"

overwrite("AFSSH.inp",24,KT_input)
overwrite("AFSSH.inp",16,nquant_input)
overwrite("AFSSH.inp",29,mol_input)
overwrite("AFSSH.inp",30,init_state_input)

