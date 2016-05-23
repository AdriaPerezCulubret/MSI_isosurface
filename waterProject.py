from htmd import *
from htmd.molecule.util import maxDistance
import itertools
import argparse
import glob




#mol = Molecule("pub.htmd.org/CXCL12-confAnalysis/1x14/structure.pdb")
#mol.read('pub.htmd.org/CXCL12-confAnalysis/1x14/traj.xtc')
#
#
#
#
##--- compute max distance ---
#from htmd.molecule.util import maxDistance
#D = maxDistance(mol, 'all')
#maxdistance = D*1.1
#maxdistance
#

#mol.wrap('protein')
#mol.align('protein')
#mol.filter('name OH2 and x^2+y^2+z^2<'+str(round(maxdistance)**2))
#mol.moveBy([maxdistance, maxdistance, maxdistance])
#
#
#cubedistance=int(round(maxdistance*2))
#
#
#
#grid = np.zeros([cubedistance,cubedistance,cubedistance])
#
#
#for atom, frame in itertools.product(*map(range, (mol.coords.shape[0], mol.coords.shape[2]))):
#    coord_list = [mol.coords[atom][0][frame], mol.coords[atom][1][frame], mol.coords[atom][2][frame]]
#    if max(coord_list) > 112:
#        continue
#    if min(coord_list) < 0:
#        continue
#    else:
#        grid[int(round(coord_list[0]))][int(round(coord_list[1]))][int(round(coord_list[2]))] += 1
#print("Done!")
#
## --- making count matrix --- END
## --- energy and writeVoxels
#grid_div = grid/7066
#grid_div.shape
#
#
#
#for x,y,z in itertools.product(*map(range, (grid_div.shape[0], grid_div.shape[1],grid_div.shape[2]))):
#    prob = grid_div[x][y][z]
#    if prob == 0:
#        continue
#    else:
#        grid_div[x][y][z] = (np.log(grid_div[x][y][z])) * 0.001987191 * 298 * -1 
#        
#
#
#from htmd.molecule.util import writeVoxels
#min_vec= np.array([0,0,0])
#max_vec=np.array([112,112,112])
#res_vec=np.array([1,1,1])
#writeVoxels(grid_div,'isosurf.cube', min_vec, max_vec, res_vec )
## --- energy and writeVoxels --- END
#
def cmdline_parser():

    parser = argparse.ArgumentParser(description="""This program gets a pdb and a cif structure and returns a pdb file with the protein placed 
                                                    within two layers of dummy atoms corresponding to a lipidic membrane.""")

    parser.add_argument('-i','--input',
                        dest = "infile",
                        action = "store",
                        required = True,
                        default = None,
                        help = "Input file or directory with the .pdb and .xtc files")

    parser.add_argument('-o','--output',
                        dest = "outfile",
                        action = "store",
                        default = "isosurface",
                        help = "Name of the .cube output file")


    options = parser.parse_args()

    return options

def maxDistance_calculator(infile):
    counter = 0
    mol = []
    distance_list=[]
    for file in glob.glob(infile +'/*/*.pdb'):
        mol.append(Molecule(file))
        mol[counter].read(glob.glob(infile+'/*/*.xtc')[counter])
        distance_list.append(maxDistance(mol[counter], 'all'))
        counter +=1
        print(file + " Done!")
    maxdist = max(distance_list)
    maxdist = D*1.1
    return maxdist

def grid_generator():
    pass

def isosurface_generator():
    pass


cmdline = cmdline_parser()
maxdistance=maxDistance_calculator(cmdline.infile)
print (maxdistance)