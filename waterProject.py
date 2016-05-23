from htmd import *
from htmd.molecule.util import maxDistance
from htmd.molecule.util import writeVoxels
import itertools
import argparse
import glob

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
    maxdist = maxdist*1.1
    return maxdist

def moving_and_filtering(mol,maxdistance):
    first = True
    for molecule in mol:
        if first:
            molecule.wrap('protein')
            molecule.align('protein')
            first = False
        else:
            molecule.wrap('protein')
            molecule.align('all',refmol=mol[0])
    print("Wrapped and aligned!")
    for molecule in mol:
        molecule.filter('name OH2 and x^2+y^2+z^2<'+str(round(maxdistance)**2))
        #Still outpoints are saved, so we will need another filtering
        molecule.moveBy([maxdistance, maxdistance, maxdistance])
    print("Filtered and moved!")


def grid_generator(mol, cubedistance):
    grid_list = []
    divider = 0
    for molecule in mol:
        grid = np.zeros([cubedistance,cubedistance,cubedistance])
        for atom, frame in itertools.product(*map(range, (molecule.coords.shape[0], molecule.coords.shape[2]))):
            coord_list = [molecule.coords[atom][0][frame], molecule.coords[atom][1][frame], molecule.coords[atom][2][frame]]
            if max(coord_list) > 112:
                continue
            if min(coord_list) < 0:
                continue
            else:
                grid[int(round(coord_list[0]))][int(round(coord_list[1]))][int(round(coord_list[2]))] += 1
        divider += molecule.coords.shape[0]
        grid_list.append(grid)
    final_grid = np.sum(grid_list, axis = 0)/divider
    print("Grid calculated!")
    return final_grid

def isosurface_generator(grid, cubedistance,output):
    for x,y,z in itertools.product(*map(range, (grid.shape[0], grid.shape[1],grid.shape[2]))):
        prob = grid[x][y][z]
        if prob == 0:
            continue
        else:
            grid[x][y][z] = (np.log(grid[x][y][z])) * 0.001987191 * 298 * -1
    
    print("Energy calculated!")
    min_vec= np.array([0,0,0])
    max_vec=np.array([cubedistance,cubedistance,cubedistance])
    res_vec=np.array([1,1,1])
    writeVoxels(grid, output+'.cube', min_vec, max_vec, res_vec) 
    print("Output file "+output+".cube generated on the current directory")     

def main():
    cmdline = cmdline_parser()
    maxdistance, mol_list = maxDistance_calculator(cmdline.infile)
    cubedistance=int(round(maxdistance*2))
    moving_and_filtering(mol_list, cubedistance)
    grid = grid_generator(mol_list, maxdistane)
    isosurface_generator(grid,cubedistance,cmdline.outfile)
    print ("Bye!")

if __name__ == '__main__':
    main()
