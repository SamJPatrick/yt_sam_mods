import sys
import os
import yt

import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
#from treelib import Tree


CLUMP_DIR = "Britton_sim_data/Clump_data/Clumps"
CLUMP_FILE = "DD0295_clump_info_pisn.h5"
FILE_DICT = {'pisn': "DD0295_clump_info_pisn.h5",
             'hn': "DD0232_clump_info_hn.h5",
             'ccsn': "DD0179_clump_info_pisn.h5"}
MATRIX_NAME = "test_matrix.png"

N_CRIT = unyt_quantity(5.0, 'cm**(-3)')
N_CELLS = 8
VOL_MAX = unyt_quantity(8.0e3, 'pc**3')


class Node:

    added_list = []
    d_matrix = None
    thresh = 0

    
    def __init__(self, clump_num):

        self.clump_num = clump_num
        self.__class__.added_list.append(clump_num)
        self.child_1 = None
        self.child_2 = None


    def add_child(self, clump_num):

        if (self.child_1 == None):
            self.child_1 = Node(clump_num)
            return self.child_1
        elif (self.child_2 == None):
            self.child_2 = Node(clump_num)
            return self.child_2
        else:
            assert RuntimeError


    #def __str__(self):

    #     print(self.clump_num)
    #    if (type(self.child_1) is not None):
    #        print("-->")
    #        print(self.child_1)
    #    if (type(self.child_2) is not None):
    #        print("---->")
    #        print(self.child_2)


    #def make_tree(self, tree):

    #    if (type(self.child_1) is not None):
    #        tree.create_node(self.child_1.clump_num, parent= self.clump_num)
    #        self.child_1.make_tree(tree)
    #    if (type(self.child_2) is not None):
    #        tree.create_node(self.child_2.clump_num, parent= self.clump_num)
    #        self.child_2.make_tree(tree)
    #    return tree


    def append_to_list(self, node_list):

        node_list.append(self.clump_num)
        try :
            node_list = self.child_1.append_to_list(node_list)
        except AttributeError:
            pass
        return node_list
                

    def search_up(self, num_pre):

        n_clumps = self.__class__.d_matrix.shape[0]
        for row in range (num_pre - 1, -1, -1):
            if (type(self).d_matrix[row][self.clump_num] < self.__class__.thresh):
                if (type(self).not_yet_added(row)):
                    node = self.add_child(row)
                    if (self.clump_num == 0):
                        if (row != 0):
                            node.search_up(self.clump_num)
                        node.search_right(self.clump_num)
                    elif (self.clump_num == n_clumps - 1):
                        if (row != 0):
                            node.search_up(self.clump_num)
                        node.search_left(self.clump_num)
                    else :
                        node.search_left(self.clump_num)
                        node.search_right(self.clump_num)
                break

    
    def search_down(self, num_pre):

        n_clumps = self.__class__.d_matrix.shape[0]
        for row in range (num_pre + 1, n_clumps, 1):
            if (type(self).d_matrix[row][self.clump_num] < self.__class__.thresh):
                if (type(self).not_yet_added(row)):
                    node = self.add_child(row)
                    if (self.clump_num == 0):
                        if (row != n_clumps - 1):
                            node.search_down(self.clump_num)
                        node.search_right(self.clump_num)
                    elif (self.clump_num == n_clumps - 1):
                        if (row != n_clumps - 1):
                            node.search_down(self.clump_num)
                        node.search_left(self.clump_num)
                    else :
                        node.search_left(self.clump_num)
                        node.search_right(self.clump_num)
                break
                
                
    def search_left(self, num_pre):

        n_clumps = self.__class__.d_matrix.shape[0]
        for column in range (num_pre - 1, -1, -1):
            if (type(self).d_matrix[self.clump_num][column] < self.__class__.thresh):
                if (type(self).not_yet_added(column)):
                    node = self.add_child(column)
                    if (self.clump_num == 0):
                        if (column != 0):
                            node.search_left(self.clump_num)
                        node.search_down(self.clump_num)
                    elif (self.clump_num == n_clumps - 1):
                        if (column != 0):
                            node.search_left(self.clump_num)
                        node.search_up(self.clump_num)
                    else :
                        node.search_up(self.clump_num)
                        node.search_down(self.clump_num)
                break

                
    def search_right(self, num_pre):

        n_clumps = self.__class__.d_matrix.shape[0]
        for column in range (num_pre + 1, n_clumps, 1):
            if (type(self).d_matrix[self.clump_num][column] < self.__class__.thresh):
                if (type(self).not_yet_added(column)):
                    node = self.add_child(column)
                    if (self.clump_num == 0):
                        if (column != n_clumps - 1):
                            node.search_right(self.clump_num)
                        node.search_down(self.clump_num)
                    elif (self.clump_num == n_clumps - 1):
                        if (column != n_clumps - 1):
                            node.search_right(self.clump_num)
                        node.search_up(self.clump_num)
                    else :
                        node.search_up(self.clump_num)
                        node.search_down(self.clump_num)
                break


    def search_first_right(self):

        n_clumps = self.__class__.d_matrix.shape[0]
        for column in range (0, n_clumps, 1):
            if (type(self).d_matrix[self.clump_num][column] < self.__class__.thresh):
                if (type(self).not_yet_added(column)):
                    node = self.add_child(column)
                    if (self.clump_num == 0):
                        if (column != n_clumps - 1):
                            node.search_right(self.clump_num)
                        node.search_down(self.clump_num)
                    elif (self.clump_num == n_clumps - 1):
                        if (column != n_clumps - 1):
                            node.search_right(self.clump_num)
                        node.search_up(self.clump_num)
                    else :
                        node.search_up(self.clump_num)
                        node.search_down(self.clump_num)
                break


    def produce_tree_list(self):

        return self.__class__.added_list


    @classmethod
    def not_yet_added(cls, number):

        if (number not in cls.added_list):
            return True
        else :
            return False


    @classmethod
    def set_threshold(cls, d_matrix):

        cls.d_matrix = d_matrix
        d_mean = np.mean(d_matrix)
        thresh = d_mean
        #thresh = np.minimum(np.abs(d_mean - np.std(d_matrix)), np.cbrt(d_mean))
        cls.thresh = thresh
        return thresh



def sanatise_clumps(clump_list):

    vol_min = (np.array([clump[('clump', 'total_cells')] for clump in clump_list]) > N_CELLS)
    vol_max = (np.array([clump[('clump', 'volume')] for clump in clump_list]) < VOL_MAX)
    vol_mask = vol_min & vol_max
    dens_mask = (np.array([clump[('clump', 'mean_number_density')] for clump in clump_list]) > N_CRIT)
    sanity_mask = vol_mask & dens_mask

    sane_list = np.array(clump_list)[sanity_mask]
    good_frac = (sum(sanity_mask) / len(sanity_mask)) * 100
    bad_clumps = len(sanity_mask) - sum(sanity_mask)
    print(f"Percentage of clumps remaining is {good_frac:.2f}% with {bad_clumps} bad clumps caught")
    return sane_list





if __name__ == '__main__':

    try :
        star_type = sys.argv[1]
        df = yt.load(os.path.join(CLUMP_DIR, FILE_DICT[star_type]))
    except KeyError, IndexError:
        df = yt.load(os.path.join(CLUMP_DIR, FILE_DICT['pisn']))
    except FileNotFoundError:
        print("Error, CLUMP_DIR needs to be redefined")

    sane_leaves = sanatise_clumps(df.leaves)
    positions = np.array([clump[('clump', 'com')] for clump in sane_leaves])
    masses = np.array([clump[('clump', 'cell_mass')] for clump in sane_leaves])
    n_clumps = len(positions)
    com = np.average(positions, axis=0, weights= masses)
    com_distances = unyt_array([np.linalg.norm(position - com) for position in positions], 'pc')
    mean_com = np.mean(com_distances)
    
    d_matrix = np.zeros((n_clumps, n_clumps))
    for i in range (n_clumps):
        for j in range (n_clumps):
            if (j > i):
                continue
            elif (i == j):
                d_matrix[i][j] = 0
            else :
                value = np.linalg.norm(positions[i] - positions[j])
                d_matrix[i][j] = value
                d_matrix[j][i] = value
                
    pos_start = np.argmin(com_distances)
    master_clump = Node(pos_start)
    mean_inter = unyt_quantity(master_clump.set_threshold(d_matrix), 'pc')
    master_clump.search_first_right()
    #node_list = np.array(master_clump.append_to_list([]))
    node_list = np.array(master_clump.produce_tree_list())
    print(node_list)
    print(com)
    print(mean_inter)
    print(mean_com)

    densities = transpose_unyt([clump[('clump', 'max_number_density')] for clump in sane_leaves])
    metallicities = transpose_unyt([clump[('clump', 'max_metallicity')] for clump in sane_leaves])
    distances = unyt_array([np.linalg.norm(clump[('clump', 'com')] - com) for clump in sane_leaves], 'pc')
    densities_slct = transpose_unyt([clump[('clump', 'max_number_density')] for clump in sane_leaves])[node_list]
    metallicities_slct = transpose_unyt([clump[('clump', 'max_metallicity')] for clump in sane_leaves])[node_list]
    distances_slct = unyt_array([np.linalg.norm(clump[('clump', 'com')] - com) for clump in sane_leaves], 'pc')[node_list]
    
    plt.figure()
    dots_all = plt.scatter(distances, densities, c='blue')
    dots_slct = plt.scatter(distances_slct, densities_slct, c='red')
    plt.xscale('log')
    plt.xlim(5e-1, 5e1)
    plt.xlabel("Distance from COM (pc)")
    plt.yscale('log')
    plt.ylim(1e2, 1e9)
    plt.ylabel("Clump density (cm$^{-3}$)")
    plt.legend([dots_slct], [f'{x:.3f}' for x in [mean_inter]], loc='upper right')
    plt.axvline(x= mean_com, linestyle='--')
    plt.savefig("clump_positions_pisn_distances.png")

    plt.figure()
    dots_all = plt.scatter(distances, metallicities, c='blue')
    dots_slct = plt.scatter(distances_slct, metallicities_slct, c='red')
    plt.xscale('log')
    plt.xlim(5e-1, 5e1)
    plt.xlabel("Distance from COM (pc)")
    plt.ylim(0.0, 7e-4)
    plt.ylabel("Clump metallicity (Z$_{\odot}$)")
    plt.legend([dots_slct], [f'{x:.3f}' for x in [mean_inter]], loc='upper right')
    plt.axvline(x= mean_com, linestyle='--')
    plt.savefig("clump_positions_pisn_metals.png")


    #tree = Tree()
    #tree.create_node(master_clump.clump_num)
    #tree = master_clump.make_tree(tree)
    #tree.show()

