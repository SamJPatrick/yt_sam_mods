import sys
import os
import yt

import numpy as np
from unyt import unyt_quantity, unyt_array
from yt.extensions.sam_mods.unyt_funcs import transpose_unyt
#from treelib import Tree


N_CRIT = unyt_quantity(5.0, 'cm**(-3)')
N_CELLS = 8
VOL_MAX = unyt_quantity(8.0e3, 'pc**3')



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



class Node:

    added_list = []
    d_matrix = None
    thresh = 0

    
    def __init__(self, clump_num):

        self.clump_num = clump_num
        self.__class__.added_list.append(clump_num)
        self.child_1 = None
        self.child_2 = None
        self.child_3 = None


    def add_child(self, clump_num):

        if (self.child_1 == None):
            self.child_1 = Node(clump_num)
            return self.child_1
        elif (self.child_2 == None):
            self.child_2 = Node(clump_num)
            return self.child_2
        elif (self.child_3 == None):
            self.child_3 = Node(clump_num)
            return self.child_3
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
                    if (self.clump_num != 0):
                        node.search_left(self.clump_num)
                    if (self.clump_num != n_clumps - 1):
                        node.search_right(self.clump_num)
                    if (row != 0):
                        node.search_up(self.clump_num)
                break
    

    
    def search_down(self, num_pre):
    
        n_clumps = self.__class__.d_matrix.shape[0]
        for row in range (num_pre + 1, n_clumps, 1):
            if (type(self).d_matrix[row][self.clump_num] < self.__class__.thresh):
                if (type(self).not_yet_added(row)):
                    node = self.add_child(row)
                    if (self.clump_num != 0):
                        node.search_left(self.clump_num)
                    if (self.clump_num != n_clumps - 1):
                        node.search_right(self.clump_num)
                    if (row != n_clumps - 1):
                        node.search_down(self.clump_num)
                break
            

    
    def search_left(self, num_pre):
    
        n_clumps = self.__class__.d_matrix.shape[0]
        for column in range (num_pre - 1, -1, -1):
            if (type(self).d_matrix[self.clump_num][column] < self.__class__.thresh):
                if (type(self).not_yet_added(column)):
                    node = self.add_child(column)
                    if (self.clump_num != 0):
                        node.search_up(self.clump_num)
                    if (self.clump_num != n_clumps - 1):
                        node.search_down(self.clump_num)
                    if (column != 0):
                        node.search_left(self.clump_num)
                break
    

    
    def search_right(self, num_pre):
    
        n_clumps = self.__class__.d_matrix.shape[0]
        for column in range (num_pre + 1, n_clumps, 1):
            if (type(self).d_matrix[self.clump_num][column] < self.__class__.thresh):
                if (type(self).not_yet_added(column)):
                    node = self.add_child(column)
                    if (self.clump_num != 0):
                        node.search_up(self.clump_num)
                    if (self.clump_num != n_clumps - 1):
                        node.search_down(self.clump_num)
                    if (column != n_clumps -1):
                        node.search_right(self.clump_num)
                break
    

    
    def search_first_right(self):
    
        n_clumps = self.__class__.d_matrix.shape[0]
        for column in range (0, n_clumps, 1):
            if (type(self).d_matrix[self.clump_num][column] < self.__class__.thresh):
                if (type(self).not_yet_added(column)):
                    node = self.add_child(column)
                    if (self.clump_num != 0):
                        node.search_up(self.clump_num)
                    if (self.clump_num != n_clumps - 1):
                        node.search_down(self.clump_num)
                    if (column != n_clumps -1):
                        node.search_right(self.clump_num)
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
    def set_threshold(cls, positions, masses):

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
        thresh = np.mean(d_matrix[np.triu_indices(n_clumps, 0)])
        cls.d_matrix = d_matrix
        cls.thresh = thresh
        cls.pos_start = pos_start

        print(com)
        print(thresh)
        print(mean_com)

        return (pos_start, com, thresh, mean_com)

