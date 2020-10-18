#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import statistics
import argparse
import os
import sys
import random
random.seed(9001)

from copy import deepcopy
from random import randint
from operator import itemgetter

import networkx as nx

__author__ = "Nicolas Bisson"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Nicolas Bisson"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Nicolas Bisson"
__email__ = "Nicolas.bisson@numericable.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fichier_entree):
    """
    Return sequences yield of a fastaq
    """
    with open(fichier_entree) as file:
        for line in file:
            yield next(file).strip()
            next(file)
            next(file)


def cut_kmer(seq, kmer_size):
    """
    Return a k_mer yield from a sequence and kmer_size
    """
    for i in range(len(seq)-kmer_size+1):
        yield seq[i:i+kmer_size]


def build_kmer_dict(fichier, kmer_size):
    """
    Return a kmer_dictionnary from a fastq and kmer_size
    """
    kmer_dict = {}
    for seq in read_fastq(fichier):
        kmers = cut_kmer(seq, kmer_size)
        for kmer in kmers:
            if kmer in seq:
                if kmer not in kmer_dict:
                    kmer_dict[kmer] = 0
                kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    """
    Build a DG graph from a kmer_dictonnary
    """
    graph = nx.DiGraph()
    ele_graph = []
    for key in kmer_dict:
        prefix = key[:-1]
        suffix = key[1:]
        ele_graph.append((prefix, suffix, kmer_dict[key]))
    graph.add_weighted_edges_from(ele_graph)
    return graph


def get_starting_nodes(graph):
    """
    Return the starting nodes of a graph
    """
    starting_nodes = []
    for node in graph.nodes:
        if list(graph.predecessors(node)) == []:  # pas de predecesseur donc noeud entrant
            starting_nodes.append(node)
    return starting_nodes


def get_sink_nodes(graph):
    """
     Return the sink nodes of a graph
    """
    sink_nodes = []
    for node in graph.nodes:
        if list(graph.successors(node)) == []:  # pas de successeur donc noeud sortant
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(graph, starting_nodes, sink_nodes):
    """
    Return contigs from paths of the graph between starting and sink nodes
    """
    contigs = []
    for start_node in starting_nodes:
        for sink_node in sink_nodes:
            for path in nx.all_simple_paths(graph, start_node, sink_node):
                contig = path[0]
                for i in range(1, len(path)):
                    contig += path[i][-1]
                contigs.append((contig, len(contig)))
    return contigs


def fill(text, width=80):
    """
    Split a text to respect the fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_tuple_list, fichier_sortie):
    """
    Save the contigs into a fasta file
    """
    with open(fichier_sortie, "w") as file:
        cpt_contig = 0
        for contig_tuple in contigs_tuple_list:
            file.write(">contig_{} len={}\n".format(cpt_contig,
                                                    contig_tuple[1]))
            file.write(fill(contig_tuple[0]) + "\n")
            cpt_contig += 1


def std(value_list):
    """
    Give the standard deviation between all values of the given list
    """
    return statistics.stdev(value_list)


def path_average_weight(graph, path):
    """
    Return the average weight of a path
    """
    weights = 0
    for i in range(len(path)-1):
        weights += graph.edges[path[i], path[i+1]]["weight"]
    average = weights/(len(path)-1)
    return average


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """
    Remove paths from the graph
    Also remove entry and sink node when parameter set to True
    """
    start = 0
    end = 0
    if delete_entry_node is False:
        start = 1
    if delete_sink_node is False:
        end = 1
    for path in path_list:
        for i in range(start, len(path)-end):
            graph.remove_node(path[i])
    return graph


def select_best_path(graph, path_list, path_length, path_weight,
                     delete_entry_node=False, delete_sink_node=False):
    """
    Selects the best path from a list of path of the graph and remove the other paths
    Parameters:
          graph: a directed graph
          path_list: a list of paths
          path_length: a list of path lengths
          path_weight: a list of average path weights
          delete_entry_node: boolean, True  to remove the entry node
          delete_sink_node: boolean, True to remove the sink node
    """
    max_weight = max(path_weight)
    max_length = max(path_length)
    pos = list(range(0, len(path_list)))
    del_path = []
    pos_tmp = deepcopy(pos)

    if std(path_weight) != 0:
        for i in pos_tmp:
            if path_weight[i] < max_weight:
                del_path.append(path_list[i])
                pos.remove(i)

    path_length = [path_length[i] for i in pos]
    pos_tmp = deepcopy(pos)

    if len(pos) > 1 and std(path_length) != 0:
        for i in pos_tmp:
            if path_length[i] < max_length:
                del_path.append(path_list[i])
                pos.remove(i)

    if len(pos) > 1:
        pos.remove(randint(0, len(pos)))
        for i in pos:
            del_path.append(path_list[i])

    return remove_paths(graph, del_path, delete_entry_node, delete_sink_node)


def solve_bubble(graph, predecessor_node, successor_node):
    """
    Remove a bubble from the graph
    """
    path_list = []
    path_length = []
    path_weight = []

    for path in nx.all_simple_paths(graph, predecessor_node, successor_node):
        path_list.append(path)
        path_weight.append(path_average_weight(graph, path))
        path_length.append(len(path))

    return select_best_path(graph, path_list, path_length, path_weight)


def simplify_bubbles(graph):
    """
    Remove all bubbles from the graph
    """
    bubbles = []
    for node in graph.nodes:
        middle = list(graph.predecessors(node))
        if len(middle) >= 2:
            c_ancestor = nx.lowest_common_ancestor(graph, middle[0], middle[1])
            bubbles.append([c_ancestor, node])

    for bubble in bubbles:
        graph = solve_bubble(graph, bubble[0], bubble[1])

    return graph


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    args = get_arguments()
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    start_nodes = get_starting_nodes(graph)
    end_nodes = get_sink_nodes(graph)
    graph = simplify_bubbles(graph)
    #graph = solve_entry_tips(graph, list_start_nodes)
    #graph = solve_out_tips(graph, list_end_nodes)
    list_contigs = get_contigs(graph, start_nodes, end_nodes)
    save_contigs(list_contigs, args.output_file)


if __name__ == '__main__':
    main()
