__author__ = 'nikita_kartashov'

from sys import argv

import dendropy as dp

from labeled_node import LabeledNode
from unrooted_tree import UnrootedBinaryTree


taxa = dp.TaxonSet()


def nt(s):
    return dp.Tree.get_from_string(s, 'newick', taxon_set=taxa)


def get_taxa(node):
    taxon_set = set()
    if node.child_nodes():
        for child in node.child_nodes():
            result = get_taxa(child)
            taxon_set = taxon_set.union(result)
    else:
        taxon_set.add(node.taxon.__str__())
    return taxon_set


def build_labeled_node(node):
    taxon_set = set()
    if node.child_nodes():
        left_child, right_child = node.child_nodes()
        left_node, right_node = build_labeled_node(left_child), build_labeled_node(right_child)
        taxon_set = taxon_set.union(left_node.label())
        taxon_set = taxon_set.union(right_node.label())
        return LabeledNode(taxon_set, left_node, right_node)
    else:
        taxon_set.add(node.taxon.__str__())
    return LabeledNode(taxon_set)


mgra_tree = nt('(T,((A,(S,(N,M))),((U,(E,P)),D)));')
real_tree = nt('(T,(D,(((U,P),E),((A,S),(M,N)))));')
test_tree1 = nt('((A, D), (B, C));')
test_tree2 = nt('((A, B), (C, D));')


def main():
    global left_tree, right_tree, root, t1, t2
    if len(argv) < 3:
        print("Need to input 2 trees in newick format")
        exit(1)
    left_tree = nt(argv[1].strip('\''))
    right_tree = nt(argv[2].strip('\''))
    root = mgra_tree.nodes()[0]
    t1 = UnrootedBinaryTree(left_tree.nodes()[0])
    t2 = UnrootedBinaryTree(right_tree.nodes()[0])
    print('\n'.join(t1.pretty_branches(t1.compare_branches(t2))))
    print('Distance difference: {0}'.format(t1.compare_distances(t2)))


if __name__ == '__main__':
    main()

