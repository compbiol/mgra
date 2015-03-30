__author__ = 'nikita_kartashov'

from Queue import Queue


class UnrootedBinaryTree(object):
    def __init__(self, artificial_root_node):
        # Vertices are stored as a list of lists, so if node x is connected to y
        # y in vertices[x] = x in vertices[y] = True
        self.__vertices = []
        # ith element is true if ith vertex is a leaf
        self.__leaves = {}
        self.__process_artificial_root(artificial_root_node)
        self.__leaf_distances = self.__calculate_pairwise_leaf_distances()

    def __calculate_pairwise_leaf_distances(self):
        distances = [[0 for _ in self.__vertices] for _ in self.__vertices]
        for start_leaf in self.__leaves.iterkeys():
            for end_leaf in self.__leaves.iterkeys():
                self.calculate_distance_between_two(distances, start_leaf, end_leaf)
        leaf_distances = {}
        for left_leaf, left_genome in self.__leaves.iteritems():
            leaf_distances[left_genome] = {}
            for right_leaf, right_genome in self.__leaves.iteritems():
                leaf_distances[left_genome][right_genome] = distances[left_leaf][right_leaf]
        return leaf_distances

    def pretty_branches(self, branches):
        complete_label = set(self.leaf_labels())
        return ['{0} + {1}'.format(', '.join(branch),
                                   ', '.join(complete_label.difference(branch))) for branch in branches]

    def branches_when_rooted_at(self, root_label):
        root_node, _ = filter(lambda (node, label): label == root_label, self.__leaves.iteritems())[0]
        labels = [None for _ in self.__vertices]
        _, tree_labels = self.__label_nodes(root_node, labels)
        return {frozenset(element for element in s) for s in tree_labels}

    def compare_branches(self, other_tree, rooted_at=None):
        if rooted_at is None:
            k, v = next(self.__leaves.iteritems())
            rooted_at = v
        return self.branches_when_rooted_at(rooted_at).difference(other_tree.branches_when_rooted_at(rooted_at))

    def compare_distances(self, other_tree):
        assert (set(self.leaf_labels()) == set(other_tree.leaf_labels()))
        distance_difference = 0
        for label in self.leaf_labels():
            for other_label in self.leaf_labels():
                distance_difference += \
                    abs(self.leaf_distances()[label][other_label] - other_tree.leaf_distances()[label][other_label])
        return distance_difference * 1.0 / len(self.leaf_labels())


    def __label_nodes(self, current_node, labels, visited=None):
        if visited is None:
            visited = [False for _ in self.__vertices]
        visited[current_node] = True
        current_node_labels = set() if len(self.__vertices[current_node]) != 1 else {self.__leaves[current_node]}
        current_tree_labels = [current_node_labels]
        nodes_to_visit = filter(lambda n: not visited[n], self.__vertices[current_node])
        for node in nodes_to_visit:
            new_node_labels, new_tree_labels = self.__label_nodes(node, labels, visited)
            current_node_labels.update(new_node_labels)
            current_tree_labels.extend(new_tree_labels)

        labels[current_node] = current_node_labels
        return current_node_labels, current_tree_labels

    def leaf_distances(self):
        return self.__leaf_distances

    def leaves(self):
        return self.__leaves.keys()

    def leaf_labels(self):
        return self.__leaves.values()

    def calculate_distance_between_two(self, distances, start, end, current=None, distance_so_far=0, visited=None):
        if current is None:
            current = start
        if visited is None:
            visited = [False for _ in self.__vertices]
        if current == end:
            distances[start][end] = distance_so_far
            distances[end][start] = distance_so_far
        visited[current] = True
        for node in filter(lambda n: not visited[n], self.__vertices[current]):
            self.calculate_distance_between_two(distances, start, end, node, distance_so_far + 1, visited)

    def __process_artificial_root(self, artificial_root_node):
        nodes_to_process = Queue()
        for node in artificial_root_node.child_nodes():
            nodes_to_process.put((node, None))

        while not nodes_to_process.empty():
            node, parent = nodes_to_process.get()
            self.__vertices.append([])
            current_node = len(self.__vertices) - 1
            for child in node.child_nodes():
                nodes_to_process.put((child, current_node))
            if not node.child_nodes():
                self.__leaves[current_node] = node.taxon.__str__()
            if parent is not None:
                self.__vertices[parent].append(current_node)
                self.__vertices[current_node].append(parent)
        self.__vertices[0].append(1)
        self.__vertices[1].append(0)