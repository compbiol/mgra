__author__ = 'nikita_kartashov'


class LabeledNode(object):
    def __init__(self, label, left=None, right=None):
        self.__label = label
        self.__left = left
        self.__right = right

    def label(self):
        return self.__label

    def left(self):
        return self.__left

    def right(self):
        return self.__right
