import csv
import random
import math


class Model:

    # instance variables
    def __init__(self):
        self.allNodes = []
        self.customers = []
        self.matrix = []
        self.capacity = -1

    def BuildModel(self):
        random.seed(1)
        depot = Node(0, 35, 35, 0, 0)
        self.allNodes.append(depot)

        self.capacity = 200
        file = open('input.csv')
        csvreader = csv.reader(file)

        for row in csvreader:
            node = Node(int(row[0]), int(row[1]), int(row[2]), int(row[3]), int(row[4]))
            self.allNodes.append(node)
            self.customers.append(node)

        rows = len(self.allNodes)
        self.matrix = [[0.0 for x in range(rows)] for y in range(rows)]

        for i in range(0, len(self.allNodes)):
            for j in range(0, len(self.allNodes)):
                a = self.allNodes[i]
                b = self.allNodes[j]
                dist = math.sqrt(math.pow(a.x - b.x, 2) + math.pow(a.y - b.y, 2))
                self.matrix[i][j] = dist


class Node:
    def __init__(self, idd, xx, yy, dem, untime):
        self.x = xx
        self.y = yy
        self.ID = idd
        self.demand = dem
        self.isRouted = False
        self.positionInRoute=0
        self.untime = untime
        self.waitingtime = 0
        self.cost_up_to_here = 0


class Route:
    def __init__(self, dp, cap):
        self.sequenceOfNodes = []
        self.sequenceOfNodes.append(dp)
        self.cost = 0
        self.capacity = cap
        self.load = 0
