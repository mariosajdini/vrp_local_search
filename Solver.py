import copy

from VRP_Model import *
from SolutionDrawer import *


class Solution:
    def __init__(self):
        self.cost = 0.0
        self.routes = []


class RelocationMove(object):
    def __init__(self):
        self.originRoutePosition = None
        self.targetRoutePosition = None
        self.originNodePosition = None
        self.targetNodePosition = None
        self.costChangeOriginRt = None
        self.costChangeTargetRt = None
        self.moveCost = None

    def Initialize(self):
        self.originRoutePosition = None
        self.targetRoutePosition = None
        self.originNodePosition = None
        self.targetNodePosition = None
        self.costChangeOriginRt = None
        self.costChangeTargetRt = None
        self.moveCost = 10 ** 9


class SwapMove(object):
    def __init__(self):
        self.positionOfFirstRoute = None
        self.positionOfSecondRoute = None
        self.positionOfFirstNode = None
        self.positionOfSecondNode = None
        self.costChangeFirstRt = None
        self.costChangeSecondRt = None
        self.moveCost = None

    def Initialize(self):
        self.positionOfFirstRoute = None
        self.positionOfSecondRoute = None
        self.positionOfFirstNode = None
        self.positionOfSecondNode = None
        self.costChangeFirstRt = None
        self.costChangeSecondRt = None
        self.moveCost = 10 ** 9


class CustomerInsertion(object):
    def __init__(self):
        self.customer = None
        self.route = None
        self.cost = 10 ** 9


class CustomerInsertionAllPositions(object):
    def __init__(self):
        self.customer = None
        self.route = None
        self.insertionPosition = None
        self.cost = 10 ** 9


class TwoOptMove(object):
    def __init__(self):
        self.positionOfFirstRoute = None
        self.positionOfSecondRoute = None
        self.positionOfFirstNode = None
        self.positionOfSecondNode = None
        self.moveCost = None

    def Initialize(self):
        self.positionOfFirstRoute = None
        self.positionOfSecondRoute = None
        self.positionOfFirstNode = None
        self.positionOfSecondNode = None
        self.moveCost = 10 ** 9


class Solver:
    def __init__(self, m):
        self.allNodes = m.allNodes
        self.customers = m.customers
        self.depot = m.allNodes[0]
        self.distanceMatrix = m.matrix
        self.capacity = m.capacity
        self.used = {self.depot.ID}
        self.sol = None
        self.bestSolution = None

    def solve(self):
        self.SetRoutedFlagToFalseForAllCustomers()
        self.ApplyNearestNeighborMethod()
        # self.MinimumInsertions()
        self.ReportSolution(self.sol)
        self.LocalSearch(0)
        return self.sol

    def SetRoutedFlagToFalseForAllCustomers(self):
        for i in range(0, len(self.customers)):
            self.customers[i].isRouted = False

    def find_node(self, r, n):
        dist = self.distanceMatrix[n.ID]
        nearest_v = 10000000
        nearest_index = 0
        for i in range(0, len(dist)):
            if dist[i] < nearest_v and i not in self.used and r.load + self.allNodes[i].demand < r.capacity:
                nearest_index = i
                nearest_v = dist[i]
        self.used.add(nearest_index)
        return nearest_index, nearest_v

    def ApplyNearestNeighborMethod(self):
        modelIsFeasible = True
        self.sol = Solution()
        for i in range(0, 14):
            self.sol.routes.append(Route(self.depot, self.capacity))
        j = 0
        for i in range(1, len(self.allNodes)):
            route = self.sol.routes[j % 14]
            node = route.sequenceOfNodes[-1]
            nearest_possible, value = self.find_node(route, node)
            self.allNodes[nearest_possible].isRouted = True
            if (len(route.sequenceOfNodes) == 1):
                self.allNodes[nearest_possible].waitingtime = value
            else:
                self.allNodes[nearest_possible].waitingtime = route.sequenceOfNodes[-1].waitingtime + value
            route.load += self.allNodes[nearest_possible].demand
            route.cost += self.allNodes[nearest_possible].waitingtime
            self.allNodes[nearest_possible].waitingtime += 10
            route.sequenceOfNodes.append(self.allNodes[nearest_possible])
            self.allNodes[nearest_possible].positionInRoute = len(route.sequenceOfNodes) - 1
            j += 1

        if (modelIsFeasible == False):
            print('FeasibilityIssue')
            # reportSolution

    def Always_keep_an_empty_route(self):
        if len(self.sol.routes) == 0:
            rt = Route(self.depot, self.capacity)
            self.sol.routes.append(rt)
        else:
            rt = self.sol.routes[-1]
            if len(rt.sequenceOfNodes) > 2:
                rt = Route(self.depot, self.capacity)
                self.sol.routes.append(rt)

    def LocalSearch(self, operator):
        self.bestSolution = self.cloneSolution(self.sol)
        terminationCondition = False
        localSearchIterator = 0

        rm = RelocationMove()
        sm = SwapMove()
        top = TwoOptMove()

        while terminationCondition is False:

            self.InitializeOperators(rm, sm, top)
            SolDrawer.draw(localSearchIterator, self.sol, self.allNodes)

            # Relocations
            if operator == 0:
                self.FindBestRelocationMove(rm)
                if rm.originRoutePosition is not None:
                    if rm.moveCost < 0:
                        self.ApplyRelocationMove(rm)
                    else:
                        terminationCondition = True
            # Swaps
            elif operator == 1:
                self.FindBestSwapMove(sm)
                if sm.positionOfFirstRoute is not None:
                    if sm.moveCost < 0:
                        self.ApplySwapMove(sm)
                    else:
                        terminationCondition = True
            # two opt
            elif operator == 2:
                self.FindBestTwoOptMove(top)
                if top.positionOfFirstRoute is not None:
                    if top.moveCost < 0:
                        self.ApplyTwoOptMove(top)
                    else:
                        terminationCondition = True

            if self.sol.cost < self.bestSolution.cost:
                self.bestSolution = self.cloneSolution(self.sol)

            localSearchIterator = localSearchIterator + 1
            print(localSearchIterator, self.sol.cost)

        self.sol = self.bestSolution

    def cloneRoute(self, rt: Route):
        cloned = Route(self.depot, self.capacity)
        cloned.cost = rt.cost
        cloned.load = rt.load
        cloned.sequenceOfNodes = rt.sequenceOfNodes.copy()
        return cloned

    def cloneSolution(self, sol: Solution):
        cloned = Solution()
        for i in range(0, len(sol.routes)):
            rt = sol.routes[i]
            clonedRoute = self.cloneRoute(rt)
            cloned.routes.append(clonedRoute)
        cloned.cost = self.sol.cost
        return cloned

    def FindBestRelocationMove(self, rm):
        for originRouteIndex in range(0, len(self.sol.routes)):
            rt1: Route = self.sol.routes[originRouteIndex]
            for originNodeIndex in range(1, len(rt1.sequenceOfNodes) - 1):
                for targetRouteIndex in range(0, len(self.sol.routes)):
                    rt2: Route = self.sol.routes[targetRouteIndex]
                    for targetNodeIndex in range(0, len(rt2.sequenceOfNodes) - 1):

                        if originRouteIndex == targetRouteIndex and (
                                targetNodeIndex == originNodeIndex or targetNodeIndex == originNodeIndex - 1):
                            continue

                        A: Node = rt1.sequenceOfNodes[originNodeIndex - 1]
                        B: Node = rt1.sequenceOfNodes[originNodeIndex]
                        cut_of_AB = len(rt1.sequenceOfNodes) - A.positionInRoute - 1
                        C: Node = rt1.sequenceOfNodes[originNodeIndex + 1]
                        cut_of_BC = len((rt1.sequenceOfNodes)) - B.positionInRoute - 1
                        F: Node = rt2.sequenceOfNodes[targetNodeIndex]
                        G: Node = rt2.sequenceOfNodes[targetNodeIndex + 1]
                        cut_of_FG = len(rt2.sequenceOfNodes) - F.positionInRoute - 1
                        oldCost = rt1.cost + rt2.cost
                        if rt1 != rt2:
                            if rt2.load + B.demand > rt2.capacity:
                                continue
                        A: Node = rt1.sequenceOfNodes[originNodeIndex - 1]
                        B: Node = rt1.sequenceOfNodes[originNodeIndex]
                        cut_of_AB = len(rt1.sequenceOfNodes) - A.positionInRoute - 1
                        C: Node = rt1.sequenceOfNodes[originNodeIndex + 1]
                        cut_of_BC = len((rt1.sequenceOfNodes)) - B.positionInRoute - 1
                        F: Node = rt2.sequenceOfNodes[targetNodeIndex]
                        G: Node = rt2.sequenceOfNodes[targetNodeIndex + 1]
                        cut_of_FG = len(rt2.sequenceOfNodes) - F.positionInRoute - 1
                        oldCost = rt1.cost + rt2.cost
                        cost_cummulative_removed_of_rt1 = 0
                        for i in range(A.positionInRoute, 0, -1):
                            pred = i - 1
                            n1 = rt1.sequenceOfNodes[i]
                            n2 = rt1.sequenceOfNodes[pred]
                            cost_cummulative_removed_of_rt1 += self.distanceMatrix[n1.ID][n2.ID]
                        cost_removed = cut_of_AB * self.distanceMatrix[A.ID][B.ID] + cut_of_BC * \
                                       self.distanceMatrix[B.ID][C.ID]
                        cost_cummulative_removed_of_rt1 += (len(rt1.sequenceOfNodes) - 2) * 10
                        originRtCostChange = -cost_cummulative_removed_of_rt1 - \
                                             cost_removed + cut_of_BC * self.distanceMatrix[A.ID][C.ID]
                        cost_cummulative_add_of_rt2 = 0
                        for i in range(F.positionInRoute, 0, -1):
                            pred = i - 1
                            n1 = rt2.sequenceOfNodes[i]
                            n2 = rt2.sequenceOfNodes[pred]
                            cost_cummulative_add_of_rt2 += self.distanceMatrix[n1.ID][n2.ID]
                        cost_cummulative_add_of_rt2 += (len(rt2.sequenceOfNodes) - 1) * 10
                        targetRtCostChange = cost_cummulative_add_of_rt2 - cut_of_FG * self.distanceMatrix[F.ID][G.ID] \
                                             + cut_of_FG * self.distanceMatrix[B.ID][G.ID] \
                                             + (cut_of_FG + 1) * self.distanceMatrix[F.ID][B.ID]
                        moveCost = targetRtCostChange + originRtCostChange

                        # # easy way
                        # test_rt1 = copy.deepcopy(rt1)
                        # test_rt2 = copy.deepcopy(rt2)
                        # if test_rt1 == test_rt2:
                        #      del test_rt1.sequenceOfNodes[originNodeIndex]
                        #      if originNodeIndex < targetNodeIndex:
                        #          test_rt2.sequenceOfNodes.insert(targetNodeIndex, B)
                        #      else:
                        #          test_rt2.sequenceOfNodes.insert(targetNodeIndex + 1, B)
                        # else:
                        #      del test_rt1.sequenceOfNodes[originNodeIndex]
                        #      test_rt2.sequenceOfNodes.insert(targetNodeIndex + 1, B)
                        #
                        # c1, l1 = self.calculate_route_details(test_rt1.sequenceOfNodes)
                        # test_rt1.cost = c1
                        # test_rt1.load = l1
                        # c2, l2 = self.calculate_route_details(test_rt2.sequenceOfNodes)
                        # test_rt2.cost = c2
                        # test_rt2.load = l2
                        # new_cost = c1 + c2
                        # originRtCostChange = c1 - rt1.cost
                        # targetRtCostChange = c2 - rt2.cost
                        # moveCost = new_cost - oldCost

                        if moveCost < rm.moveCost:
                            self.StoreBestRelocationMove(originRouteIndex, targetRouteIndex, originNodeIndex,
                                                         targetNodeIndex, moveCost, originRtCostChange,
                                                         targetRtCostChange, rm)

    def FindBestSwapMove(self, sm):
        for firstRouteIndex in range(0, len(self.sol.routes)):
            rt1: Route = self.sol.routes[firstRouteIndex]
            for secondRouteIndex in range(firstRouteIndex, len(self.sol.routes)):
                rt2: Route = self.sol.routes[secondRouteIndex]
                for firstNodeIndex in range(1, len(rt1.sequenceOfNodes) - 1):
                    startOfSecondNodeIndex = 1
                    if rt1 == rt2:
                        startOfSecondNodeIndex = firstNodeIndex + 1
                    for secondNodeIndex in range(startOfSecondNodeIndex, len(rt2.sequenceOfNodes) - 1):
                        a1: Node = rt1.sequenceOfNodes[firstNodeIndex - 1]
                        cut_of_a1_b1 = len(rt1.sequenceOfNodes) - a1.positionInRoute - 1
                        b1 = rt1.sequenceOfNodes[firstNodeIndex]
                        c1 = rt1.sequenceOfNodes[firstNodeIndex + 1]
                        cut_of_b1_c1 = len(rt1.sequenceOfNodes) - b1.positionInRoute - 1

                        a2 = rt2.sequenceOfNodes[secondNodeIndex - 1]
                        cut_of_a2_b2 = len(rt2.sequenceOfNodes) - a2.positionInRoute - 1
                        b2 = rt2.sequenceOfNodes[secondNodeIndex]
                        cut_of_b2_c2 = len(rt2.sequenceOfNodes) - b2.positionInRoute - 1
                        c2 = rt2.sequenceOfNodes[secondNodeIndex + 1]
                        old_cst = None
                        moveCost = None
                        costChangeFirstRoute = None
                        costChangeSecondRoute = None
                        if rt1 == rt2:
                            old_cst = rt1
                            if firstNodeIndex == secondNodeIndex - 1:
                                # case of consecutive nodes swap
                                costRemoved = cut_of_a1_b1 * self.distanceMatrix[a1.ID][b1.ID] + cut_of_b1_c1 * \
                                              self.distanceMatrix[b1.ID][b2.ID] + \
                                              cut_of_b2_c2 * self.distanceMatrix[b2.ID][c2.ID]
                                costAdded = cut_of_a1_b1 * self.distanceMatrix[a1.ID][b2.ID] + cut_of_b1_c1 * \
                                            self.distanceMatrix[b2.ID][b1.ID] + \
                                            cut_of_b2_c2 * self.distanceMatrix[b1.ID][c2.ID]
                                moveCost = costAdded - costRemoved
                            else:

                                costRemoved1 = cut_of_a1_b1 * self.distanceMatrix[a1.ID][b1.ID] + cut_of_b1_c1 * \
                                               self.distanceMatrix[b1.ID][c1.ID]
                                costAdded1 = cut_of_a1_b1 * self.distanceMatrix[a1.ID][b2.ID] + cut_of_a1_b1 * \
                                             self.distanceMatrix[b2.ID][c1.ID]
                                costRemoved2 = cut_of_a2_b2 * self.distanceMatrix[a2.ID][b2.ID] + cut_of_b2_c2 * \
                                               self.distanceMatrix[b2.ID][c2.ID]
                                costAdded2 = cut_of_a2_b2 * self.distanceMatrix[a2.ID][b1.ID] + cut_of_b2_c2 * \
                                             self.distanceMatrix[b1.ID][c2.ID]
                                moveCost = costAdded1 + costAdded2 - (costRemoved1 + costRemoved2)
                        else:
                            if rt1.load - b1.demand + b2.demand > self.capacity:
                                continue
                            if rt2.load - b2.demand + b1.demand > self.capacity:
                                continue
                            old_cst = rt1.cost + rt2.cost
                            costRemoved1 = cut_of_a1_b1 * self.distanceMatrix[a1.ID][b1.ID] + cut_of_b1_c1 * \
                                           self.distanceMatrix[b1.ID][c1.ID]
                            costAdded1 = cut_of_a1_b1 * self.distanceMatrix[a1.ID][b2.ID] + cut_of_b1_c1 * \
                                         self.distanceMatrix[b2.ID][c1.ID]
                            costRemoved2 = cut_of_a2_b2 * self.distanceMatrix[a2.ID][b2.ID] + cut_of_b2_c2 * \
                                           self.distanceMatrix[b2.ID][c2.ID]
                            costAdded2 = cut_of_a2_b2 * self.distanceMatrix[a2.ID][b1.ID] + cut_of_b2_c2 * \
                                         self.distanceMatrix[b1.ID][c2.ID]

                            costChangeFirstRoute = costAdded1 - costRemoved1
                            costChangeSecondRoute = costAdded2 - costRemoved2

                            moveCost = costAdded1 + costAdded2 - (costRemoved1 + costRemoved2)
                        # test_rt1 = copy.deepcopy(rt1)
                        # test_rt2 = copy.deepcopy(rt2)
                        # test_rt1.sequenceOfNodes[sm.positionOfFirstNode] = b2
                        # rt2.sequenceOfNodes[sm.positionOfSecondNode] = b1
                        # test_rt2.sequenceOfNodes[sm.positionOfFirstNode].positionInRoute = sm.positionOfFirstNode
                        # rt1.sequenceOfNodes[sm.positionOfSecondNode].positionInRoute = sm.positionOfSecondNode
                        #
                        if moveCost < sm.moveCost:
                            self.StoreBestSwapMove(firstRouteIndex, secondRouteIndex, firstNodeIndex, secondNodeIndex,
                                                   moveCost, costChangeFirstRoute, costChangeSecondRoute, sm)

    def ApplyRelocationMove(self, rm: RelocationMove):
        oldCost = self.CalculateTotalCost(self.sol)

        originRt = self.sol.routes[rm.originRoutePosition]
        targetRt = self.sol.routes[rm.targetRoutePosition]

        B = originRt.sequenceOfNodes[rm.originNodePosition]

        if originRt == targetRt:
            del originRt.sequenceOfNodes[rm.originNodePosition]
            if (rm.originNodePosition < rm.targetNodePosition):
                targetRt.sequenceOfNodes.insert(rm.targetNodePosition, B)
            else:
                targetRt.sequenceOfNodes.insert(rm.targetNodePosition + 1, B)

            originRt.cost += rm.moveCost
        else:
            del originRt.sequenceOfNodes[rm.originNodePosition]
            targetRt.sequenceOfNodes.insert(rm.targetNodePosition + 1, B)
            originRt.cost += rm.costChangeOriginRt
            targetRt.cost += rm.costChangeTargetRt
            originRt.load -= B.demand
            targetRt.load += B.demand

        self.sol.cost += rm.moveCost

        newCost = self.CalculateTotalCost(self.sol)
        self.Update_Positions_In_Routes(originRt, targetRt)
        # debuggingOnly
        if abs((newCost - oldCost) - rm.moveCost) > 0.0001:
            print("Tryed to do relocation from Node " + str(rm.originNodePosition) + " of Route " + str(
                rm.originRoutePosition) + " to Node " + str(rm.targetNodePosition) + " of Route " + str(
                rm.targetRoutePosition))
            self.ReportSolution(self.sol)
            print('Cost Issue')

    def Update_Positions_In_Routes(self, rt1, rt2):
        for i in range(0, len(rt1.sequenceOfNodes)):
            n1: Node = rt1.sequenceOfNodes[i]
            n1.positionInRoute = i
        for i in range(0, len(rt2.sequenceOfNodes)):
            n1: Node = rt2.sequenceOfNodes[i]
            n1.positionInRoute = i

    def ApplySwapMove(self, sm):
        print("swaped the " + str(sm.positionOfFirstNode) + " node of the " +
              str(sm.positionOfFirstRoute) + " route " + " with the " + str(
            sm.positionOfSecondNode) + " node of the " + str(sm.positionOfSecondRoute) + " route")
        oldCost = self.CalculateTotalCost(self.sol)
        rt1 = self.sol.routes[sm.positionOfFirstRoute]
        rt2 = self.sol.routes[sm.positionOfSecondRoute]
        b1: Node = rt1.sequenceOfNodes[sm.positionOfFirstNode]
        b2: Node = rt2.sequenceOfNodes[sm.positionOfSecondNode]
        x1 = b1.positionInRoute
        x2 = b2.positionInRoute
        b1.positionInRoute = x1
        b2.positionInRoute = x2
        rt1.sequenceOfNodes[sm.positionOfFirstNode] = b2
        rt2.sequenceOfNodes[sm.positionOfSecondNode] = b1
        rt1.sequenceOfNodes[sm.positionOfFirstNode].positionInRoute = sm.positionOfFirstNode
        rt2.sequenceOfNodes[sm.positionOfSecondNode].positionInRoute = sm.positionOfSecondNode
        if (rt1 == rt2):
            rt1.cost += sm.moveCost
        else:
            rt1.cost += sm.costChangeFirstRt
            rt2.cost += sm.costChangeSecondRt
            rt1.load = rt1.load - b1.demand + b2.demand
            rt2.load = rt2.load + b1.demand - b2.demand

        self.sol.cost += sm.moveCost

        newCost = self.CalculateTotalCost(self.sol)
        # debuggingOnly
        if abs((newCost - oldCost) - sm.moveCost) > 0.0001:
            print(newCost - oldCost)
            print(sm.moveCost)
            print('Cost Issue')

    def ReportSolution(self, sol):
        tc = 0
        for i in range(0, len(sol.routes)):
            rt = sol.routes[i]
            for j in range(0, len(rt.sequenceOfNodes)):
                if j == 0:
                    print(rt.sequenceOfNodes[j].ID, end=' ')
                else:
                    print("," + str(rt.sequenceOfNodes[j].ID), end=' ')
            print(rt.cost)
            tc += rt.cost
            sol.cost += rt.cost
        print(tc)

    def StoreBestRelocationMove(self, originRouteIndex, targetRouteIndex, originNodeIndex, targetNodeIndex, moveCost,
                                originRtCostChange, targetRtCostChange, rm: RelocationMove):
        rm.originRoutePosition = originRouteIndex
        rm.originNodePosition = originNodeIndex
        rm.targetRoutePosition = targetRouteIndex
        rm.targetNodePosition = targetNodeIndex
        rm.costChangeOriginRt = originRtCostChange
        rm.costChangeTargetRt = targetRtCostChange
        rm.moveCost = moveCost

    def StoreBestSwapMove(self, firstRouteIndex, secondRouteIndex, firstNodeIndex, secondNodeIndex, moveCost,
                          costChangeFirstRoute, costChangeSecondRoute, sm):
        sm.positionOfFirstRoute = firstRouteIndex
        sm.positionOfSecondRoute = secondRouteIndex
        sm.positionOfFirstNode = firstNodeIndex
        sm.positionOfSecondNode = secondNodeIndex
        sm.costChangeFirstRt = costChangeFirstRoute
        sm.costChangeSecondRt = costChangeSecondRoute
        sm.moveCost = moveCost

    def CalculateTotalCost(self, sol):
        c = 0
        for i in range(0, len(sol.routes)):
            rt = sol.routes[i]
            ci, l = self.calculate_route_details(rt.sequenceOfNodes)
            c += ci
        return c

    def InitializeOperators(self, rm, sm, top):
        rm.Initialize()
        sm.Initialize()
        top.Initialize()

    def FindBestTwoOptMove(self, top):
        for rtInd1 in range(0, len(self.sol.routes)):
            rt1: Route = self.sol.routes[rtInd1]
            for rtInd2 in range(rtInd1, len(self.sol.routes)):
                rt2: Route = self.sol.routes[rtInd2]
                for nodeInd1 in range(0, len(rt1.sequenceOfNodes) - 1):
                    start2 = 0
                    if (rt1 == rt2):
                        start2 = nodeInd1 + 2
                    for nodeInd2 in range(start2, len(rt2.sequenceOfNodes) - 1):
                        moveCost = 10 ** 9

                        A = rt1.sequenceOfNodes[nodeInd1]
                        B = rt1.sequenceOfNodes[nodeInd1 + 1]
                        K = rt2.sequenceOfNodes[nodeInd2]
                        L = rt2.sequenceOfNodes[nodeInd2 + 1]

                        if rt1 == rt2:
                            if nodeInd1 == 0 and nodeInd2 == len(rt1.sequenceOfNodes) - 2:
                                continue
                            cost_now = rt1.cost
                            testrt = copy.deepcopy(rt1)
                            reversedSegment = reversed(
                                testrt.sequenceOfNodes[nodeInd1 + 1: nodeInd2 + 1])
                            # lst = list(reversedSegment)
                            # lst2 = list(reversedSegment)
                            testrt.sequenceOfNodes[nodeInd1 + 1: nodeInd2 + 1] = reversedSegment
                            cost_new = self.calculate_cost(testrt)
                            moveCost = cost_new - cost_now
                        else:
                            if nodeInd1 == 0 and nodeInd2 == 0:
                                continue
                            if nodeInd1 == len(rt1.sequenceOfNodes) - 2 and nodeInd2 == len(rt2.sequenceOfNodes) - 2:
                                continue

                            if self.CapacityIsViolated(rt1, nodeInd1, rt2, nodeInd2):
                                continue
                            cost_old = rt1.cost + rt2.cost
                            test1 = copy.deepcopy(rt1)
                            test2 = copy.deepcopy(rt2)
                            relocatedSegmentOfRt1 = test1.sequenceOfNodes[nodeInd1 + 1:]

                            # slice with the nodes from position top.positionOfFirstNode + 1 onwards
                            relocatedSegmentOfRt2 = test2.sequenceOfNodes[nodeInd2 + 1:]

                            del test1.sequenceOfNodes[nodeInd1 + 1:]
                            del test2.sequenceOfNodes[nodeInd2 + 1:]

                            test1.sequenceOfNodes.extend(relocatedSegmentOfRt2)
                            test2.sequenceOfNodes.extend(relocatedSegmentOfRt1)
                            # cost added not calculated this way
                            c1 = self.calculate_cost(test1)
                            c2 = self.calculate_cost(test2)
                            cost_new = c1 + c2
                            moveCost = cost_new - cost_old
                        if moveCost < top.moveCost:
                            self.StoreBestTwoOptMove(rtInd1, rtInd2, nodeInd1, nodeInd2, moveCost, top)

    def CapacityIsViolated(self, rt1, nodeInd1, rt2, nodeInd2):

        rt1FirstSegmentLoad = 0
        for i in range(0, nodeInd1 + 1):
            n = rt1.sequenceOfNodes[i]
            rt1FirstSegmentLoad += n.demand
        rt1SecondSegmentLoad = rt1.load - rt1FirstSegmentLoad

        rt2FirstSegmentLoad = 0
        for i in range(0, nodeInd2 + 1):
            n = rt2.sequenceOfNodes[i]
            rt2FirstSegmentLoad += n.demand
        rt2SecondSegmentLoad = rt2.load - rt2FirstSegmentLoad

        if (rt1FirstSegmentLoad + rt2SecondSegmentLoad > rt1.capacity):
            return True
        if (rt2FirstSegmentLoad + rt1SecondSegmentLoad > rt2.capacity):
            return True

        return False

    def StoreBestTwoOptMove(self, rtInd1, rtInd2, nodeInd1, nodeInd2, moveCost, top):
        top.positionOfFirstRoute = rtInd1
        top.positionOfSecondRoute = rtInd2
        top.positionOfFirstNode = nodeInd1
        top.positionOfSecondNode = nodeInd2
        top.moveCost = moveCost

    def calculate_route_details(self, nodes_sequence):
        rt_load = 0
        rt_cumulative_cost = 0
        tot_time = 0
        for i in range(len(nodes_sequence) - 1):
            from_node = nodes_sequence[i]
            to_node = nodes_sequence[i + 1]
            tot_time += self.distanceMatrix[from_node.ID][to_node.ID]
            rt_cumulative_cost += tot_time
            tot_time += 10
            rt_load += from_node.demand
        return rt_cumulative_cost, rt_load

    def ApplyTwoOptMove(self, top):
        rt1: Route = self.sol.routes[top.positionOfFirstRoute]
        rt2: Route = self.sol.routes[top.positionOfSecondRoute]

        if rt1 == rt2:
            # reverses the nodes in the segment [positionOfFirstNode + 1,  top.positionOfSecondNode]
            reversedSegment = reversed(rt1.sequenceOfNodes[top.positionOfFirstNode + 1: top.positionOfSecondNode + 1])
            # lst = list(reversedSegment)
            # lst2 = list(reversedSegment)
            rt1.sequenceOfNodes[top.positionOfFirstNode + 1: top.positionOfSecondNode + 1] = reversedSegment

            # reversedSegmentList = list(reversed(rt1.sequenceOfNodes[top.positionOfFirstNode + 1: top.positionOfSecondNode + 1]))
            # rt1.sequenceOfNodes[top.positionOfFirstNode + 1: top.positionOfSecondNode + 1] = reversedSegmentList

            rt1.cost += top.moveCost

        else:
            # slice with the nodes from position top.positionOfFirstNode + 1 onwards
            relocatedSegmentOfRt1 = rt1.sequenceOfNodes[top.positionOfFirstNode + 1:]

            # slice with the nodes from position top.positionOfFirstNode + 1 onwards
            relocatedSegmentOfRt2 = rt2.sequenceOfNodes[top.positionOfSecondNode + 1:]

            del rt1.sequenceOfNodes[top.positionOfFirstNode + 1:]
            del rt2.sequenceOfNodes[top.positionOfSecondNode + 1:]

            rt1.sequenceOfNodes.extend(relocatedSegmentOfRt2)
            rt2.sequenceOfNodes.extend(relocatedSegmentOfRt1)

            self.UpdateRouteCostAndLoad(rt1)
            self.UpdateRouteCostAndLoad(rt2)
            self.Update_Positions_In_Routes(rt1,rt2)

        self.sol.cost += top.moveCost

    def UpdateRouteCostAndLoad(self, rt: Route):
        c, li = self.calculate_route_details(rt.sequenceOfNodes)
        tc = c
        tl = li
        rt.load = li
        rt.cost = c

    def calculate_cost(self, rt: Route):
        cost = 0
        for i in range(0, len(rt.sequenceOfNodes) - 1):
            n = len(rt.sequenceOfNodes) - i - 1
            cost += n * self.distanceMatrix[rt.sequenceOfNodes[i].ID][rt.sequenceOfNodes[i + 1].ID]
        n = len(rt.sequenceOfNodes) - 2
        cummulative = n * (n + 1) / 2
        cost += cummulative * 10
        return cost

    def IdentifyMinimumCostInsertion(self, best_insertion):
        for i in range(0, len(self.customers)):
            candidateCust: Node = self.customers[i]
            if candidateCust.isRouted is False:
                for rt in self.sol.routes:
                    if rt.load + candidateCust.demand <= rt.capacity:
                        for j in range(0, len(rt.sequenceOfNodes) - 1):
                            A = rt.sequenceOfNodes[j]
                            B = rt.sequenceOfNodes[j + 1]
                            costAdded = self.distanceMatrix[A.ID][candidateCust.ID] + \
                                        self.distanceMatrix[candidateCust.ID][
                                            B.ID]
                            costRemoved = self.distanceMatrix[A.ID][B.ID]
                            trialCost = costAdded - costRemoved
                            if trialCost < best_insertion.cost:
                                best_insertion.customer = candidateCust
                                best_insertion.route = rt
                                best_insertion.insertionPosition = j
                                best_insertion.cost = trialCost
                    else:
                        continue

    def ApplyCustomerInsertionAllPositions(self, insertion):
        insCustomer = insertion.customer
        rt = insertion.route
        # before the second depot occurrence
        insIndex = insertion.insertionPosition
        rt.sequenceOfNodes.insert(insIndex + 1, insCustomer)
        rt.cost += insertion.cost
        self.sol.cost += insertion.cost
        rt.load += insCustomer.demand
        insCustomer.isRouted = True
