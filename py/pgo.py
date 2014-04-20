from goatools import obo_parser
from collections import defaultdict

def getGOParents(goGraph, term):
    acc = {}
    parents = goGraph[term].parents
    if not parents:
        return acc
    else:
        for parent in parents:
            acc[parent.id] = True
            for parent2 in getGOParents(goGraph, parent.id):
                acc[parent2] = True
    return acc

if __name__ == "__main__":
    pass
