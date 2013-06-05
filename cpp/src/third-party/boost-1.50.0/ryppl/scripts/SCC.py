# Shamelessly lifted from Tim Peters' email 
# http://mail.python.org/pipermail/python-dev/2003-April/034656.html

# This implements Tarjan's linear-time algorithm for finding the maximal
# strongly connected components.  It takes time proportional to the sum
# of the number of nodes and arcs.
#
# Two functions must be passed to the constructor:
#     node2id     graph node -> a unique integer
#     successors  graph node -> sequence of immediate successor graph nodes
#
# Call method getsccs() with an iterable producing the root nodes of the graph.
# The result is a list of SCCs, each of which is a list of graph nodes.
# This is a partitioning of all graph nodes reachable from the roots,
# where each SCC is a maximal subset such that each node in an SCC is
# reachable from all other nodes in the SCC.  Note that the derived graph
# where each SCC is a single "supernode" is necessarily acyclic (else if
# SCC1 and SCC2 were in a cycle, each node in SCC1 would be reachable from
# each node in SCC1 and SCC2, contradicting that SCC1 is a maximal subset).
# The list of SCCs returned by getsccs() is in a topological sort order wrt
# this derived DAG.

class SCC(object):
    def __init__(self, node2id, successors):
        self.node2id = node2id
        self.successors = successors

    def getsccs(self, roots):
        import sys

        node2id, successors = self.node2id, self.successors
        get_dfsnum = iter(xrange(sys.maxint)).next
        id2dfsnum = {}
        id2lowest = {}
        stack = []
        id2stacki = {}
        sccs = []

        def dfs(v, v_id):
            id2dfsnum[v_id] = id2lowest[v_id] = v_dfsnum = get_dfsnum()
            id2stacki[v_id] = len(stack)
            stack.append((v, v_id))
            for w in successors(v):
                w_id = node2id(w)
                if w_id not in id2dfsnum:   # first time we saw w
                    dfs(w, w_id)
                    id2lowest[v_id] = min(id2lowest[v_id], id2lowest[w_id])
                else:
                    w_dfsnum = id2dfsnum[w_id]
                    if w_dfsnum < v_dfsnum and w_id in id2stacki:
                        id2lowest[v_id] = min(id2lowest[v_id], w_dfsnum)

            if id2lowest[v_id] == v_dfsnum:
                i = id2stacki[v_id]
                scc = []
                for w, w_id in stack[i:]:
                    del id2stacki[w_id]
                    scc.append(w)
                del stack[i:]
                sccs.append(scc)

        for v in roots:
            v_id = node2id(v)
            if v_id not in id2dfsnum:
                dfs(v, v_id)
        sccs.reverse()
        return sccs

_basic_tests = """
>>> succs = {1: [2], 2: []}
>>> s = SCC(int, lambda i: succs[i])

The order in which the roots are listed doesn't matter:  we get the unique
topsort regardless.

>>> s.getsccs([1])
[[1], [2]]
>>> s.getsccs([1, 2])
[[1], [2]]
>>> s.getsccs([2, 1])
[[1], [2]]

But note that 1 isn't reachable from 2, so giving 2 as the only root won't
find 1.

>>> s.getsccs([2])
[[2]]

>>> succs = {1: [2],
...          2: [3, 5],
...          3: [2, 4],
...          4: [3],
...          5: [2]}
>>> s = SCC(int, lambda i: succs[i])
>>> s.getsccs([1])
[[1], [2, 3, 4, 5]]
>>> s.getsccs(range(1, 6))
[[1], [2, 3, 4, 5]]

Break the link from 4 back to 2.
>>> succs[4] = []
>>> s.getsccs([1])
[[1], [2, 3, 5], [4]]

You can do the same thing with different vertex types

>>> succs = {'a': ['b'],
...          'b': ['c', 'e'],
...          'c': ['b', 'd'],
...          'd': ['c'],
...          'e': ['b']}
>>> s = SCC(lambda x:x, lambda i: succs[i])
>>> s.getsccs('abcde')
[['a'], ['b', 'c', 'd', 'e']]

Break the link from d back to b.
>>> succs['d'] = []
>>> s.getsccs('a')
[['a'], ['b', 'c', 'e'], ['d']]
"""

__test__ = {'basic': _basic_tests}

def _test():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _test()

def trygc():
    import gc
    gc.collect()
    s = SCC(id, gc.get_referents)
    for scc in s.getsccs(gc.get_objects()):
        if len(scc) == 1:
            continue
        print "SCC w/", len(scc), "objects"
        for x in scc:
            print "   ", hex(id(x)), type(x),
            if hasattr(x, "__name__"):
                print x.__name__,
            print
