from xml.etree.cElementTree import Element
from lazydict import lazydict
from transitive import *

def direct_successors(dumps, v):
    ret = set(fp.findtext('arg') for fp in dumps[v].findall('find-package'))
    ret.discard(v)
    return ret

class UsageSuccessors(object):
    def __init__(self, dumps):
        # Cache the transitive closure of all usage dependencies
        self.usage_dependencies = to_mutable_graph(
            dumps, 
            lambda dumps,v: 
            set(d.text for d in dumps[v].findall('depends/dependency') ))
        inplace_transitive_closure(self.usage_dependencies)

    def __call__(self, dumps, v):
        return set(d for s in direct_successors(dumps,v)
                   for d in self.usage_dependencies[s])

class Successors(object):
    def __init__(self, dumps):
        self.usage_successors = UsageSuccessors(dumps)
    def __call__(self, dumps, v):
        return direct_successors(dumps,v) | self.usage_successors(dumps, v)

def to_mutable_graph(dumps, successor_function=None, vertex_filter=lambda x:True):
    successor_function = successor_function or Successors(dumps)

    return lazydict(
        set 
      , ((v, set(x for x in successor_function(dumps,v) 
                 if vertex_filter(x) 
                 and x != v    # filter out self-loops
                 ))
         for v in dumps if vertex_filter(v)))

def run(dump_dir=None):
    from read_dumps import read_dumps
    from display_graph import show_digraph, show_digraph2
    from transitive import inplace_transitive_reduction
    dumps = read_dumps(dump_dir)
    
    direct = to_mutable_graph(dumps, direct_successors)
    usage = to_mutable_graph(dumps, usage_successors)

    inplace_transitive_reduction(direct)
    inplace_transitive_reduction(usage)
    
    show_digraph2(direct, usage, layout='neato', overlap='false', splines='True')

    # from pprint import pprint
    # pprint(direct)

    from SCC import SCC
    sccs = SCC(str, lambda i: successors(dumps, i)).getsccs(dumps)
    long_sccs = [s for s in sccs if len(s) > 1]
    assert len(long_sccs) == 0, str(long_sccs)

if __name__ == '__main__':
    import sys
    run(None if len(sys.argv) == 1 else sys.argv[1])
