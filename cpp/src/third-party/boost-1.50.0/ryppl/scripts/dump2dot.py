# Copyright (C) 2012 Dave Abrahams <dave@boostpro.com>
#
# Distributed under the Boost Software License, Version 1.0.  See
# accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt

import sys_path_setup
import glob, sys
from subprocess import check_output, check_call, Popen, PIPE
from xml.etree.cElementTree import ElementTree, Element
from path import Path
from read_dumps import read_dumps
from depgraph import *
from transitive import *
from display_graph import *

colors=('red','green','orange', 'blue', 'indigo', 'violet')

def first(seq):
    return iter(seq).next()

def run(dump_dir=None):
    g = dumps = read_dumps(dump_dir)

    # find all strongly-connected components
    from SCC import SCC
    successors = Successors(g)
    sccs = SCC(str, lambda i: successors(g, i)).getsccs(g)
    # map each vertex to a scc set
    scc = {}
    for component in sccs:
        s = set(component)
        for u in s:
            scc[u] = s

    long_sccs = [s for s in sccs if len(s) > 1]
    print 'long_sccs=', long_sccs

    # color each vertex in an SCC of size > 1 according to its SCC
    color = {}
    for i,s in enumerate(long_sccs):
        for u in s:
            color[u] = colors[i]

    if '--all' in sys.argv[1:]:
        V = g
    else:
        V = set(u for u in g if successors(g,u))

    direct_graph = to_mutable_graph(g, direct_successors, V.__contains__)

    t_redux = to_mutable_graph(g, UsageSuccessors(g), vertex_filter=V.__contains__)
    inplace_transitive_reduction(t_redux)

    class Format(object):
        def vertex_attributes(self, s):
            ret = ['color='+color[s]] if s in color else []
            ret += ['fontsize=9', 'fontname=GillSans']
            if dumps[s].find('libraries/library') is not None:
                ret+=['shape=box3d','style=bold']
            if s.startswith('Boost') and len(s) > len('Boost'):
                ret +=['label=".'+s[len('Boost'):]+'"']
            return ret

        def edge_attributes(self, s, t):
            if t in direct_graph[s]:
                return ['style=bold']
            elif t in t_redux[s]:
                return ['style=dashed','color=blue']
            else:
                return ['style=dashed','color=gray']

    if '--all' in sys.argv[1:]:
        # full_graph = dict((v, direct_graph[v] | t_redux[v]) for v in g)
        full_graph = to_mutable_graph(g, vertex_filter=V.__contains__)
        show_digraph(full_graph, formatter=Format(), ranksep=0.5, nodesep=0.05, splines='true', layout='dot')
        # show_digraph(full_graph, formatter=Format(), layout='neato', overlap='false', ordering='out', nodesep=0.05, splines='true')
    else:
        full_graph = to_mutable_graph(g, vertex_filter=V.__contains__)
        show_digraph(full_graph, formatter=Format(), ranksep=0.5, nodesep=0.05, splines='true', layout='dot')

if __name__ == '__main__':
    argv = set(sys.argv[1:])
    argv.discard('--all')
    run(dump_dir=Path(argv.pop()) if len(argv) > 0 else None)
