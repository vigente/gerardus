import tempfile
import subprocess

class Formatter:
    def vertex_attributes(self, s):
        return None
    def edge_attributes(self, s, t):
        return None

def digraph(g, formatter=Formatter(), **kw):
    result = ['digraph G {']
    result += ['    %s=%s;' % kv for kv in kw.items()]

    for s in g:
        line = '    %s' % s
        a = formatter.vertex_attributes(s)
        if a and len(a):
            line += ' [%s]' % ','.join(a)
        result.append(line)

        for t in g[s]:
            line = '        %s -> %s' % (s,t)
            a = formatter.edge_attributes(s,t)
            if a and len(a):
                line += ' [%s]' % ','.join(a)
            result.append(line)

    result.append('}')
    return '\n'.join(result)

def show_graphviz(text, format='pdf'):
    graph = tempfile.NamedTemporaryFile(suffix='.gv',delete=False)
    graph.write(text)
    graph.flush()

    result = tempfile.NamedTemporaryFile(suffix='.'+format, delete=False)
    result.write(
        subprocess.check_output(['dot','-T'+format,graph.name]))
    result.flush()
    subprocess.check_call(['open',result.name])
    
def show_digraph(g, **kw):
    show_graphviz(digraph(g, **kw))
    
def show_digraph2(g0, g1, colors=['red','green','black'], **kw):
    class Format(Formatter):
        def edge_attributes(self,s,t):
            index = (1 if t in g0.get(s,[]) else 0) + (2 if t in g1.get(s,[]) else 0)
            return ['color='+colors[index-1]]
    g3 = dict(
        (s, set(g0.get(s,[]))|set(g1.get(s,[])) )
        for s in set(g0)|set(g1))

    show_digraph(g3, formatter=Format(), **kw)


if __name__ == '__main__':

    def random_graph(max, a, b):
        g = {}
        for x in range(1,max):
            g[x] = set(y for y in range(1,max) if (x*x+y) % a == b)
        return g
    
    g0 = random_graph(15,11,5)
    import pprint
    pprint.pprint(g0)
    show_digraph(g0,layout='dot')
    g1 = random_graph(11,7,5)
    show_digraph(g1)
    show_digraph2(g0,g1)
