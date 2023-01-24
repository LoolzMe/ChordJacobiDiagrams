from graph_tool.all import *
import json

# g = Graph(directed=False)
# v = []

# v = list(g.add_vertex(4))

# g.add_edge(v[0], v[1])
# g.add_edge(v[0], v[1])
# g.add_edge(v[0], v[2])
# g.add_edge(v[1], v[2])
# g.add_edge(v[2], v[3])

# graph_draw(g, pos=sfdp_layout(g), vertex_text=g.vertex_index)

filenames = ["jacobiDiagrams.json", "notReducedGraphs.json", "IHXchanged.json"]


with open(filenames[1]) as d_file:
    data = json.load(d_file)

    for graph in data:
        vp = graph['vp']
        ep = graph['ep']
        vl = graph['vl']
        nv = graph['nv']

        g = Graph(directed=False)
        v = g.add_vertex(nv)

        for i in range(len(ep) // 2):
            d1 = 2*i
            d2 = 2*i + 1
            print(d1, d2)

            g.add_edge(vl[d1], vl[d2])

        graph_draw(g, pos=sfdp_layout(g), vertex_text=g.vertex_index, vertex_size=30) 