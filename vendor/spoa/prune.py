'''
This script is used to prune a messy variation graph generated from noisy long reads directly,
which is based on some statistical tricks.
'''
import copy
import sys
import random
import networkx as nx
from scipy.stats import chi2_contingency


def is_connected(gfa_file):
    links = []
    with open(gfa_file) as  fr:
        for line in fr:
            if line.startswith('S'):
                pass
            elif line.startswith('L'):
                a = line.split()
                links.append((a[1], a[3]))
    G = nx.Graph()
    G.add_edges_from(links)
    return nx.is_connected(G)


def read_gfa(gfa_file):
    nodes = {}  # node to node record
    links = {}  # link to link record
    paths = []
    path_nodes = {}  # nodes involved in all paths
    path_links2 = {}  # links to weight, double nodes
    path_links3 = {}  # links to weight, triple nodes
    flag = True
    with open(gfa_file) as  fr:
        for line in fr:
            if line.startswith('S'):
                nodes[line.split()[1]] = line
            elif line.startswith('L'):
                a = line.split()
                links[a[1] + a[2] + ':' + a[3] + a[4]] = line
            elif line.startswith('P'):
                a = line.split()
                paths.append(line)
                pre_node = ''
                for i, node in enumerate(a[2].split(',')):
                    # node = node.replace('+','').replace('-','')  # ignore +/-
                    if node[-1] == '-':
                        raise Exception("Check path direction, only + is allowed.")
                    if node in path_nodes:
                        path_nodes[node] += 1
                    else:
                        path_nodes[node] = 1
                    if i == 0:
                        pre_node = node
                        pre_pre_node = node
                    else:
                        link = pre_node + ':' + node

                        if i >= 2:
                            link3 = pre_pre_node + ':' + pre_node + ':' + node
                            if link3 in path_links3:
                                path_links3[link3] += 1
                            else:
                                path_links3[link3] = 1
                            pre_pre_node = pre_node

                        if link in path_links2:
                            path_links2[link] += 1
                        else:
                            path_links2[link] = 1
                        pre_node = node
            else:
                continue
    return nodes, links, path_nodes, path_links2, path_links3, paths


def prune_edges(in_links, path_links, min_weight, min_confidence, min_support):
    '''prune edges by three metrics'''
    links = copy.deepcopy(in_links)
    del_links = {}
    G = nx.DiGraph()
    G.add_edges_from([link.split(':') for link in in_links])
    print('Graph is connected: {}'.format(nx.is_connected(G.to_undirected())))
    edge2weight = {(link.split(':')[0], link.split(':')[1]): weight for link, weight in path_links.items()}
    v2cov = {}  # vertex to coverage, for computing support
    if min_support:  # if min_support==0, skip support computation
        v2layer = longest_path_layering_vg(G)
        print('layering graph has been finished...')
        v2cov = get_vertex_cov(v2layer, edge2weight)  # the overall sequencing coverage at this position

    for edge in G.edges():
        node1, node2 = edge  # use the source node (source --> target)
        weight = 0
        if edge in edge2weight:
            weight = edge2weight[edge]
        else:
            del_links[node1 + ':' + node2] = 1
            continue

        support = 1
        if min_support:
            support = weight * 1.0 / v2cov[node1]
        confidence_out = 0  # confidence(X—>Y),out-edges of X
        confidence_in = 0  # confidence(Y—>X),in-edges of Y

        # confidence(X—>Y),out-edges of X
        cov_hash = {}
        for n in G.successors(node1):
            if (node1 + ':' + n) in path_links:
                cov_hash[n] = path_links[node1 + ':' + n]
            else:
                cov_hash[n] = 0
        cov_sum = sum(cov_hash.values())
        if cov_sum > 0:
            confidence_out = weight * 1.0 / cov_sum

        # confidence(Y—>X),in-edges of Y
        cov_hash = {}
        for n in G.predecessors(node2):
            if (n + ':' + node2) in path_links:
                cov_hash[n] = path_links[n + ':' + node2]
            else:
                cov_hash[n] = 0
        cov_sum = sum(cov_hash.values())
        if cov_sum > 0:
            confidence_in = weight * 1.0 / cov_sum

        # remove edges via min_weight, min_confidence, min_support
        # if weight >= min_weight and confidence_out >= min_confidence and confidence_in >= min_confidence:
        if weight >= min_weight and support >= min_support and confidence_out >= min_confidence and confidence_in >= min_confidence:
            continue
        else:
            del_links[node1 + ':' + node2] = 1

    for link in del_links.keys():
        del links[link]
    return links


def prune_edges_test(in_links, path_links, min_weight, min_confidence, min_support,
                     p_cutoff, weight_cutoff):
    '''prune edges by three metrics, also use chisquare test'''
    links = copy.deepcopy(in_links)
    del_links = {}
    G = nx.DiGraph()
    G.add_edges_from([link.split(':') for link in in_links])
    print('Graph is connected: {}'.format(nx.is_connected(G.to_undirected())))
    edge2weight = {(link.split(':')[0], link.split(':')[1]): weight for link, weight in path_links.items()}
    v2cov = {}  # vertex to coverage, for computing support
    if min_support:  # if min_support==0, skip support computation
        v2layer = longest_path_layering_vg(G)
        print('layering graph has been finished...')
        v2cov = get_vertex_cov(v2layer, edge2weight)  # the overall sequencing coverage at this position

    for edge in G.edges():
        node1, node2 = edge  # use the source node (source --> target)
        weight = 0
        if edge in edge2weight:
            weight = edge2weight[edge]
        else:
            del_links[node1 + ':' + node2] = 1
            continue

        # w_XY   w_NX_Y
        # w_X_NY  w_NX_NY
        w_XY = weight

        support = 1
        if min_support:
            support = weight * 1.0 / v2cov[node1]
        confidence_out = 0  # confidence(X—>Y),out-edges of X
        confidence_in = 0  # confidence(Y—>X),in-edges of Y

        # confidence(X—>Y),out-edges of X
        cov_hash = {}
        for n in G.successors(node1):
            if (node1 + ':' + n) in path_links:
                cov_hash[n] = path_links[node1 + ':' + n]
            else:
                cov_hash[n] = 0
        cov_sum = sum(cov_hash.values())
        if cov_sum > 0:
            confidence_out = weight * 1.0 / cov_sum
        w_X_NY = cov_sum - w_XY

        # confidence(Y—>X),in-edges of Y
        cov_hash = {}
        for n in G.predecessors(node2):
            if (n + ':' + node2) in path_links:
                cov_hash[n] = path_links[n + ':' + node2]
            else:
                cov_hash[n] = 0
        cov_sum = sum(cov_hash.values())
        if cov_sum > 0:
            confidence_in = weight * 1.0 / cov_sum
        w_NX_Y = cov_sum - w_XY
        w_NX_NY = v2cov[node1] - w_XY - w_X_NY - w_NX_Y
        assert (w_X_NY >= 0)
        assert (w_XY >= 0)
        assert (w_NX_Y >= 0)
        # assert(w_NX_NY>=0)
        # print('w_NX_NY:{}'.format(w_NX_NY))
        w_NX_NY = max(1, w_NX_NY)

        # statistical test
        g, pval, dof, expctd = chi2_contingency([[w_XY, w_NX_Y], [w_X_NY, w_NX_NY]], correction=True)
        # print('pval:{}'.format(pval))

        # remove edges via min_weight, min_confidence, min_support
        if weight >= min_weight and support >= min_support and confidence_out >= min_confidence and confidence_in >= min_confidence:
            if pval > p_cutoff and weight <= weight_cutoff:
                del_links[node1 + ':' + node2] = 1
            else:
                continue
            continue
        else:
            del_links[node1 + ':' + node2] = 1

    for link in del_links.keys():
        del links[link]
    return links


def prune_nodes(in_nodes, path_nodes, min_cov):
    nodes = copy.deepcopy(in_nodes)
    for node in in_nodes.keys():
        if node not in path_nodes:
            del nodes[node]
        elif path_nodes[node] < min_cov:
            del nodes[node]
    return nodes


def write_gfa(nodes, links, paths, gfa_file):
    '''output gfa file according to links'''
    G = nx.Graph()
    G.add_edges_from([(x.split(':')[0][:-1], x.split(':')[1][:-1]) for x in links.keys()])
    print('The size of connected components is:{}'.
          format([len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]))
    G1 = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
    print('The largest component is connected:{}, size:{}'.format(nx.is_connected(G1), len(G1)))

    G1_nodes = {v for v in G1.nodes()}  # only keep nodes involved in links
    G1_edges = {(e1, e2) for e1, e2 in G1.edges()}
    with open(gfa_file, 'w') as fw:
        fw.write('H\tVN:Z:1.0\n')
        fw.write(''.join([v for k, v in nodes.items() if k in G1_nodes]))
        for link, link_info in links.items():
            node1, node2 = [a[:-1] for a in link.split(':')]
            if ((node1, node2) in G1.edges) or ((node2, node1) in G1_edges):
                fw.write(link_info)
        # write paths
        # fw.write(''.join(paths))
    return


def longest_path_layering_vg(G):
    """
    Modified version, See page 420-1 of Chapter 13: Hierarchical Graph Drawing by Healy and Nikolov.
    Book: http://cs.brown.edu/people/rtamassi/gdhandbook/chapters/hierarchical.pdf

    U = vertices assigned to any layer
    V = Vertex Set
    v is selected from v in V \ U
    """
    V = G.nodes()
    E = G.edges()
    VS = {v: {'s': set(), 'p': set()} for v in
          V}  # s is short for successors (i.e. v --> s) and p for predecessors (p --> s)
    for s, t in E:  # for source target
        VS[t]['p'].add(s)  # target <-- source
        VS[s]['s'].add(t)  # source --> target
    vertex_set = VS

    longest_path = nx.algorithms.dag.dag_longest_path(G)
    print('longest path computation is finished.')
    v2layer = {v: (i + 1) for i, v in enumerate(longest_path)}
    # print(v2layer)
    U = set([v for v in longest_path])
    V = set(vertex_set.keys())
    while U != V:
        for v in V - U:
            # print('U:{}, V:{}, choice:{}'.format(len(U), len(V), len(V-U)))
            sv_layers = []
            for sv in vertex_set[v]['s']:
                if sv in v2layer:
                    sv_layers.append(v2layer[sv])
                else:
                    continue
            if sv_layers:
                v2layer[v] = min(sv_layers) - 1
                U.add(v)
            else:
                pv_layers = []
                for pv in vertex_set[v]['p']:
                    if pv in v2layer:
                        pv_layers.append(v2layer[pv])
                    else:
                        continue
                if pv_layers:
                    v2layer[v] = max(pv_layers) + 1
                    U.add(v)
    return v2layer


def get_vertex_cov(v2layer, edge2weight):
    '''calculate the overall sequencing coverage at the position of vertex v,
    which is equal to the support of v: S(v)
    Note that it is not the number of reads passing through v.
    '''
    v2weight = {}  # weight means the overall sequencing coverage at this position

    layer2vertexes = {}
    for v, layer in v2layer.items():
        if layer in layer2vertexes:
            layer2vertexes[layer].append(v)
        else:
            layer2vertexes[layer] = [v]
    for e in edge2weight.keys():
        v1, v2 = e
        # e = (v1, v2)
        layer1 = v2layer[v1]
        layer2 = v2layer[v2]
        for layer in range(layer1, layer2):
            for v in layer2vertexes[layer]:
                if v in v2weight:
                    v2weight[v] += edge2weight[e]
                else:
                    v2weight[v] = edge2weight[e]
    # there is no need to commpute the weight of the sink node
    # vertex at the same layer have the same weight
    len(v2weight)
    return v2weight


def prune_edges_by_triple_nodes(in_links, path_links, min_weight):
    '''prune edges which are insufficiently supported by walk of triple nodes(A->B->C)'''
    min_weight = int(0.8 *min_weight)
    links = copy.deepcopy(in_links)
    link2weight = {link: 0 for link in links}  # support by triple nodes
    del_links = {}
    G = nx.DiGraph()
    G.add_edges_from([link.split(':') for link in in_links])
    print('Graph is connected: {}'.format(nx.is_connected(G.to_undirected())))
    for node in G.nodes():
        if len(list(G.predecessors(node))) == 0 or len(list(G.successors(node))) == 0:
            continue
        else:
            for pnode in G.predecessors(node):
                for snode in G.successors(node):
                    triple_node = pnode + ':' + node + ':' + snode
                    if triple_node in path_links and path_links[triple_node] >= min_weight:
                        link2weight[pnode + ':' + node] += 1
                        link2weight[node + ':' + snode] += 1

    with open("link2weight_triplenodes.txt", 'w') as fw:
        fw.write('\n'.join([k + '\t' + str(v) for k, v in link2weight.items()]))

    i=0
    for link, weight in link2weight.items():
        if weight < 1:
            del links[link]
            i+=1
    print("## Removing {} edges based on triple-nodes stratergy.".format(i))

    return links


if __name__ == '__main__':
    # gfa_file = './data/test.china100x.gfa' # must have path information (like output of vg)
    # gfa_file = './data/test.china100x.ovlpbased.gfa'
    # gfa_file = './data/test.china100x.acc100.gfa'

    # node IDs in the gfa_file should be already sorted and relabeled from 1..N,
    # because the position of node in the genome will be determined through node id.
    gfa_file, outfile, min_weight, min_confidence, min_support, p_cutoff, weight_cutoff, chi2 = sys.argv[1:]
    chi2 = True if chi2.lower() == 'true' else False

    nodes, links, path_nodes, path_links2, path_links3, paths = read_gfa(gfa_file)

    with open('path_links.out', 'w') as fw:
        fw.write('\n'.join([k.replace(':', '\t').replace('+', '') + '\t' + str(v)
                            for k, v in sorted(path_links2.items(), key=lambda h: h[1])]))

    # print(path_nodes)
    print('# original nodes: {}'.format(len(nodes)))
    # nodes = prune_nodes(nodes, path_nodes, int(min_cov))
    if not chi2:
        # print('Do not use chi-square test\n'*10)
        links = prune_edges(links, path_links2, int(min_weight), float(min_confidence), float(min_support))
        links = prune_edges_by_triple_nodes(links, path_links3, int(min_weight))
    else:
        links = prune_edges_test(links, path_links2, int(min_weight), float(min_confidence),
                                 float(min_support), float(p_cutoff), int(weight_cutoff))
    print('# pruned nodes: {}'.format(len(nodes)))

    # write_gfa(nodes, links, './data/out.gfa')
    write_gfa(nodes, links, paths, outfile)
    print("The pruned graph is connected : {}".format(is_connected(outfile)))
