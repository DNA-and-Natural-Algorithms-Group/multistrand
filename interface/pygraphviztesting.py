####################################################################
#                                                                  #
#  Copyright (c) 2010-2015 California Institute of Technology.     #
#  Distributed under the MIT License.                              #
#  (See accompanying file LICENSE or copy at                       #
#  http://opensource.org/licenses/MIT)                             #
#                                                                  #
####################################################################
#                                                                  #
#   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)         #
#                                                                  #
####################################################################

import pygraphviz as pygv
#
# a = pygv.AGraph()
# a
# node_ids = range(30)
# random
# import random
# random.choice 
# help(random.choice)
# node_data = [(i,random.choice('AGCT')) for i in node_ids]
# node_data


def add_backbone( sequence, start_idx, graph, draw_labels = True, draw_3prime=False ):
    # nodes
    indices = range(start_idx, start_idx + len(sequence))
    A_bases = ["{0:04d}".format(idx) for idx,base in zip(indices,sequence) if base == 'A' or base == 'a']
    T_bases = ["{0:04d}".format(idx) for idx,base in zip(indices,sequence) if base == 'T' or base == 't']
    G_bases = ["{0:04d}".format(idx) for idx,base in zip(indices,sequence) if base == 'G' or base == 'g']
    C_bases = ["{0:04d}".format(idx) for idx,base in zip(indices,sequence) if base == 'C' or base == 'c']

    graph.node_attr.update(shape = 'circle', width = 0.3, fontsize=10, height = 0.3 )
    if draw_labels:
        A_dict = {'color':'red','label':'A'}
        T_dict = {'color':'purple','label':'T'}
        G_dict = {'color':'blue','label':'G'}
        C_dict = {'color':'green','label':'C'}
    else:
        A_dict = {'color':'red','label':' ','fillcolor':'red','style':'filled'}
        T_dict = {'color':'purple','label':' ','fillcolor':'purple','style':'filled'}
        G_dict = {'color':'blue','label':' ','fillcolor':'blue','style':'filled'}
        C_dict = {'color':'green','label':' ','fillcolor':'green','style':'filled'}

    nodes = []
    for n in A_bases:
        nodes.append((n,A_dict))
    for n in T_bases:
        nodes.append((n,T_dict))
    for n in G_bases:
        nodes.append((n,G_dict))
    for n in C_bases:
        nodes.append((n,C_dict))

    if draw_3prime:
        newnodename = "{0:04d}".format(start_idx+len(sequence)-1) + "_3'"
        nodes.append((newnodename,{'style':"invis",'label':" "}))
        
    nodes.sort( key=lambda x:x[0], reverse=True)

    for n,d in nodes:
        graph.add_node( n, **d )

    #edges
    backbone = [("{0:04d}".format(i),"{0:04d}".format(i+1)) for i in indices[:-1]]
    
    graph.add_edges_from(backbone, len=.5, weight=1.0, color='white')
    if draw_3prime:
        graph.add_edge("{0:04d}".format(start_idx+len(sequence)-1),newnodename,arrowhead="rvee",dir="forward",len=.6,headlabel="3'",labelangle="90",labeldistance="1.0")

def add_basepairs( bp_list, graph ):
    """ bp_list must be pairs of indices """
    graph.add_edges_from(bp_list, len=.7, weight=16.0,color='crimson')#,style='dotted')

def add_3prime_label( final_node_name, graph ):
    newnodename = str(final_node_name) + "_3'"
    graph.add_node(newnodename,style="invis",label=" ")
    graph.add_edge(final_node_name,newnodename,arrowhead="rvee",dir="forward",len=.6,headlabel="3'",labelangle="90",labeldistance="1.0")

# def get_coordinates_from_nodes( graph ):
#     nodes = graph.nodes()
#     nodes.sort()
#     pass

def backbone_elements( sequence, start_idx, draw_labels = True, draw_3prime=False ):
    # nodes
    indices = range(start_idx, start_idx + len(sequence))
    A_bases = ["{0:04d}".format(idx) for idx,base in zip(indices,sequence) if base == 'A' or base == 'a']
    T_bases = ["{0:04d}".format(idx) for idx,base in zip(indices,sequence) if base == 'T' or base == 't']
    G_bases = ["{0:04d}".format(idx) for idx,base in zip(indices,sequence) if base == 'G' or base == 'g']
    C_bases = ["{0:04d}".format(idx) for idx,base in zip(indices,sequence) if base == 'C' or base == 'c']

    if draw_labels:
        A_dict = {'color':'red','label':'A','fontcolor':'white','fillcolor':'red','style':'filled'}
        T_dict = {'color':'purple','label':'T','fontcolor':'white','fillcolor':'purple','style':'filled'}
        G_dict = {'color':'blue','label':'G','fontcolor':'white','fillcolor':'blue','style':'filled'}
        C_dict = {'color':'green','label':'C','fontcolor':'white','fillcolor':'green','style':'filled'}
    else:
        A_dict = {'color':'red','label':' ','fillcolor':'red','style':'filled'}
        T_dict = {'color':'purple','label':' ','fillcolor':'purple','style':'filled'}
        G_dict = {'color':'blue','label':' ','fillcolor':'blue','style':'filled'}
        C_dict = {'color':'green','label':' ','fillcolor':'green','style':'filled'}

    nodes = []
    for n in A_bases:
        nodes.append((n,A_dict))
    for n in T_bases:
        nodes.append((n,T_dict))
    for n in G_bases:
        nodes.append((n,G_dict))
    for n in C_bases:
        nodes.append((n,C_dict))

    if draw_3prime:
        newnodename = "{0:04d}".format(start_idx+len(sequence)-1) + "_3'"
        nodes.append((newnodename,{'style':"invis",'label':" "}))
        
    # nodes.sort( key=lambda x:x[0], reverse=True)

    # for n,d in nodes:
    #     graph.add_node( n, **d )

    #edges
    edges = []
    backbone = [("{0:04d}".format(i),"{0:04d}".format(i+1)) for i in indices[:-1]]

    for e in backbone:
        edges.append((e,{'len':.5,'weight':1.0,'color':'white'}))
        
    # graph.add_edges_from(backbone, len=.5, weight=1.0) 
    if draw_3prime:
        edges.append((("{0:04d}".format(start_idx+len(sequence)-1),newnodename),{'arrowhead':"lvee",'dir':"forward",'len':.6,'headlabel':"3'",'labelangle':"90",'labeldistance':"1.0",'color':'white','fontcolor':'white'}))
    return (nodes,edges)

    
def create_graph( input_seq, input_struc, start_idx ):
    g = pygv.AGraph(normalize=True,start='regular')
    initial_idx = start_idx
    g.node_attr.update(shape = 'circle', width = 0.3, fontsize=10, height = 0.3 )
    g.graph_attr['viewport'] = '400.0,400.0,1.0'
    g.graph_attr['bgcolor'] = 'transparent'
    nodes = []
    edges = []
    for seq in input_seq.split("+"):
        #seq   = seq
        new_nodes,new_edges = backbone_elements( seq, start_idx, draw_labels=True, draw_3prime=True)
        nodes = nodes + new_nodes
        edges = edges + new_edges
        start_idx += len(seq)

    nodes.sort( key=lambda x:x[0], reverse=True)
    for n,d in nodes:
        g.add_node( n, **d )
    for e,d in edges:
        g.add_edge( e, **d )

    struc = input_struc.replace("+","") # remove + strand breaks
    left = []
    basepairs = []
    for s,idx in zip( struc, range(len(struc))):
        if s == '(':
            left.append("{0:04d}".format(idx))
        if s == ')':
            basepairs.append( (left.pop(), "{0:04d}".format(idx)) )
    add_basepairs( basepairs, g )
    if not g.has_edge( "{0:04d}".format(initial_idx), "{0:04d}".format(len(struc) + initial_idx - 1) ):
        g.add_edge( "{0:04d}".format(initial_idx), "{0:04d}".format(len(struc) + initial_idx - 1), style="invis", len=1.0, weight = 1)
    g.layout(prog='neato')
    g2 = g.copy()
    bounds = g2.graph_attr['bb'].split(",")
    bounds = [float(x) for x in bounds]
    if bounds[2] > 400.0 or bounds[3] > 400.0:
        if bounds[2] > bounds[3]:
            scaling = 400.0 / bounds[2]
        else:
            scaling = 400.0 / bounds[3]
        g = g.copy()
        g.graph_attr['viewport'] = '400.0,400.0,{0:.2f}'.format(scaling)
        g.layout(prog='neato',args="-s")

    return g,basepairs

def update_graph( input_graph, old_basepairs, newstructure ):
    #g = pygv.AGraph(normalize=True,start='regular')
    g = input_graph.copy()
    g.graph_attr['start'] = 'random'
    g.graph_attr['viewport'] = '600.0,600.0,1'
    removed_basepairs = old_basepairs[:]
    
    #g.node_attr.update(shape = 'circle', width = 0.3, fontsize=10, height = 0.3 )

    struc = newstructure.replace("+","") # remove + strand breaks
    left = []
    basepairs = []
    new_basepairs = []
    for s,idx in zip( struc, range(len(struc))):
        if s == '(':
            left.append("{0:04d}".format(idx))
        if s == ')':
            bp = (left.pop(),"{0:04d}".format(idx))
            basepairs.append(bp)
            if bp not in old_basepairs or bp == ( "{0:04d}".format(0), "{0:04d}".format(len(struc)-1) ):
                new_basepairs.append( bp )
            else:
                removed_basepairs.remove( bp )

    for edge in removed_basepairs:
        g.remove_edge( edge )

    add_basepairs( new_basepairs, g )
    if ( "{0:04d}".format(0), "{0:04d}".format(len(struc)-1) ) not in basepairs:
        g.add_edge( "{0:04d}".format(0), "{0:04d}".format(len(struc)-1), style="invis", len=1.0, weight = 1)

    g.layout(prog='neato',args="-s")
    bounds = g.graph_attr['bb'].split(",")
    bounds = [float(x) for x in bounds]
    if bounds[2] > 600.0 or bounds[3] > 600.0:
        if bounds[2] > bounds[3]:
            scaling = 600.0 / bounds[2]
        else:
            scaling = 600.0 / bounds[3]
        g = g.copy()
        g.graph_attr['viewport'] = '600.0,600.0,{0:.2f}'.format(scaling)
        g.layout(prog='neato',args="-s")
    return g, basepairs

def intermediate_frame( old_graph, new_graph, increment ):
    old_nodes = old_graph.nodes()
    new_nodes = new_graph.nodes()
    intermediate_graph = new_graph.copy()
    intermediate_graph.graph_attr['start'] = 'random'
    for a,b in zip(old_nodes,new_nodes):
        coord_a = [float(item) for item in a.attr['pos'].split(',')]
        coord_b = [float(item) for item in b.attr['pos'].split(',')]
        diff = (coord_b[0] - coord_a[0], coord_b[1] - coord_a[1])
        new_coord = (coord_a[0] + diff[0] * increment, coord_a[1] + diff[1] * increment)
        n = intermediate_graph.get_node(str(a))
        n.attr['pos'] = "{0[0]},{0[1]}!".format( new_coord )
    intermediate_graph.layout(prog='neato',args='-s')
    return intermediate_graph

def create_movie( trajectory, intermediate_count=0 ):
    g,bp = create_graph( trajectory[0][0][3], trajectory[0][0][4], 0 )
    g.draw('test0000.png')
    counter = 1
    for state in trajectory[1:]:
        new_g, new_bp = update_graph( g, bp, state[0][4])
        if intermediate_count > 0:
            for i in range(intermediate_count):
                intermediate = intermediate_frame(g,new_g, (i+1.0) / (intermediate_count+1))
                intermediate.draw('test{0:04d}.png'.format(counter))
                counter += 1
        new_g.draw('test{0:04d}.png'.format(counter))
        print "Frame {0} Completed".format(counter)
        counter += 1
        g = new_g
        bp = new_bp
               
        
