#!/usr/bin/env python2
# coding: utf8

import argparse

import re
stem_re  = r'List of (\d+) stems'
stack_re = r'List of (\d+) stacks'
neighbor_re = r'Summary of structural features of (\d+) nucleotides'

def get_matrix(data):


    m_stem = re.search(stem_re, data)
    data_stems = ''
    stems = []
    if m_stem is not None:
        data_stems = data[m_stem.end():].split('****************************************************************************')[0]
        stems = data_stems.split('--------------------------------------------------------------------------')
    stems_output = []
    for stem in stems:
        stack1 = []
        stack2 = []
        for line in stem.split('\n'):
            line = line.split()
            if len(line) > 0:
                if line[0].isdigit():
                    stack1.append(line[1])
                    stack2.append(line[2])
        for i in range(len(stack1)-1):
            stems_output.append({'edge':(stack1[i], stack1[i+1]), 'type':'stem'} )
        for i in range(len(stack2)-1):
            stems_output.append({'edge':(stack2[i], stack2[i+1]), 'type':'stem'})



    m_stack = re.search(stack_re, data)
    data_stacks=''
    if m_stack is not None:
        data_stacks = data[m_stack.end():].split('****************************************************************************')[0]
    stacks_output = []
    for stack in data_stacks.split('\n'):
        stack = stack.strip().split()
        if len(stack)>0:
            if stack[0].isdigit():
                stack = stack[3].split(',')
                for i in range(len(stack)-1):
                    stacks_output.append({'edge':(stack[i], stack[i+1]), 'type':'stack'} )


    m_neighbor = re.search(neighbor_re, data) 
    data_neighbors=''
    if m_neighbor is not None:
        data_neighbors = data[m_neighbor.end():].split('****************************************************************************')[0]
    neighbor_output = []
    last_item = None
    for stack in data_neighbors.split('\n'):
        stack = stack.strip().split()
        if len(stack)>0:
            if stack[0].isdigit():
                if last_item is None:
                    last_item = stack[3]
                else:
                    if last_item.split('.')[:-3]==stack[3].split('.')[:-3]:
                        neighbor_output.append({'edge':( last_item, stack[3]), 'type':'neighbor'} )
                        last_item = stack[3]
                    else:
                        last_item = stack[3]


    edges = stacks_output+stems_output + neighbor_output


    a_minor_re = r'List of (\d+) A-minor motifs'
    m_a_minor = re.search(a_minor_re, data) 
    a_minor_data=''
    
    if m_a_minor is not None:
        a_minor_data = data[m_a_minor.end():].split('****************************************************************************')[0]

            
            
    

    a_minors = {}
    a_minor_vert = set()
    last_item = None
    for stack in a_minor_data.split('\n'):
        stack = stack.strip().split()
        if len(stack)>0:
            if stack[0].isdigit():
                a_minor_vert.add(stack[3])
                f_stack, mid_stack = stack[3].split('|')
                s_stack, th_stack  = mid_stack.split(',')
                if f_stack not in a_minors:
                    a_minors[f_stack] = []
                a_minors[f_stack].append(stack[3])
                if s_stack not in a_minors:
                    a_minors[s_stack] = []
                a_minors[s_stack].append(stack[3])
                if th_stack not in a_minors:
                    a_minors[th_stack] = []
                a_minors[th_stack].append(stack[3])
    matrix = {}
    matrix_type = {}
    for i, x in enumerate(a_minor_vert):
        matrix[x]={}
        matrix_type[x] = {}
        for y in a_minor_vert:
            if x<y:
                matrix[x][y]=[]
                matrix_type[x][y]=[]
                

    edges = stacks_output+stems_output + neighbor_output 


    for e in edges:
        a, b = min(e['edge']), max(e['edge'])
        if a in a_minors and b in a_minors:
            for x in a_minors[a]:
                if x in a_minor_vert:
                    for y in a_minors[b]:
                        if y in a_minor_vert:
                            if x < y:
                                matrix[x][y].append(e['type']+'_'+a+'_'+b)
                            if x > y:
                                matrix[y][x].append(e['type']+'_'+a+'_'+b)


    for x in matrix:
        matrix[x][x]=[]                                
                                
    for c in a_minors:
        if len(a_minors[c])>1:
            for i in range(len(a_minors[c])):
                if a_minors[c][i] in a_minor_vert:
                    for j in range(len(a_minors[c])):
                        if a_minors[c][j] in a_minor_vert:
                            if a_minors[c][i]<a_minors[c][j]:
                                matrix[a_minors[c][i]][a_minors[c][j]].append('a_minor_'+c)
                            else: 
                                matrix[a_minors[c][j]][a_minors[c][i]].append('a_minor_'+c)

    for x in matrix:
        for y in matrix:
            if x < y:
                matrix[y][x]=matrix[x][y]            
        matrix[x][x]=[]
    
    for x in matrix:
        for y in matrix:
            for t in matrix[x][y]:
                if t.find('a_minor_')!= -1:
                    matrix[y][x] = []
                    matrix[x][y] = []
    
    
    for x in matrix_type:
        for y in matrix_type:
            if x<y:
                eat = {}
                for e in matrix[x][y]:
                    et, e1, e2 = e.split('_')
                    ax, nnx = x.split('|')
                    lnx, rnx  = nnx.split(',')
                    ay, nny = y.split('|')
                    lny, rny  = nny.split(',')

                    type_str = ''
                    if lnx in (e1, e2):
                        type_str += 'L'
                    elif rnx in (e1, e2):
                        type_str += 'R'
                    elif ax in (e1, e2):
                        type_str += 'A'
                    if lny in (e1, e2):
                        type_str += 'L'
                    elif rny in (e1, e2):
                        type_str += 'R'
                    elif ay in (e1, e2):
                        type_str += 'A'

                    if type_str not in eat:
                        eat[type_str] = []
                    eat[type_str].append(et[0])
                    eat[type_str].sort()
                for e in eat:
                    new_t = ''.join(eat[e])+e
                    if x<y and new_t not in matrix_type[x][y]:
                        matrix_type[x][y].append(new_t)
                        matrix_type[x][y].sort()
                    if y<x:
                        matrix_type[y][x].append(new_t)
                        matrix_type[y][x].sort()
                        
                        
    for x in matrix_type:
        for y in matrix_type:
            if x < y:
                matrix_type[y][x]=[item[:-2]+item[-1]+ item[-2] for item in matrix_type[x][y]]
        matrix_type[x][x]=[]

    return matrix, matrix_type, len(a_minor_vert)




has_cycle = 0

order = 0

def sub_dfs(u, matrix, visited, component):
    global has_cycle
    global order
    visited[u][0] = 1
    visited[u][2] = order
    order += 1
    numerator = {}
    
    for v in matrix[u]:
        if len(matrix[u][v])>0:
            a, b = min((u, v)), max((u, v))
            if {'edge': (a, b), 'links': matrix[a][b], 'edge_pattern':visited[a][1] + '_'+visited[b][1]} not in component:
                if visited[v][0] == 1:
                    has_cycle+=1
                component.append({'edge': (a, b), 'links': matrix[a][b], 'edge_pattern': visited[a][1] + '_'+visited[b][1]})
            if visited[v][0] == 0:
                sub_dfs(v, matrix, visited, component)
    visited[u][0] = 2
    

def dfs(matrix):
    global has_cycle
    global order
    text=[]
    cycle_text = []
    line_components = []
    cycle_components = []
    triv_components =[]
    visited = {ver:[0, 'v'+str(i), 0] for i, ver in enumerate(matrix)}
    component = []
    for v in matrix:            
        if not visited[v][0]:
            component = []
            has_cycle = 0
            sub_dfs(v, matrix, visited, component)
            if has_cycle > 0:
                cycle_components.append(component)
            elif len(component)>0:
                line_components.append(component)
            if len(component) == 0:
                triv_components.append(v)
            
    order = 0
    has_cycle = 0
    return  cycle_components, line_components, triv_components
   
def get_type(c):
    ver = set()
    for e in c:
        ver.add(e['edge'][0])
        ver.add(e['edge'][1])
    ad = set()
    pr = set()
    for v in ver:
        a, p = v.split('|')
        ad.add(a)
        pr.add(p)
    return  len(ver), len(ad), len(pr)


def get_comp(file_name):
    results = []
    my_matrix = get_matrix(open(file_name).read())
    cycle_components, line_components, triv_components = dfs(my_matrix[0])
    for c in  cycle_components:
        for i, e in enumerate(c):
            e['type_pattern'] = [t for t in set(my_matrix[1][e['edge'][0]][e['edge'][1]])]
            c[i]=e
        comp = {}
        comp['comp'] = c
        comp['type'] = get_type(c)
        comp['is_cycle'] = True
        results.append(comp)
    for c in  line_components:
        for i, e in enumerate(c):
            e['type_pattern'] = [t for t in set(my_matrix[1][e['edge'][0]][e['edge'][1]])]
            c[i]=e
        comp = {}
        comp['comp'] = c
        comp['type'] = get_type(c)
        comp['is_cycle'] = False
        results.append(comp)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_name', help='file_name')
    args = parser.parse_args()
    print args.file_name
    print get_comp(args.file_name)


