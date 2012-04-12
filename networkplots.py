#! /usr/bin/env python
from zdcpf import *
import networkx as nx

def AtoKh_graph(N,pathadmat='./settings/admatold.txt'):
    Ad=np.genfromtxt(pathadmat,dtype='d')
    L=0
    G=nx.Graph()
    listFlows=[]
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    L+=1
    h=np.zeros(L*2)
    h=np.append(h,np.zeros(3*len(Ad)))
    L=0
    j=0
    i=0
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    h[2*L]=Ad[i,j]
                    h[2*L+1]=Ad[j,i]
                    if L>0: listFlows.append([str(N[j-1].label)+" to " +str(N[i-1].label), L-1])
                    G.add_edge(str(N[j-1].label),str(N[i-1].label),weight=max(h[2*L],h[2*L+1]))
                    L+=1
    return listFlows,G 


def drawgrid(N):
    c,G=AtoKh_graph(N)
    pos=nx.spring_layout(G)
    pos['FI']=[1.0,1.0]
    pos['SE']=[0.75,1.0]
    pos['NO']=[0.5,1.0]
    pos['DK']=[0.5,0.875]
    pos['PL']=[0.75,0.8]
    pos['GR']=[0.7,0.0]
    pos['BG']=[0.9,0.0]
    pos['SK']=[0.90,0.55]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RO']=[1.0,0.15]
    pos['RS']=[0.85,0.15]
    pos['BA']=[0.65,0.15]
    pos['HR']=[0.75,0.3]
    pos['IE']=[0.0,0.95]
    pos['GB']=[0.15,0.85]
    pos['FR']=[0.15,0.60]
    pos['ES']=[0.15,0.35]
    pos['PT']=[0.0,0.15]
    pos['BE']=[0.3,0.8]
    pos['NL']=[0.40,0.85]
    pos['LU']=[0.325,0.575]
    pos['DE']=[0.45,0.7]
    pos['CH']=[0.4,0.45]
    pos['IT']=[0.4,0.2]
    pos['AT']=[0.55,0.45]
    pos['SI']=[0.55,0.3]

    for i,j in pos.items():
        for l in j:
            l*=1.8



    nx.draw_networkx_nodes(G,pos,node_size=500,node_color='b',facecolor=(1,1,1))
    esmall=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=500]
    emid=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>500 and d['weight']<=1500]
    elarge=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1500]

    nx.draw_networkx_edges(G,pos,edgelist=esmall,width=2)#1
    nx.draw_networkx_edges(G,pos,edgelist=emid,width=2)#3
    nx.draw_networkx_edges(G,pos,edgelist=elarge,width=2)#6
    nx.draw_networkx_labels(G,pos,font_size=24,font_color='w',font_family='sans-serif')
    axis('off')
    savefig("weighted_graph.png") # save as png
    show() # display

