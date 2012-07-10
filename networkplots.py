#! /usr/bin/env python
import networkx as nx
from zdcpf import *
from logistic_gammas import *
import matplotlib.colors as col
import numpy as np
import matplotlib.pyplot as plt

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


def get_node_colors(year=2035,step=2,combifit=False):
    ''' Depending on gamma values from logfit, assign a color from red
    (gamma=0) through yellow (gamma=0.5) to green (gamma=1) to each
    node => create a movie of development in time.'''
    gamma = get_basepath_gamma(year,step=step,combifit=combifit)
    # mix the color
    cl = []
    for gam in gamma:
        rp = int(round(get_positive(255.-get_positive(gam-0.5)*2.*256.))) # red
        gp = int(round(get_positive(255.-get_positive(0.5-gam)*2.*256.))) # green
        cl.append('#{0:02X}{1:02X}{2:02X}'.format(rp,gp,10))
    labels = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']
    color_dict = dict(zip(labels,cl))
    return color_dict


def get_color_map():
    cl = []
    for i in range(256):
        rp = int(get_positive(255-get_positive(i-128)*2)) # red
        gp = int(get_positive(255-get_positive(128-i)*2)) # green
        cl.append('#{0:02X}{1:02X}{2:02X}'.format(rp,gp,70))
    cmap3 = col.ListedColormap(cl,'indexed')
    # a=np.outer(np.arange(0,1,0.01),np.ones(10))   # pseudo image data
    # plt.imshow(a,aspect='auto',cmap=cm.get_cmap(cmap3),origin="lower")
    # plt.show()
    return cm.get_cmap(cmap3)


def drawgrid(N,node_colors=None,figname='weighted_graph.png',path='./figures/',title_='Weighted network'):
    c,G=AtoKh_graph(N)
    pos=nx.spring_layout(G)
    pos['FI']=[1.0,1.2]
    pos['SE']=[0.75,1.1]
    pos['NO']=[0.5,1.2]
    pos['DK']=[0.5,0.95]
    pos['PL']=[0.75,0.8]
    pos['GR']=[0.78,-0.05]
    pos['BG']=[0.98,0.03]
    pos['SK']=[0.95,0.59]
    pos['CZ']=[0.7,0.63]
    pos['HU']=[0.92,0.4]
    pos['RO']=[1.05,0.2]
    pos['RS']=[0.85,0.15]
    pos['BA']=[0.65,0.15]
    pos['HR']=[0.75,0.3]
    pos['IE']=[-0.15,0.925]
    pos['GB']=[0.05,0.9]
    pos['FR']=[0.06,0.60]
    pos['ES']=[-0.0,0.25]
    pos['PT']=[-0.15,0.15]
    pos['BE']=[0.25,0.75]
    pos['NL']=[0.37,0.85]
    pos['LU']=[0.325,0.575]
    pos['DE']=[0.45,0.7]
    pos['CH']=[0.4,0.45]
    pos['IT']=[0.3,0.2]
    pos['AT']=[0.58,0.45]
    pos['SI']=[0.55,0.28]

    for i,j in pos.items():
        for l in j:
            l*=1.8

    # hard work just to switch the stupid labels off in the first draw...
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']
    labels=dict(zip(ISO, ['' for name in ISO]))

    fig=figure(1); clf()
    if (node_colors == None):
        nx.draw_networkx(G,pos,node_size=500,node_color='b',facecolor=(1,1,1),labels=labels)
    else:
        colors = [node_colors[v] for v in G] # put colors in right order
        k=0
        # check that node order matches color order
        # check=True
        # for v in G:
        #     for l in range(len(ISO)):
        #         if (ISO[l] == v and node_colors[ISO[l]] != colors[k]):
        #             print ISO[l],v,node_colors[ISO[l]],colors[k]
        #             check=False
        #     k +=1
        # print check
        nx.draw_networkx(G,pos,node_size=500,node_color=colors,facecolor=(1,1,1),labels=labels)
    esmall=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=500]
    emid=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>500 and d['weight']<=1500]
    elarge=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1500]

    nx.draw_networkx_edges(G,pos,edgelist=esmall,width=1)#1
    nx.draw_networkx_edges(G,pos,edgelist=emid,width=3)#3
    nx.draw_networkx_edges(G,pos,edgelist=elarge,width=6)#6
    # now draw custom labels
    nx.draw_networkx_labels(G,pos,font_size=16,font_color='w',font_family='sans-serif',font_style='bold')
    axis('off')
    fig.suptitle(title_,fontsize=20)
    # add a color bar if colors were given
    if (node_colors != None):
        ax1 = fig.add_axes([0.86, 0.18, 0.05, 0.64])
        tcks = np.arange(0.,1.2,0.2)
        lbls = ['{0:.2f}'.format(t) for t in tcks]
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=get_color_map())
        cb1.set_ticks(tcks)
        cb1.set_ticklabels(lbls)
        cb1.set_label(r'VRES penetration $\gamma$',fontsize=14)
        # color_map = get_color_map()
        # pp=imshow([[],[]],cmap=color_map) # empty plot to get the color map
        # x = colorbar(pp, ticks = tcks,ticklabels = lbls)
    save_figure(figname,path=path) # save as png
    #show() # display


def logfit_movie(years=np.arange(2010,2050+1,1),step=2,combifit=False):
    # make pictures for each year
    N=Nodes()
    for year in years:
        color_dict=get_node_colors(year=year,step=step,combifit=combifit)
        figname = 'weighted_graph_year_{0:d}_step_{1:d}.png'.format(year,step)
        title_='Reference year {0:d}'.format(year)
        drawgrid(N,node_colors=color_dict,figname=figname,path='./figures/video_step_{0:d}/'.format(step),title_=title_)
    # glue them together to make a movie
    command = ('mencoder',
           'mf://figures/video_step_{0:d}/'.format(step)+'*.png',
           '-mf',
           'type=png:w=800:h=600:fps=5',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'video.avi')
    os.spawnvp(os.P_WAIT, 'mencoder', command)
    return

def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
	
    figure(fignumber)
    savefig(path + figname, dpi=dpi)
    print 'Saved figure:',path + figname
    sys.stdout.flush()

