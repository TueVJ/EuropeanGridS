from zdcpf import *
from pylab import *
from scipy import *
import scipy.optimize as optimize
from numpy import concatenate as conc
import numpy as np
from cvxopt import matrix,solvers
from time import time
import os
from copy import deepcopy
import pickle

def copper_flow():
    N = Nodes()
    N.set_alphas(0.7)
    N.set_gammas(1.0)
    N,F,lF = zdcpf(N,coop=0,copper=1)
    N.save_nodes('copper_nodes')
    save('./results/'+'copper_flows',F)

def gamma_copper():
    N = Nodes()
    N.set_alphas(0.7)
    for i in range(101):
        gamma = i*0.01
        print "Now calculating for gamma = ",gamma
        N.set_gammas(gamma)
        N,F,lF = zdcpf(N,coop=0,copper=1)
        name = 'copper_gamma_%.2f' % gamma
        print name
        N.save_nodes(name+'_nodes')
        save('./results/'+name+'_flows',F)

def gamma_today(start=None):
    N = Nodes()
    N.set_alphas(0.7)
    if start != None:
        skip = start
    else:
        skip = 0
    for i in range(101):
        if (i < skip):
            continue
        gamma = i*0.01
        print "Now calculating for gamma = ",gamma
        N.set_gammas(gamma)
        N,F,lF = zdcpf(N,coop=0,copper=0)
        name = 'today_linecap_gamma_%.2f' % gamma
        print name
        N.save_nodes(name+'_nodes')
        save('./results/'+name+'_flows',F)

def gamma_quant(quant=0.90,start=None):
    h0=get_quant(quant)
    N = Nodes()
    N.set_alphas(0.7)
    if start != None:
        skip = start
    else:
        skip = 0
    for i in range(101):
        if (i < skip):
            continue
        gamma = i*0.01
        print "Now calculating for gamma = ",gamma
        N.set_gammas(gamma)
        N,F,lF = zdcpf(N,coop=0,h0=h0)
        name = 'quant_%.2f_gamma_%.2f' % (quant,gamma)
        print name
        N.save_nodes(name+'_nodes')
        save('./results/'+name+'_flows',F)

def get_balancing_vs_gamma(filename='quant_0.90_gamma'):
    gammas=arange(0,1.01,0.01)
    Bvsg=np.zeros((len(gammas),2))
    j=0
    for gamma in gammas:
        Bvsg[j,0]=gamma
        load_filename=(filename+'_%.2f_nodes.npz' % gamma)
        print load_filename
        N=Nodes(load_filename=load_filename)
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        Bvsg[j,1]=c
        j+=1
    save('./results/Bvsg_'+filename,Bvsg)
    return Bvsg

def plot_balancing_vs_gamma(filenames=['today_linecap_gamma','quant_0.40_gamma','quant_0.90_gamma','copper_gamma'],title_=r"homogenous increase in $\gamma\,$; $ \alpha_{\rm W}=0.7 $",label=['line capacities as of today',r'40$\,$% quantile line capacities',r'90$\,$% quantile line capacities','copper plate'],picname='balancing_vs_gamma.png'):
    figure(1); clf()
    cl = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers
    pp = []
    i=0
    for name in filenames:
        fname = './results/Bvsg_'+name+'.npy'
        print 'Processing '+fname
        if not os.path.exists(fname):
            data = get_balancing_vs_gamma(filename=name)
        else:
            data = np.load(fname)
        pp_x = array(data[:,0])
        pp_y = array(data[:,1]) # balancing
        pp_y = pp_y - (1-pp_x)  # excess balancing
        pp_ = plot(pp_x,pp_y,lw=1.5,color=cl[i])
        pp.extend(pp_)
        i += 1

    title(title_)
    axis(xmin=0,xmax=1,ymin=0,ymax=0.3)
    xlabel('share $\gamma$ of VRES in total electricity production')
    ylabel(r'excess balancing/av.h.l. (bal.$ - (1-\gamma)$)')

    pp_label = label
    leg = legend(pp,pp_label,loc='upper left');
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    save_figure(picname)

def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
	
    figure(fignumber)
    savefig(path + figname, dpi=dpi)
    print 'Saved figure:',path + figname
    sys.stdout.flush()
