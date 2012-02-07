from zdcpf import *

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

def gamma_today():
    N = Nodes()
    N.set_alphas(0.7)
    for i in range(101):
        gamma = i*0.01
        print "Now calculating for gamma = ",gamma
        N.set_gammas(gamma)
        N,F,lF = zdcpf(N,coop=0,copper=0)
        name = 'today_linecap_gamma_%.2f' % gamma
        print name
        N.save_nodes(name+'_nodes')
        save('./results/'+name+'_flows',F)

def gamma_quant(quant=0.90):
    h0=get_quant(quant)
    N = Nodes()
    N.set_alphas(0.7)
    for i in range(101):
        gamma = i*0.01
        print "Now calculating for gamma = ",gamma
        N.set_gammas(gamma)
        N,F,lF = zdcpf(N,coop=0,h0=h0)
        name = 'quant_%.2f_gamma_%.2f' % (quant,gamma)
        print name
        N.save_nodes(name+'_nodes')
        save('./results/'+name+'_flows',F)
