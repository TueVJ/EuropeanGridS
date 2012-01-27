from zdcpf import *

def copper_flow():
    N = Nodes()
    N,F,lF = zdcpf(N,coop=0,copper=1)
    N.save_nodes('copper_nodes')
    save('./results/'+'copper_flows',F)
