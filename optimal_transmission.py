from zdcpf import *

def get_line_overload(copperfile='./results/copper_flows.npy',h=None,overloadfile='./results/overload.npy'):
    '''For each point in the time series, get the total line overload
    the copper flow would produce when line capacities have the values
    h.'''
    N = Nodes()
    K, h0, lF = AtoKh(N)
    h0 = h0[2:90] #because of dummy node etc.
    h0 = h0[::2]
    del N
    if (h == None):
        h=get_quant(0.82) # first quantile where no h is zero
        h = h[::2] # every line appears in two directions in h, pick one
    flows = np.load(copperfile)
    nhours = shape(flows)[1]
    nlinks = shape(flows)[0]
    if (len(h) != nlinks):
        print 'Error: Number of links = ',nlinks, 'number of capacities = ',len(h)
        return
    overload = []
    for i in range(nlinks):
        overload_i = []
        print 'Now looking at link '+lF[i][0]
        for j in range (nhours):
            overload_i.append(get_positive(flows[i][j]-h[i]))
        overload.append(overload_i)
    save(overloadfile,overload)
    return
        
def get_critical_times(overloadfile='./results/overload.npy',h=None,threshold=2.,numlines=5,critfile='./results/critical_times.npy'):
    '''Return the hours at which total overload is larger than <treshold>
    times the line capacity <h> for at least <numlines> lines.'''
    if (not os.path.exists(overloadfile)):
        get_line_overload(h=h,overloadfile=overloadfile)
    N = Nodes()
    K, h0, lF = AtoKh(N)
    h0 = h0[2:90] #because of dummy node etc.
    h0 = h0[::2]
    del N
    if (h == None):
        h=get_quant(0.82) # first quantile where no h is zero
        h = h[::2] # every line appears in two directions in h, pick one
    overload = np.load(overloadfile)
    ctimes = []
    nlinks = shape(overload)[0]
    nhours = shape(overload)[1]
    if (len(h) != nlinks):
        print 'Error: Number of links = ',nlinks, 'number of capacities = ',len(h)
        return
    for i in range(nhours):
        overlines = 0
        for j in range(nlinks):
            if (overload[j][i]/h[j] >= threshold):
                overlines += 1
        if (overlines >= numlines):
            print i,overlines
            ctimes.append(i)
    save(critfile,ctimes)
    return len(ctimes)
    
    
