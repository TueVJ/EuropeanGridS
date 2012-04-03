#! /usr/bin/env python
import gc
from pylab import *
from scipy import *
from scipy.sparse import coo_matrix # analog of prev. used cvxopt.spmatrix
from numpy import concatenate as conc
import numpy as np
from time import time
import sys, os
from copy import deepcopy
import ctypes as ct

colors_countries = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers.

def get_positive(x):
    return x*(x>0.)  #Possibly it has to be x>1e-10.

class node:
    def __init__(self,path,fileName,ID):
        self.id = ID
        data = load(path + fileName)
        self.gamma = 1.0
        self.alpha = 0.7
        self.load = 1000*array(map(double,data['L']))
        self.nhours = len(self.load)
        self.normwind = array(map(double,data['Gw']))
        self.normsolar = array(map(double,data['Gs']))
        self.mean = mean(self.load)
        self.balancing = np.zeros(self.nhours)
        self.curtailment = np.zeros(self.nhours)
        self.label = data['datalabel']
        self.mismatch = None
        self.colored_import = None #Set using self.set_colored_i_import()
        del data.f
        data.close()
        self._update_()
        self.gen=np.zeros(self.nhours)        

    def _update_(self):
        self.mismatch=(self.get_wind()+self.get_solar())-self.load
    
    def get_import(self):
        """Returns import power time series in units of MW."""
        return get_positive(get_positive(-self.mismatch) - self.balancing) #Balancing is exported if it exceeds the local residual load.
        
    def get_export(self):	
        """Returns export power time series in units of MW."""
        return get_positive(self.mismatch) - self.curtailment + get_positive(self.balancing - get_positive(-self.mismatch))

    def get_localRES(self):
        """Returns the local use of RES power time series in units of MW."""
        return self.get_wind() + self.get_solar() - self.curtailment  - self.get_export()
		
    def get_localBalancing(self):
        """Returns the local use of balancing power time series in units of MW."""
        return get_positive(-self.mismatch) - self.get_import()
		
    def get_wind(self):
        """Returns wind power time series in units of MW."""
        return self.mean*self.gamma*self.alpha*self.normwind
	
    def get_solar(self):
        """Returns solar power time series in units of MW."""
        return self.mean*self.gamma*(1.-self.alpha)*self.normsolar
			
    def set_gamma(self,gamma,operation='='):
        if operation == '=':
            self.gamma = gamma
        else:
            self.gamma *= gamma
        self._update_()
    
    def set_alpha(self,alpha):
        self.alpha=alpha
        self._update_()

    def set_colored_import_i(self,i,colored_import_i):
        if self.colored_import == None:
            self.colored_import = np.zeros((len(colored_import_i),len(self.load)))
			
        self.colored_import.transpose()[i] = colored_import_i

class Nodes:
    def __init__(self,path='./data/',files=['ISET_country_AT.npz','ISET_country_FI.npz','ISET_country_NL.npz','ISET_country_BA.npz','ISET_country_FR.npz','ISET_country_NO.npz','ISET_country_BE.npz','ISET_country_GB.npz','ISET_country_PL.npz','ISET_country_BG.npz','ISET_country_GR.npz','ISET_country_PT.npz','ISET_country_CH.npz','ISET_country_HR.npz','ISET_country_RO.npz','ISET_country_CZ.npz','ISET_country_HU.npz','ISET_country_RS.npz','ISET_country_DE.npz','ISET_country_IE.npz','ISET_country_SE.npz','ISET_country_DK.npz','ISET_country_IT.npz','ISET_country_SI.npz','ISET_country_ES.npz','ISET_country_LU.npz','ISET_country_SK.npz'],load_filename=None):
        self.cache=[]
        for i in range(len(files)):
            n=node(path,files[i],i)
            self.cache=append(self.cache,n)
        F=np.zeros((size(files),self.cache[0].nhours))
        
        if load_filename != None:
            self._load_nodes_(load_filename,path='./results/')

    def __getitem__(self,x):
        return self.cache[x]
        
    def __len__(self):
        return len(self.cache)

    def set_gammas(self,value):
        # to change a single node's gamma, just write XX.set_gamma(yy)
        # 'value' can be a single number or a vector
        if size(value)==1:
            for i in self.cache: i.set_gamma(value)
        elif size(value)!=np.size(self.cache):
            print "Wrong gamma vector size. ", np.size(value,0)," were  received, ",size(self.cache)," were expected."
        else:
            for i in self.cache:
                i.set_gamma(value[i.id])

    def set_alphas(self,value):
        # to change a single node's alpha, just write XX.set_gamma(yy)
        # 'value' can be a single number or a vector
        if size(value)==1:
            for i in self.cache: i.set_alpha(value)
        elif size(value)!=size(self.cache):
            print "Wrong gamma vector size. ", size(value,0)," were  received, ",size(self.cache)," were expected."
        else:
            for i in self.cache:
                i.set_alpha(value[i.id])

    def save_nodes(self,filename,path='./results/'):
        """Saves the contents of a Nodes instance to a npz file."""
        
        attribute = dir(self[0])
        save_str = []
        #Determine which attributes to be saved
        for attribute in dir(self[0]):
            if attribute[0]=='_':
                continue
            elif is_numlike(getattr(self[0],attribute)) or is_string_like(getattr(self[0],attribute)):
                print attribute
                save_str.append(attribute + '=' + 'array([self[i].'+attribute+' for i in arange(len(self))])')

        #Write save file
        eval('savez(path+filename,'+','.join(save_str)+')')

        print 'Saved nodes to file: ', path+filename
        sys.stdout.flush()
        
    def save_nodes_small(self,filename,path='./results/'):
        """Saves the non-redundant contents of a Nodes instance to a
        npz file. A little less convenient, but much smaller."""
        
        attribute = dir(self[0])
        save_str = []
        #Determine which attributes to be saved
        for attribute in dir(self[0]):
            if attribute[0]=='_':
                continue
            if (attribute=='load' or attribute=='normwind' or attribute=='normsolar' or attribute=='mismatch' or attribute=='gen'):
                continue
            elif is_numlike(getattr(self[0],attribute)) or is_string_like(getattr(self[0],attribute)):
                print attribute
                save_str.append(attribute + '=' + 'array([self[i].'+attribute+' for i in arange(len(self))])')

        #Write save file
        eval('savez(path+filename,'+','.join(save_str)+')')

        print 'Saved nodes to file: ', path+filename
        sys.stdout.flush()
        
    def _load_nodes_(self,load_filename,path='./results/'):
        """Loads a Nodes instance from an npz file."""

        npzobj = load(path+load_filename)
        
        # for attribute in npzobj.files:
        #    for i in arange(len(self)):
        #        print self.cache[i],attribute,npzobj[attribute][i]
        #        setattr(self.cache[i],attribute,npzobj[attribute][i])

        for i in arange(len(self)):
            setattr(self.cache[i],'balancing',npzobj['balancing'][i])
            setattr(self.cache[i],'curtailment',npzobj['curtailment'][i])

        del npzobj.f
        npzobj.close()

        for n in self.cache:
            n._update_()
    
        print 'Loaded nodes from file: ', path+load_filename
        sys.stdout.flush()

    def add_colored_import(self, F, node_id=None, incidence_matrix='incidence.txt', path='./settings/',lapse=None):
	"""Type N.add_colored_stuff(F,lapse=x)"""
        if lapse==None:
            lapse=self.cache[0].mismatch.shape[0]
    	
        if type(incidence_matrix)==str:
            K = genfromtxt(path+incidence_matrix,delimiter=",",dtype='d') #Incidence matrix.
            K = K[1:].transpose()[1:].transpose() #Remove dummy row/column.
        else:
            K = incidence_matrix
    			
        for t in arange(lapse):
    	
            export_ = array([self.cache[i].get_export()[t] for i in arange(len(self.cache))])
            import_ = array([self.cache[i].get_import()[t] for i in arange(len(self.cache))])
    	
            FF, C = get_colored_flow(F.transpose()[t], copy(export_), incidence_matrix=K)
    	
            CC = C*kron(ones((K.shape[0],1)),import_)
    	
            #Update Node(s)
            if node_id == None:
                for node_id_ in arange(len(self.cache)):
                    self.cache[node_id_].set_colored_import_i(t,CC.transpose()[node_id_])
            else:
                self.cache[node_id].set_colored_import_i(t,CC.transpose()[node_id])

def get_colored_flow(flow, export, incidence_matrix='incidence.txt'):
    """flow: vector of flows at time t. export: vector of export at each node at time t."""
    if type(incidence_matrix)==str:
        K = np.genfromtxt(incidence_matrix,delimiter=",",dtype='d') #Incidence matrix.
        K = K[1:].transpose()[1:].transpose() #Remove dummy row/column.
    else:
        K = incidence_matrix
    Kf = K*kron(ones((K.shape[0],1)),-flow) #Modified incidence matrix that has positive values for the flow into a node.
    FF = array(mat(abs(K))*mat(Kf).transpose()) #Flow matrix with positive values indicating the positive flow into a node.
    FF = get_positive(floor(FF))
        #"Column sum" = 1
    for i in arange(FF.shape[1]):
        sum_ = (FF.transpose()[i].sum() + export[i])
        FF.transpose()[i] = FF.transpose()[i]/sum_
        export[i] = export[i]/sum_
    	
    	#Calculate color matrix	
    try:	
        C = -mat(diag(export))*inv(mat(FF)-mat(eye(FF.shape[0])))	
    except LinAlgError:
        print "Error (dfkln387c): Singular matrix"
        print mat(FF)-mat(eye(FF.shape[0]))
    
        C = np.zeros(FF.shape)

    return array(FF), array(C)   

def AtoKh(N,pathadmat='./settings/admat.txt'):
    Ad=np.genfromtxt(pathadmat,dtype='d')
    L=0
    listFlows=[]
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0: 
                    L+=1
    K_values=[]
    K_column_indices=[]
    K_row_indices=[]
    h=np.zeros(L*2)
    h=np.append(h,np.zeros(3*len(Ad)))
    L=0
    j=0
    i=0
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0: 
                    K_values.extend([1,-1])
                    K_column_indices.extend([L,L])
                    K_row_indices.extend([j,i])
                    h[2*L]=Ad[i,j]
                    h[2*L+1]=Ad[j,i]
                    listFlows.append([str(N[j].label)+" to " +str(N[i].label), L])
                    L+=1
    K=coo_matrix((K_values,(K_row_indices,K_column_indices)))
    # rowstring='{'
    # for i in range(len(K_values)):
    #     rowstring += str(K_row_indices[i]+1)+','+str(K_column_indices[i]+1)+' '
    # rowstring += '}'
    # print rowstring
    return K, K_values, h, listFlows               

def sdcpf(N,admat='admat.txt',path='./settings/',copper=0,lapse=None,b=None,h0=None):
    Nlinks=44
    Nnodes=27
    firststep=ct.CDLL('./balmin/libbalmin.so')
    firststep.balmin.restype=ct.c_double # default return type is int
    ufirststep=ct.CDLL('./ubalmin/libubalmin.so')
    ufirststep.ubalmin.restype=ct.c_double # default return type is int
    secondstep=ct.CDLL('./flowmin/libflowmin.so')
    secondstep.flowmin.restype=ct.c_int # just for fun
    usecondstep=ct.CDLL('./uflowmin/libuflowmin.so')
    usecondstep.uflowmin.restype=ct.c_int # just for fun

    if lapse == None:
        lapse=N[0].nhours
    km,kv,H,Lf=AtoKh(N) # dummy node has been deleted from admat.txt!!!
    km = np.matrix(km.todense())
    if (h0 != None):
        H=h0
    h_neg=-H[1:88:2]
    h_pos=H[0:88:2]
    #print 'h_pos: ',shape(h_neg),h_neg
    if (copper == 1):
        h_neg=-1.e6*np.ones(Nlinks)
        h_pos=1.e6*np.ones(Nlinks)
    #print 'h_pos: ',shape(h_neg),h_neg
    flw=np.zeros(Nlinks)
    k=array([float(i) for i in kv])
    K=k.ctypes.data_as(ct.c_void_p)
    H_neg=h_neg.ctypes.data_as(ct.c_void_p)
    H_pos=h_pos.ctypes.data_as(ct.c_void_p)
    Flw=flw.ctypes.data_as(ct.POINTER(ct.c_double)) # for return by ref

    delta=np.zeros(Nnodes)

    F=np.zeros((Nlinks,lapse)) # save flows and deltas to calc bal and curt later
    deltas=np.zeros((Nnodes,lapse))
    eps=1e-1
    start=time()
    if (copper == 0):
        for t in range(lapse):
            for i in N:
                delta[i.id]=i.mismatch[t]
            deltas[:,t]=delta
            Delta=delta.ctypes.data_as(ct.c_void_p)
            MinBal=firststep.balmin(Delta,K,H_neg,H_pos)
            # print "MinBal is ", MinBal
            minbal=ct.c_double(MinBal+eps)
            dummy=secondstep.flowmin(Delta,K,H_neg,H_pos,minbal,Flw)
            for i in range(Nlinks):
                F[i,t]=Flw[i]
            # print "MinFlows are "
            # print F[:,t]
            end=time()
            if (np.mod(t,2073)==0) and t>0:
                 print "Elapsed time is %3.1f seconds. t = %u out of %u" % ((end-start), t, lapse)
                 sys.stdout.flush()
        end=time()
        print "Calculation took %3.1f seconds." % (end-start)
        sys.stdout.flush()
    else: # use special unbounded C-fctns for unbounded copper flow
        for t in range(lapse):
            for i in N:
                delta[i.id]=i.mismatch[t]
            deltas[:,t]=delta
            Delta=delta.ctypes.data_as(ct.c_void_p)
            MinBal=ufirststep.ubalmin(Delta,K)
            # print "MinBal is ", MinBal
            minbal=ct.c_double(MinBal+eps)
            dummy=usecondstep.uflowmin(Delta,K,minbal,Flw)
            for i in range(Nlinks):
                F[i,t]=Flw[i]
            # print "MinFlows are "
            # print F[:,t]
            end=time()
            if (np.mod(t,2073)==0) and t>0:
                print "Elapsed time is %3.1f seconds. t = %u out of %u" % ((end-start), t, lapse)
                sys.stdout.flush()
        end=time()
        print "Calculation took %3.1f seconds." % (end-start)
        sys.stdout.flush()

    # no matter how the flows were obtained, we still have to calc bal and curt
    start2=time()
    for t in range(lapse):
        tmp=np.array(F[:,t])
        tmp2=np.array(deltas[:,t])
        tmp3=transpose(tmp2-dot(km,tmp))
        balancing=[-tmp3[i,0] if tmp3[i,0]<0 else 0. for i in range(len(tmp3))]
        curtailment=[tmp3[i,0] if tmp3[i,0]>0 else 0. for i in range(len(tmp3))]
        for i in N:
            i.balancing[t] = balancing[i.id]
            i.curtailment[t] = curtailment[i.id]
    end=time()
    print "Assigning balancing and curtailment took %3.1f seconds." % (end-start2)
    print "Complete calculation took %3.1f seconds." % (end-start)
    return F

def get_quant(quant=0.99,filename='results/copper_flows.npy'):
    outfile = 'results/linecap_quant_%.4f.npy' % quant
    if os.path.exists(outfile):
        hs = np.load(outfile)
        return hs
    f=np.load(filename)
    flows=[]
    for i in f:
        flows.append(i*(i>0.))
        flows.append(-i*(i<0.))
    a=np.zeros(len(flows))
    hs=np.zeros(len(flows))
    for i in range(len(flows)):
        a=hist(flows[i],cumulative=True,bins=500,normed=True)
        for j in range(len(a[0])):
            if (a[0][j]>=quant):
                hs[i]=a[1][j]
                break
        clf()
    np.save(outfile,hs)
    return hs

def show_hist(link,filename='results/copper_flows.npy',e=1,b=500):
    f=np.load(filename)
    flows=[]
    for i in f:
        flows.append(i)
#        flows.append(i*(i>e))
#        flows.append(i*(i<-e))
#    a=hist(flows[link],bins=b,range=[0.1,flows[link].max()],normed=1,histtype='stepfilled')
#    b=hist(flows[link+1],bins=b,range=[flows[link+1].min(),-0.1],normed=1,histtype='stepfilled')
    a=hist(flows[link],bins=b,normed=1,histtype='stepfilled')
    show()


def find_balancing_reduction_quantiles(reduction=[0.50,0.90],eps=1.e-3,guess=[0.885,0.98],stepsize=0.01,file_copper='copper_nodes.npz',file_notrans='homogenous_gamma_1.00_linecap_0.40Q_nodes.npz',gamma=1.,alpha=None,save_filename=None):

    '''Loop over different quantile line capacities until the quantile
    is found that leads to a reduction of balancing by <reduction>
    times what is possible, with a relative uncertainty of
    <eps>. <guess> specifies your first guess for the quantile.'''

    try:
        dmy=len(reduction)
        dmy=len(guess)
    except TypeError:
        print "reduction and guess must be (possibly one-element) lists!"
        return
    if (len(reduction) != len(guess)):
        print "Number of reduction quantile is %u, number of first guesses is %u!" % (len(reduction),len(guess))
        return

    N=Nodes(load_filename=file_notrans)
    a=0.; b=0.
    for i in N:
        a+=sum(i.balancing)
        b+=i.mean*i.nhours
    balmax=a/b
    del N
    N=Nodes(load_filename=file_copper)
    a=0.; b=0.
    for i in N:
        a+=sum(i.balancing)
        b+=i.mean*i.nhours
    balmin=a/b
    del N

    # check if there is a significant reduction possible at all
    if (2.*(balmax-balmin)/(balmax+balmin)<eps):
        if (save_filename != None):
            for i in range(len(reduction)):
                # print save_filename[i]
                # sys.stdout.flush()
                link_flows=save_filename[i]+'_flows.npy'
                link_nodes=save_filename[i]+'_nodes.npz'
                target_flows=file_notrans[:-10]+'_flows.npy'
                # print link_flows, target_flows
                # print link_nodes, file_notrans
                os.symlink(file_notrans,'./results/'+link_nodes)
                os.symlink(target_flows,'./results/'+link_flows)
        quant=np.zeros(len(reduction))
        return quant

    quant=guess # initial guess
    for i in range(len(reduction)):
        baltarget=balmin+(1.-reduction[i])*(balmax-balmin)
        print '%9s %10s %10s %10s' % ('reduction','balmin','balmax','baltarget')
        print '%9.4f %10.7f %10.7f %10.7f' % (reduction[i],balmin,balmax,baltarget)
        step=stepsize
        olddist=0.
        balreal=0.
        while True:
            h=get_quant(quant[i])
            # memory leak when the same N is used over and
            # over again (why?) -> del and new all the time
            N=Nodes()
            if (alpha == None):
                N.set_alphas(0.7)
            else:
                N.set_alphas(alpha)
            N.set_gammas(gamma)
            F=sdcpf(N,h0=h)
            a=0.; b=0.
            for j in N:
                a+=sum(j.balancing)
                b+=j.mean*j.nhours
            balreal=a/b
            reldist=abs(1.-balreal/baltarget)
            dist=baltarget-balreal
            if (reldist < eps or step<0.00125):
                print '%12s %13s %14s %9s %9s %9s' % ('distance','old distance','relative dist.','quantile','stepsize','balreal')
                print '%12.8f %13.8f %14.4f %9.4f %9.6f %9.7f' % (dist, olddist, reldist,quant[i],step,balreal)
                # if a filename is provided, save minimal nodes and flows
                if (save_filename != None): 
                    N.save_nodes_small(save_filename[i]+'_nodes')
                    np.save('./results/'+save_filename[i]+'_flows',F)
                print 'gc.get_count before collection ',gc.get_count()
                gc.collect()
                print 'gc.get_count after collection ',gc.get_count()
                del N #not here, but in same scope where it was constructed!
                break
            if (dist*olddist<0.): # sign change = we passed the perfect point! now reduce step size
                step=step/2.
            if dist<0:
                quant[i] +=step
            if dist>0:
                quant[i] -=step
            if (quant[i]>=1.): # we are clearly overshooting -> go back
                step=step/2.
                quant[i]=0.9999
            print '%12s %13s %14s %9s %9s %9s' % ('distance','old distance','relative dist.','quantile','stepsize','balreal')
            print '%12.8f %13.8f %14.4f %9.4f %9.6f %9.7f' % (dist, olddist, reldist,quant[i],step,balreal)
            olddist=dist
            print 'gc.get_count before collection ',gc.get_count()
            gc.collect()
            print 'gc.get_count after collection ',gc.get_count()
            del N
    return quant #, 1.-(balreal-balmin)/(balmax-balmin)
