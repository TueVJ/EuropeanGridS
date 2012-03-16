#! /usr/bin/env python
from pylab import *
from scipy import *
from numpy import concatenate as conc
import numpy as np
from cvxopt import matrix,solvers,spmatrix,sparse,spdiag
from time import time
import sys, os
from copy import deepcopy
import ctypes as ct

colors_countries = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers.
solvers.options['show_progress']=False

def get_positive(x):
    return x*(x>0.)  #Possibly it has to be x>1e-10.

class node:
    def __init__(self,path,fileName,ID):
        self.id = ID
        data = load(path + fileName)
        self.gamma = 1.0#float(setup[ID][0])
        self.alpha = 0.725#float(setup[ID][1])  #Alpha should be expanded to a vector.  completalpha() can be applied in update()
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
        
    def _load_nodes_(self,load_filename,path='./results/'):
        """Loads a Nodes instance from an npz file."""

        npzobj = load(path+load_filename)
        
        # for attribute in npzobj.files:
        #    for i in arange(len(self)):
        #        print self.cache[i],attribute,npzobj[attribute][i]
        #        setattr(self.cache[i],attribute,npzobj[attribute][i])

        for i in arange(len(self)):
            setattr(self.cache[i],'balancing',npzobj['balancing'][i])

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

def dcpowerflow(c,P,q,G,h,A,b):
    eps=1e-2 # small tolerance
    Nnodes=np.size(matrix(A),0)+1 # includes dummy node
    Nlinks=np.size(matrix(P),0)-2*Nnodes # includes dummy link
    # minimize balancing/curtailment first
    sol1=solvers.lp(c,G,h,A,b)
    x1=sol1['x']
    sumbal=sum(x1[Nlinks+1:Nlinks+Nnodes]) # minimal balancing
    sumcurt=sum(x1[Nlinks+Nnodes+1:Nlinks+2*Nnodes]) # minimal curtailment
    # add condition that balancing/curtailment be minimal to
    # inequality constraint (in equality constraint it causes sing. KKT mat.)
    h1=concatenate((h,matrix([sumbal+eps,sumcurt+eps])),0)
    h1=matrix(h1,tc='d')
    colind=arange(Nlinks+1,Nlinks+Nnodes)
    rowind=ndarray((len(colind),),int)*0
    balones=spmatrix(1.,rowind,colind,(1,Nlinks+2*Nnodes))
    colind=arange(Nlinks+Nnodes+1,Nlinks+2*Nnodes)
    curtones=spmatrix(1.,rowind,colind,(1,Nlinks+2*Nnodes))
    G1=sparse([G,balones,curtones],tc='d')
    P1=P
    #print P1[Nlinks+2,Nlinks+2]
    P1[Nlinks+1:,:]=P1[Nlinks+1:,:]*1e0
    q1=np.zeros(Nlinks+2*Nnodes)
    q1[Nlinks+1:]=1. # skip the dummy node
    #print q1[Nlinks]
    q1[Nlinks+Nnodes]=0. # skip the dummy node
    q1=2.e1*q1
    q1=matrix(q1,tc='d')
    sol=solvers.qp(P,q1,G1,h1,A,b)
    return sol['x']

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
                    if L>0: listFlows.append([str(N[j-1].label)+" to " +str(N[i-1].label), L-1])
                    L+=1
    K=spmatrix(K_values,K_row_indices,K_column_indices)
    # rowstring='{'
    # for i in range(len(K_values)):
    #     rowstring += str(K_row_indices[i]+1)+','+str(K_column_indices[i]+1)+' '
    # rowstring += '}'
    # print rowstring
    return K_values,h, listFlows               

def generatemat(N,admat='admat.txt',b=None,path='./settings/',copper=0,h0=None):
    K,h, listFlows=AtoKh(N,path+admat)
    if h0 != None: 
        h[2:90]=h0
    if b != None:
        for i in range(2*np.size(listFlows,0)): h[i]*=b
    Nnodes=np.size(matrix(K),0)
    Nlinks=np.size(matrix(K),1)
    # These numbers include the dummy node and link
    # With this info, we create the P matrix, sized (Nlinks+2Nnodes)^2
    # because a row is needed for each flow, and two for each node
    P1=spmatrix(1,range(Nlinks),range(Nlinks))
    P2=spmatrix(1e-6,range(2*Nnodes),range(2*Nnodes)) # and the bal/cur part has dif. coeffs
    P= spdiag([P1,P2])
    # Then we make the q vector, whose values will be changed all the time
    q=np.zeros(Nlinks+2*Nnodes)  # q has the same size and structure as the solution 'x'
    # Then we build the equality constraint matrix A
    # The stucture is more or less [ K | -I | I ]
    A1=spmatrix(1.,range(Nnodes),range(Nnodes))
    A=sparse([[K],[-A1],[A1]])
    # First row is cut so that A is nonsingular (see doc)
    A=A[1:,:]
    # b vector will be defined by the mismatches, in MAIN
    # Finally, the inequality matrix and vector, G and h.
    # Refer to doc to understand what the hell I'm doing, as I build
    # G...
    G1_values=[]
    G1_row_ind=[]
    G1_col_ind=[]
    for i in range(Nlinks):
      G1_values.extend([1,-1])
      G1_col_ind.extend([i,i])
    G1_row_ind=range(2*Nlinks)
    G1=spmatrix(G1_values,G1_row_ind,G1_col_ind)
    # to model copper plate, we forget about the effect of G matrix on the flows
    if copper == 1:
        G1=spmatrix([1,-1],[0,1],[0,0],(2*Nlinks,Nlinks))
   # G1 is ready, now we make G2
    G2=spmatrix([],[],[],(2*Nlinks,2*Nnodes))
    # G3 is built as [ 0 | -I | 0 ]
    G3=spmatrix(-1.,range(Nnodes),range(Nnodes))
    G3_hlp1=spmatrix([],[],[],(Nnodes,Nlinks))
    G3_hlp2=spmatrix([],[],[],(Nnodes,Nnodes))
    G3=sparse([[G3_hlp1],[G3],[G3_hlp2]])
    G4_values=[]
    G4_row_ind=[]
    G4_col_ind=[]
    for i in range(Nnodes):
      G4_values.extend([1,-1])
      G4_col_ind.extend([i,i])
    G4_row_ind=range(2*(Nnodes))
    G4=spmatrix(G4_values,G4_row_ind,G4_col_ind)
    G5=spmatrix([],[],[],(2*(Nnodes),Nlinks+Nnodes))
    G5=sparse([[G5],[G4]])
    G=sparse([[G1],[G2]])
    G=sparse([G,G3,G5])
    return P,q,G,h,A,K, listFlows

def runtimeseries(N,F,c,P,q,G,h,A,coop,lapse):
    if lapse==None:
        lapse=Nodes[0].mismatch.shape[0]
    Nlinks=np.size(matrix(F),0)
    Nnodes=np.size(matrix(A),0)
    start=time()
    b=np.zeros(Nnodes)
    b=matrix(b,tc='d')
    P_b=P[Nlinks+2:Nlinks+Nnodes+2,:]*1e6
    for t in range(lapse):
        for i in N:
            b[i.id]=i.mismatch[t]
            # from default, both curtailment and balancing have a minimum of 0.
            # in order to prevent export of curtailment, max curtailment is set to b
            h[2*Nlinks+Nnodes+5+2*i.id]=0
            if b[i.id]>0:
                h[2*Nlinks+Nnodes+5+2*i.id]=b[i.id]
        # then, we set the values of q_b and q_r for bal and cur, according to doc.
        # for Gorm's inequalities, we need f,L,delta
        f=P[0,0]
        L=Nnodes-1
        d=np.array(b)
        #d=b
        excess=np.dot(d.T,d>0)#[0][0]
        deficit=abs(np.dot(d.T,d<0))#[0][0])
        delta=min(excess,deficit)
        q_r=L*f*2*delta*0.5
        q_b=L*f*2*delta+q_r*(1.5)
        q[Nlinks+2:Nlinks+Nnodes+2]=q_b
        q[Nlinks+Nnodes+2:]=q_r
        if coop==1:
            P[Nlinks+2:Nlinks+Nnodes+2,:]=P_b*L*f*deficit*.99
        opt=dcpowerflow(c,P,q,G,h,A,b)   ########### Save relevant solution as flows
        for j in range(Nlinks):
            F[j][t]=opt[j+1]           
        for k in N:              ########### Save balancing at each node
            k.balancing[t]=opt[2+Nlinks+k.id]
            k.curtailment[t]=opt[3+Nlinks+Nnodes+k.id]  
        end=time()
        if (np.mod(t,1000)==0) and t>0:
            print "Elapsed time is ",round(end-start)," seconds. t = ",t," out of ",lapse
            sys.stdout.flush()
    end=time()
    #sys.stdout.flush()
    print "Calculation took ",round(end-start)," seconds."
    return N,F

def get_quant(quant=0.99,filename='results/copper_flows.npy'):
    outfile = 'results/linecap_quant_%.2f.npy' % quant
    if os.path.exists(outfile):
        hs = np.load(outfile)
        return hs
    f=np.load(filename)
    flows=[]
    for i in f:
        flows.append(i*(i>0.))
        flows.append(-i*(i<0.))
    a=np.zeros(len(flows))
    b=np.zeros(len(flows))
    hs=np.zeros(len(flows))
    for i in range(len(flows)):
        a=hist(flows[i],cumulative=True,bins=100,normed=True)
        for j in range(len(a[0])):
            if (a[0][j]>=quant):
                hs[i]=a[1][j]
                break
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

def zdcpf(N,admat='admat.txt',path='./settings/',coop=0,copper=0,lapse=None,b=None,h0=None):
    if lapse == None:
        lapse=N[0].nhours
    P,q,G,h,A,K, listFlows = generatemat(N,admat,b,path,copper,h0)
    Nnodes=np.size(matrix(K),0)-1
    Nlinks=np.size(matrix(K),1)-1
    F=np.zeros((Nlinks,lapse))
    P=sparse(P,tc='d')
    q=matrix(q,tc='d')
    G=sparse(G,tc='d')
    h=matrix(h,tc='d')
    A=sparse(A,tc='d')
    c=concatenate((1e-6*np.ones(Nlinks+1),np.ones(2*Nnodes+2)))
    c=matrix(c,tc='d')
    #print shape(matrix(G))
    N,F=runtimeseries(N,F,c,P,q,G,h,A,coop,lapse)
    return N,F, listFlows

def sdcpf(admat='admat.txt',path='./settings/',copper=0,lapse=None,b=None,h0=None):
    eps=1e-1
    N=Nodes()
    firststep=ct.CDLL('./balmin/libbalmin.so')
    firststep.balmin.restype=ct.c_double # default return type is int
    secondstep=ct.CDLL('./flowmin/libflowmin.so')
    secondstep.flowmin.restype=ct.POINTER(ct.c_double)

    if lapse == None:
        lapse=N[0].nhours
    kv,H,Lf=AtoKh(N) # dummy node has been deleted from admat.txt!!!
    h_neg=-H[1:88:2]
    h_pos=H[0:88:2]
    k=array([float(i) for i in kv])
    Nlinks=44
    Nnodes=27
    K=k.ctypes.data_as(ct.c_void_p)
    H_neg=h_neg.ctypes.data_as(ct.c_void_p)
    H_pos=h_pos.ctypes.data_as(ct.c_void_p)

    delta=np.zeros(Nnodes)

    start=time()
    for t in range(lapse):
        for i in N:
            delta[i.id]=i.mismatch[t]
        Delta=delta.ctypes.data_as(ct.c_void_p)
        MinBal=firststep.balmin(Delta,K,H_neg,H_pos)
        #print "MinBal is ", MinBal
        minbal=ct.c_double(MinBal+eps)
        MinFlow=secondstep.flowmin(Delta,K,H_neg,H_pos,minbal)
        minflows=np.zeros(Nlinks)
        for i in range(Nlinks):
            minflows[i]=MinFlow[i]
        #print "MinFlows are "
        #print minflows
        end=time()
        if (np.mod(t,2000)==0) and t>0:
            print "Elapsed time is %3.1f seconds. t = %u out of %u" % ((end-start), t, lapse)
            sys.stdout.flush()
    end=time()
    print "Calculation took %3.1f seconds." % (end-start)


