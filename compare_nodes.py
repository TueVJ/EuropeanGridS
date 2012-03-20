from  zdcpf import *

def compare_all_nodes(name1,name2,N=None,epsrel=1e-5,epsabs=5,silent=True):
    if N==None:
        N=[]
        for i in [1,2]:
            if (i == 1):
                N.append(Nodes(load_filename=name1))
        else:
            N.append(Nodes(load_filename=name2))

    balcheck=True
    curtcheck=True
    for i in N[0]:
        for j in N[1]:
            if (i.label == j.label):
                bcheck,ccheck = compare_2_nodes(i,j,epsrel=epsrel,epsabs=epsabs,silent=silent)
                if ~bcheck:
                    balcheck=False
                if ~ccheck:
                    curtcheck=False

    nhours = N[0][0].nhours
    Nnodes = len(N[0])
    for i in range(nhours):
        sumbal=zeros(2)
        sumcurt=zeros(2)
        for j in range(2):
            for n in N[j]:
                sumbal[j] += n.balancing[i]
                sumcurt[j] += n.curtailment[i]
        summe=sumbal[0]+sumbal[1]
        diff=2*abs(sumbal[0]-sumbal[1])/summe
        if (diff > epsrel and Nnodes*summe > epsabs):
            print 'Sum of balancing at time ',i,' differs!'
            print name1,':',sumbal[0]
            print name2,':',sumbal[1]
        summe=sumcurt[0]+sumcurt[1]
        diff=2*abs(sumcurt[0]-sumcurt[1])/summe
        if (diff > epsrel and Nnodes*summe > epsabs):
            print 'Sum of curtailment at time ',i,' differs!'
            print name1,':',sumcurt[0]
            print name2,':',sumcurt[1]

    del N
    return balcheck,curtcheck

def compare_2_nodes(n1,n2,epsrel=1e-5,epsabs=5,silent=True):
    if (n1.label != n2.label):
        print 'Nodes ', n1.label, ' and ', n2.label, 'differ in label!'
        print n1.label , '!=', n2.label
        return -1
    if (n1.id != n2.id):
        print 'Nodes ', n1.label, ' and ', n2.label, 'differ in id!'
        print n1.id , '!=', n2.id
        return -1
    if (n1.gamma != n2.gamma):
        print 'Nodes ', n1.label, ' and ', n2.label, 'differ in gamma!'
        print n1.gamma , '!=', n2.gamma
        return -1
    if (n1.alpha != n2.alpha):
        print 'Nodes ', n1.label, ' and ', n2.label, 'differ in alpha!'
        print n1.alpha , '!=', n2.alpha
        return -1
    if (n1.nhours != n2.nhours):
        print 'Nodes ', n1.label, ' and ', n2.label, 'differ in nhours!'
        print n1.nhours , '!=', n2.nhours
        return -1
    if (n1.mean != n2.mean):
        print 'Nodes ', n1.label, ' and ', n2.label, 'differ in mean!'
        print n1.mean , '!=', n2.mean
    for i in range(n1.nhours):
        if (n1.load[i] != n2.load[i]):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in load at time',i, '!'
            print n1.load[i] , '!=', n2.load[i]

    for i in range(n1.nhours):
        if (n1.normwind[i] != n2.normwind[i]):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in normwind at time',i, '!'
            print n1.normwind[i] , '!=', n2.normwind[i]

    for i in range(n1.nhours):
        if (n1.normsolar[i] != n2.normsolar[i]):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in normsolar at time',i, '!'
            print n1.normsolar[i] , '!=', n2.normsolar[i]

    for i in range(n1.nhours):
        if (n1.mismatch[i] != n2.mismatch[i]):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in mismatch at time',i, '!'
            print n1.mismatch[i] , '!=', n2.mismatch[i]

    balcheck=True
    for i in range(n1.nhours):
        diff = 2.*abs(n1.balancing[i] - n2.balancing[i])/(n1.balancing[i] + n2.balancing[i])
        summe=(n1.balancing[i] + n2.balancing[i])
        if ((diff > epsrel) and (summe > epsabs)):
            balcheck=False
            if not silent:
                print 'Nodes ', n1.label, ' and ', n2.label, 'differ in balancing at time',i, '!'
                print n1.balancing[i] , '!=', n2.balancing[i]
                print 'Relative difference: ',diff

    curtcheck=True
    for i in range(n1.nhours):
        diff = 2*abs(n1.curtailment[i] - n2.curtailment[i])/(n1.curtailment[i] + n2.curtailment[i])
        summe = n1.curtailment[i] + n2.curtailment[i]
        if (diff > epsrel and summe > epsabs):
            curtcheck=False
            if not silent:
                print 'Nodes ', n1.label, ' and ', n2.label, 'differ in curtailment at time',i, '!'
                print n1.curtailment[i] , '!=', n2.curtailment[i]
                print 'Relative difference: ',diff

    return balcheck,curtcheck

def compare_flows(name1,name2,epsrel=1.e-3,epsabs=20,silent=False):
    f1=np.load('./results/'+name1+'.npy')
    f2=np.load('./results/'+name2+'.npy')
    if (len(f1) != len(f2)):
        print 'Different number of transmission lines in files ',name1,' and ',name2,'!'
        return False
    nlines=len(f1)
    nhours=len(f1[0])
    flsum=np.zeros((2,nhours))
    flowcheck=True
    for i in range(nhours):
        for j in range(nlines):
            flsum[0,i] += (f1[j][i])*(f1[j][i])
            flsum[1,i] += (f2[j][i])*(f2[j][i])
            summe=f1[j][i]+f2[j][i]
            reldiff=2*abs(f1[j][i]-f2[j][i])/summe
            if (reldiff>epsrel and summe>epsabs):
                flowcheck = False
                if not silent:
                    print 'At time ',i, ' large difference in flows at index ',j
                    print name1,': ',f1[j][i]
                    print name2,': ',f2[j][i]
                    print 'Relative difference: ', reldiff
        summe=flsum[0,i]+flsum[1,i]
        reldiff=2*abs(flsum[0,i]-flsum[1,i])/summe
        if (reldiff>epsrel and summe>epsabs):
            print 'Sum of squared flows at time ',i,' differs'
            print name1,': ',flsum[0,i]
            print name2,': ',flsum[1,i]
            print 'Relative difference: ', reldiff
    return flowcheck

                        

