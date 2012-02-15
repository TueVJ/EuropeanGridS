from  zdcpf import *

def compare_all_nodes(name1,name2,N=None,eps=1e-5):
    if N==None:
        N=[]
        for i in [1,2]:
            if (i == 1):
                N.append(Nodes(load_filename=name1))
        else:
            N.append(Nodes(load_filename=name2))

    for i in N[0]:
        for j in N[1]:
            if (i.label == j.label):
                test = compare_2_nodes(i,j,eps=eps)
                if (test == -1):
                    return -1

    return 1
            

def compare_2_nodes(n1,n2,eps=1e-5):
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
        return -1
    for i in range(n1.nhours):
        if (n1.load[i] != n2.load[i]):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in load at time',i, '!'
            print n1.load[i] , '!=', n2.load[i]
            return -1

    for i in range(n1.nhours):
        if (n1.normwind[i] != n2.normwind[i]):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in normwind at time',i, '!'
            print n1.normwind[i] , '!=', n2.normwind[i]
            return -1

    for i in range(n1.nhours):
        if (n1.normsolar[i] != n2.normsolar[i]):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in normsolar at time',i, '!'
            print n1.normsolar[i] , '!=', n2.normsolar[i]
            return -1

    for i in range(n1.nhours):
        if (n1.mismatch[i] != n2.mismatch[i]):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in mismatch at time',i, '!'
            print n1.mismatch[i] , '!=', n2.mismatch[i]
            return -1

    for i in range(n1.nhours):
        if (abs(n1.balancing[i] - n2.balancing[i]) >= eps):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in balancing at time',i, '!'
            print n1.balancing[i] , '!=', n2.balancing[i]
            return -1

    for i in range(n1.nhours):
        if (abs(n1.curtailment[i] - n2.curtailment[i]) >= eps):
            print 'Nodes ', n1.label, ' and ', n2.label, 'differ in curtailment at time',i, '!'
            print n1.curtailment[i] , '!=', n2.curtailment[i]
            return -1

    return 1

                        

                        

