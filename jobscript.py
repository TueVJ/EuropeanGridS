from zdcpf import *
import scipy.optimize as optimize
from logistic_gammas import *
from Database_v2 import * # only needed for optimal alphas
from mpl_toolkits.axes_grid1 import host_subplot
import matplotlib.pyplot as plt


######################################################
########### Generate the results #####################
######################################################

def copper_flow():
    N = Nodes()
    N.set_alphas(0.7)
    N.set_gammas(1.0)
    F = sdcpf(N,copper=1)
    N.save_nodes('copper_nodes')
    np.save('./results/'+'copper_flows',F)


def gamma_homogenous(linecap='copper',start=None,stop=None,alpha=None,gpercent=None):
    if (linecap == 'copper'):
        copper = 1
        h0 = None
    elif (linecap == '0.40Q'):
        h0 = get_quant(0.40)
        copper = 0
    elif (linecap == '0.90Q'):
        h0 = get_quant(0.90)
        copper = 0
    elif (linecap == 'today'):
        h0 = None
        copper = 0
    elif (linecap == '0.6_times_0.90Q'):
        h0 = 0.6*get_quant(0.90)
        copper = 0
    elif (linecap == '0.8_times_0.90Q'):
        h0 = 0.8*get_quant(0.90)
        copper = 0
    else:
        print 'invalid linecap. abort'
        return
    print 'linecap: ',linecap

    if start != None:
        skip = start
    else:
        skip = 0
    if stop != None:
        skip_end = stop
    else:
        skip_end=101

    N = Nodes()
    if (alpha==None):
        N.set_alphas(0.7)
    else:
        N.set_alphas(alpha)
    if (gpercent==None):
        gvals= arange(skip,skip_end,1)
    else:
        gvals=gpercent
    for gval in gvals:
        gamma = gval*0.01
        print "Now calculating for gamma = ",gamma
        N.set_gammas(gamma)
        F = sdcpf(N,copper=copper,h0=h0)
        name = ('homogenous_gamma_%.2f' % gamma)+'_linecap_'+linecap
        if (alpha != None):
            name += '_alpha_%.2f' % alpha
        print name
        N.save_nodes_small(name+'_nodes')
        np.save('./results/'+name+'_flows',F)


def gamma_homogenous_balred(red=[0.50,0.90],path='./results/',guessfile='homogenous_gamma_balred_quantiles',notrans_file='Bvsg_homogenous_gamma_linecap_0.40Q',copper_file='Bvsg_homogenous_gamma_linecap_copper',start=None,stop=None):
    if start != None:
        skip = start
    else:
        skip = 0
    if stop != None:
        skip_end = stop
    else:
        skip_end=101

    gvals= arange(skip,skip_end,1)
    quantiles=[]
    # guessfile += '.npy'
    # guesses = np.load(path+guessfile)
    cfile = copper_file + '.npy'
    balmins=np.load(path+cfile)[:,1]
    #print balmins
    nfile = notrans_file + '.npy'
    balmaxes = np.load(path+nfile)[:,1]
    #print balmaxes
    last_quant=np.zeros(2)
    for gval in gvals:
        gamma = gval*0.01
        print "Now calculating for gamma = ",gamma
        save_filename = []
        for i in range(len(red)):
            save_filename.append(('homogenous_gamma_%.2f_balred_%.2f' % (gamma,red[i])))
        balmax=balmaxes[gval]
        balmin=balmins[gval]
        guess=[0.0,0.0]
        f_notrans='homogenous_gamma_%.2f_linecap_0.40Q_nodes.npz' % gamma
        if (gval == 26): guess=[0.75,0.86]
        if (start == 50): guess=[0.8075,0.9075]
        if (last_quant[0] !=0.): guess=last_quant
        last_quant=find_balancing_reduction_quantiles(reduction=red,eps=1.e-4,guess=guess,stepsize=0.0025,copper=balmin,notrans=balmax,file_notrans=f_notrans,gamma=gamma,alpha=0.7,save_filename=save_filename)
        quantiles.append(np.array(last_quant).copy())
    qsave_file='homogenous_gamma'
    if (start != None):
        qsave_file += '_from_%.2f' % min(gvals)
    if (stop != None):
        qsave_file += '_to_%.2f' % max(gvals)
    qsave_file += '_balred_quantiles_refined'
    np.save('./results/'+qsave_file,quantiles)


def gamma_logfit(linecap='copper',step=2,start=None,stop=None):
    if (linecap == 'copper'):
        copper = 1
        h0 = None
    elif (linecap == '0.40Q'):
        h0 = get_quant(0.40)
        copper = 0
    elif (linecap == '0.90Q'):
        h0 = get_quant(0.90)
        copper = 0
    elif (linecap == 'today'):
        h0 = None
        copper = 0
    elif (linecap == '0.6_times_0.90Q'):
        h0 = 0.6*get_quant(0.90)
        copper = 0
    elif (linecap == '0.8_times_0.90Q'):
        h0 = 0.8*get_quant(0.90)
        copper = 0
    else:
        print 'invalid linecap. abort'
        return
    print 'linecap: ',linecap
        
    #generate_basepath_gamma_alpha(step=step)
    if start != None:
        skip = start
    else:
        skip = 0
    if stop != None:
        skip_end = stop
    else:
        skip_end=60
    years= arange(1990,2050+1,1)
    for year in years:
        if (year<1990+skip or year > 1990+skip_end):
            continue
        gammas=array(get_basepath_gamma(year,step=step))
        alphas=array(get_basepath_alpha(year,step=step))
        N = Nodes()
        N.set_alphas(alphas)
        N.set_gammas(gammas)
        print "Now calculating for year = ",year
        F = sdcpf(N,h0=h0,copper=copper)
        name = 'logfit_gamma_year_%u_linecap_%s_step_%u' % (year,linecap,step)
        print name
        N.save_nodes_small(name+'_nodes')
        np.save('./results/'+name+'_flows',F)
        del N


def gamma_logfit_balred(red=[0.50,0.90],path='./results/',guessfile='logistic_gamma_balred_quantiles',notrans_file='Bvsg_logfit_gamma_linecap_0.40Q',copper_file='Bvsg_logfit_gamma_linecap_copper',step=2,start=None,stop=None):
        
    #generate_basepath_gamma_alpha(step=step)
    if start != None:
        skip = start
    else:
        skip = 0
    if stop != None:
        skip_end = stop
    else:
        skip_end=60+1
    years= arange(1990+skip,1990+skip_end,1)
    quantiles = []
    guessfile += ('_step_%u' % step) + '.npy'
    guesses = np.load(path+guessfile)
    cfile = copper_file + ('_step_%u.npy' % step)
    balmins=np.load(path+cfile)[:,1]
    nfile = notrans_file + ('_step_%u.npy' % step)
    balmaxes = np.load(path+nfile)[:,1]
    last_quant = np.zeros(2)
    for year in years:
        gammas=array(get_basepath_gamma(year,step=step))
        alphas=array(get_basepath_alpha(year,step=step))
        print "Now calculating for year = ",year
        save_filename = []
        for i in range(len(red)):
            save_filename.append(('logfit_gamma_year_%u_balred_%.2f_step_%u' % (year,red[i],step)))
        f_notrans='logfit_gamma_year_%u_linecap_0.40Q_step_%u_nodes.npz' % (year,step)
        balmax=balmaxes[year-1990]
        balmin=balmins[year-1990]
        guess=guesses[year-1990,:]
        if (2014 <= year and year <= 2017): guess=[0.80,0.88]
        if (last_quant[0] !=0): guess=last_quant
        last_quant=find_balancing_reduction_quantiles(reduction=red,eps=1.e-4,guess=guess,stepsize=0.0025,copper=balmin,notrans=balmax,file_notrans=f_notrans,gamma=gammas,alpha=alphas,save_filename=save_filename,stepname=step)
        quantiles.append(np.array(last_quant).copy())
        print quantiles
    qsave_file='logistic_gamma'
    if (start != None):
        qsave_file += '_from_%u' % min(years)
    if (stop != None):
        qsave_file += '_to_%u' % max(years)
    qsave_file += '_balred_quantiles_refined'
    qsave_file += '_step_%u' % step
    np.save('./results/'+qsave_file,quantiles)


def gamma_logfit_balred_capped(path='./results/',step=2,start=None,stop=None):
        
    red=[0.50,0.90]
    #generate_basepath_gamma_alpha(step=step)
    if start != None:
        skip = start
    else:
        skip = 0
    if stop != None:
        skip_end = stop
    else:
        skip_end=60+1
    years= arange(1990+skip,1990+skip_end,1)

    quant_file = path + ('logistic_gamma_balred_quantiles_refined_step_%u.npy' % step)
    quantiles = np.load(quant_file)
    quant_end = quantiles[-1]
    copper_file = ''
    if (step == 2):
        copper_file = 'results/copper_flows.npy'
    else:
        copper_file = ('results/logfit_gamma_year_2050_linecap_copper_step_%u_flows.npy' % step)
    Nlinks = 44
    h_end = np.zeros((2*Nlinks,2))
    for i in range(2):
        h_end[:,i] = get_quant(quant_end[i],filename=copper_file)

    for year in years:
        for i in range(len(red)):
            save_filename = 'logfit_gamma_year_%u_balred_%.2f_step_%u_capped_investment' % (year,red[i],step)
            # if no overcapacities present: just link to old files
            if (quantiles[year-1990,i] <= quant_end[i]):
                link_flows=save_filename+'_flows.npy'
                link_nodes=save_filename+'_nodes.npz'
                target = 'logfit_gamma_year_%u_balred_%.2f_step_%u' % (year,red[i],step)
                target_flows = target+'_flows.npy'
                target_nodes = target+'_nodes.npz'
                os.symlink(target_nodes,path+link_nodes)
                os.symlink(target_flows,path+link_flows)
                continue
            # if overcapacities present: cap them and calculate with capped caps
            print "Now calculating for year = ",year,", reduction = ",red[i]
            quantiles[year-1990,i]=quant_end[i]
            N=Nodes()                
            gammas=array(get_basepath_gamma(year,step=step))
            alphas=array(get_basepath_alpha(year,step=step))
            N.set_alphas(alphas)
            N.set_gammas(gammas)
            F = sdcpf(N,h0=h_end[:,i],copper=0)
            N.save_nodes_small(save_filename+'_nodes')
            np.save('./results/'+save_filename+'_flows',F)
            del N
                

######################################################
########### Process the results ######################
######################################################

def get_balancing_vs_gamma(path='./results/',prefix='homogenous_gamma',linecap='copper'):
    filename='Bvsg_'+prefix
    if (linecap.find('balred') >= 0): filename += '_'+linecap + '.npy'
    else: filename += ('_linecap_%s.npy' % linecap)
    if os.path.exists(path+filename):
        Bvsg=np.load(path+filename)
        return Bvsg
    gammas=arange(0,1.01,0.01)
    Bvsg=np.zeros((len(gammas),2))
    j=0
    for gamma in gammas:
        Bvsg[j,0]=gamma
        load_filename=prefix+'_%.2f' % gamma
        if (linecap.find('balred') >= 0): load_filename += '_'+linecap
        else: load_filename += ('_linecap_%s' % linecap)
        load_filename += '_nodes.npz'
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
    np.save(path+filename,Bvsg)
    return Bvsg


def get_flows_vs_gamma(prefix='homogenous_gamma',path='./results/',linecap='copper'):
    filename='Fvsg_'+prefix
    if (linecap.find('balred') >= 0): filename += '_'+linecap + '.npy'
    else: filename += ('_linecap_%s.npy' % linecap)
    if os.path.exists(path+filename):
        Fvsg=np.load(path+filename)
        return Fvsg
    gammas=arange(0,1.01,0.01)
    Fvsg=np.zeros((len(gammas),3))
    j=0
    for gamma in gammas:
        Fvsg[j,0]=gamma
        load_filename=prefix+'_%.2f' % gamma
        if (linecap.find('balred') >= 0): load_filename += '_'+linecap
        else: load_filename += ('_linecap_%s' % linecap)
        load_filename += '_flows.npy'
        flows=np.load(path+load_filename)
        a=0.
        b=0.
        d=0.
        for i in flows:
            a+=sum(abs(i))
            b+=sum(i*i)
        print "Absolute flows %.4e, squared flows %.4e" % (a,b)
        Fvsg[j,1]=a
        Fvsg[j,2]=b
        j+=1
    np.save(path+filename,Fvsg)
    return Fvsg


def get_balancing_vs_year(prefix='logfit_gamma',linecap='copper',step=2,capped_inv=False):
    filename = prefix
    if (linecap.find('balred')>=0):
        filename += '_'+linecap
    else:
        filename += '_linecap_'+linecap
    filename +=('_step_%u' % step)
    if (linecap.find('balred')>=0):
        if (capped_inv): filename += '_capped_investment'
    print filename
    if os.path.exists('./results/Bvsg_'+filename+'.npy'):
        Bvsg=np.load('./results/Bvsg_'+filename+'.npy')
        return Bvsg
    years=arange(1990,2050+1,1)
    Bvsg=np.zeros((len(years),2))
    j=0
    for year in years:
        Bvsg[j,0]=year
        load_filename = prefix+('_year_%u' %year)
        if (linecap.find('balred')>=0):
            load_filename += '_'+linecap
        else:
            load_filename += '_linecap_'+linecap
        load_filename += ('_step_%u' % step)
        if (linecap.find('balred')>=0):
            if (capped_inv): load_filename += '_capped_investment'
        load_filename += '_nodes.npz'
        print load_filename
        N=Nodes(load_filename=load_filename)
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        del N
        c=a/d
        Bvsg[j,1]=c
        j+=1
    np.save('./results/Bvsg_'+filename,Bvsg)
    return Bvsg


def get_gamma_vs_year(path='./results/',prefix='logfit_gamma',step=2):
    filename = prefix+'_vs_year_step_%u.npy' % step
    if os.path.exists(path+filename):
        gamma_vs_year = np.load(path+filename)
        return gamma_vs_year
    N = Nodes()
    weight = []
    for i in N:
        weight.append(i.mean)
    del N
    gamma_vs_year = []
    years=arange(1990,2050+1,1)
    for year in years:
        gamma=array(get_basepath_gamma(year,step=step))
        counter = 0.
        denom = 0.
        for i in range(len(weight)):
            counter += weight[i]*gamma[i]
            denom += weight[i]
        gamma_vs_year.append(counter/denom)
        np.save(path+filename,gamma_vs_year)
    return gamma_vs_year


def get_linecaps_vs_gamma(path='./results/',prefix='homogenous_gamma_balred_quantiles_refined'):
    name = prefix + '.npy'
    quantiles=transpose(np.load(path+name))
    #print shape(quantiles)
    N=Nodes()
    km,kv,h0,lf=AtoKh(N) # h0: linecap today
    del N
    linecaps=np.zeros((2,100+1))
    for i in range(2):
         for j in range(100+1):
            if (quantiles[i,j] == 0):
                linecaps[i,j]=0
            else:
                if (step == 2):
                    copper_file = 'results/copper_flows.npy'
                else:
                    copper_file = ('results/logfit_gamma_year_2050_linecap_copper_step_%u_flows.npy' % step)
                h=get_quant(quantiles[i,j],filename=copper_file)
                inv=0.
                for l in range(len(h)):
                    inv += get_positive(h[l]-h0[l])
                linecaps[i,j]=inv
       # for j in range(100+1):
       #      if (quantiles[i,j] == 0):
       #          h=h0
       #      else:
       #          h=get_quant(quantiles[i,j])
       #      inv=0.
       #      for l in range(len(h)):
       #          inv += h[l]
       #      linecaps[i,j]=inv
    return linecaps


def get_alpha_vs_year(path='./results/',prefix='logfit_gamma',step=2):
    filename = prefix+'_alpha_vs_year_step_%u.npy' % step
    if os.path.exists(path+filename):
        alpha_vs_year = np.load(path+filename)
        return alpha_vs_year
    N = Nodes()
    weight = []
    for i in N:
        weight.append(i.mean)
    del N
    alpha_vs_year = []
    years=arange(1990,2050+1,1)
    for year in years:
        alpha=array(get_basepath_alpha(year,step=step))
        counter = 0.
        denom = 0.
        for i in range(len(weight)):
            counter += weight[i]*alpha[i]
            denom += weight[i]
        alpha_vs_year.append(counter/denom)
        np.save(path+filename,alpha_vs_year)
    return alpha_vs_year


def get_flows_vs_year(prefix='logfit_gamma',path='./results/',linecap='copper',step=2,capped_inv=False):
    filename=path+'Fvsg_'+prefix
    if (linecap.find('balred')>=0):
        filename += '_'+linecap
    else:
        filename += '_linecap_'+linecap
    filename += '_step_%u' % step
    if (linecap.find('balred')>=0):
        if (capped_inv): filename += '_capped_investment'
    if os.path.exists(filename+'.npy'):
        Fvsg=np.load(filename+'.npy')
        return Fvsg
    years=arange(1990,2050+1,1)
    Fvsg=np.zeros((len(years),3))
    j=0
    for year in years:
        Fvsg[j,0]=year
        load_filename = prefix+('_year_%u' %year)
        if (linecap.find('balred')>=0):
            load_filename += '_'+linecap
        else:
            load_filename += '_linecap_'+linecap
        load_filename += ('_step_%u' % step)
        if (linecap.find('balred')>=0):
            if (capped_inv): load_filename += '_capped_investment'
        load_filename += '_flows.npy'
        print load_filename
        flows=np.load(path+load_filename)
        a=0.
        b=0.
        for i in flows:
            a+=sum(abs(i))
            b+=sum(i*i)
        print "Absolute flows %.4e, squared flows %.4e" % (a,b)
        Fvsg[j,1]=a
        Fvsg[j,2]=b
        j+=1
    np.save(filename,Fvsg)
    return Fvsg


def get_optimal_alphas(txtfile='../DataAndPredictionsGammaAlpha/gamma.csv',step=2):
    data = np.genfromtxt(txtfile,delimiter=',',skip_header=0)

    ISO = ['AT', 'BE', 'BG', 'BA', 'CZ', 'CH', 'DE', 'DK', 'ES', 
    'FR', 'FI', 'GB', 'GR', 'HU', 'IT', 'IE', 'HR', 'LU', 'NO', 
    'NL', 'PT', 'PL', 'RO', 'SE', 'SK', 'SI', 'RS']

    alpha = []    
    for i in arange(1,data.shape[0]-1,2): # loop over countries
        wind = array(data[i][2:-4])
        solar = array(data[i+1][2:-4])
        if (step == 3):
            wind[-1] = data[i][-4] # take alternative values for 2050 in step 3
            solar[-1] = data[i+1][-4]
        elif (step == 4):
            wind[-1] = data[i][-3] # take alternative values for 2050 in step 4
            solar[-1] = data[i+1][-3]
        gamma = solar[-1] + wind[-1]
        t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO[(i-1)/2])    
        alpha_opt = get_optimal_mix_balancing(L,Gw,Gs,gamma,CS=None)
        alpha.append(alpha_opt)
        print ISO[(i-1)/2], alpha_opt
    return alpha


def get_linecaps_vs_year(path='./results/',prefix='logistic_gamma_balred_quantiles',step=2,capped_inv=False):
    name = prefix + ('_step_%u' % step) + '.npy'
    quantiles=transpose(np.load(path+name))
    if (capped_inv):
        quant_end = quantiles[:,-1]
        for i in range(2):
            for j in range(60+1):
                if quantiles[i,j] >= quant_end[i]:
                    quantiles[i,j] = quant_end[i]
    #print shape(quantiles)
    N=Nodes()
    km,kv,h0,lf=AtoKh(N) # h0: linecap today
    del N
    linecaps=np.zeros((2,60+1))
    for i in range(2):
        for j in range(60+1):
            if (quantiles[i,j] == 0):
                linecaps[i,j]=0
            else:
                if (step != 2):
                    copper_flow_file=('results/logfit_gamma_year_2050_linecap_copper_step_%u_flows.npy' % step)
                else:
                    copper_flow_file='results/copper_flows.npy'
                h=get_quant(quantiles[i,j],filename=copper_flow_file)
                inv=0.
                for l in range(len(h)):
                    inv += get_positive(h[l]-h0[l])
                linecaps[i,j]=inv
            # if (quantiles[i,j] == 0):
            #     h=h0
            # else:
            #     if (step == 3):
            #         copper_flow_file='results/logfit_gamma_year_2050_linecap_copper_step_3_flows.npy'
            #     else:
            #         copper_flow_file='results/copper_flows.npy'
            #     h=get_quant(quantiles[i,j],filename=copper_flow_file)
            # inv=0.
            # for l in range(len(h)):
            #     inv += h[l]
            # linecaps[i,j]=inv
    return linecaps


def get_balancing_and_deficit(path='./results/',datpath='./data/',step=2,linecap='copper'):
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']
        
    years = arange(1990,2050+1,1)
    balancing = np.zeros((len(years),len(ISO)))
    deficit = np.zeros((len(years),len(ISO)))
    N_ref = Nodes()
    gamma_vs_year= array(get_gamma_vs_year(step=step))
    for year in years:
        fname = 'logfit_gamma_year_%u' % year
        if ('balred' in linecap): fname += '_'+linecap
        else: fname += '_linecap_'+linecap
        fname += '_step_%u' % step
        if ('balred' in linecap):
            fname += '_capped_investment'
        fname += '_nodes.npz'
        N = Nodes(load_filename=fname)
        gamma = np.zeros(len(N))
        alpha = np.zeros(len(N))
        for n in N:
            load = n.mean*n.nhours
            # count only excess balancing
            bal = sum(n.balancing)/load - (1.-gamma_vs_year[year-1990])
            balancing[year-1990,n.id] = bal
            alpha[n.id] = n.alpha
            gamma[n.id] = n.gamma
        del N
        N_ref.set_alphas(alpha)
        N_ref.set_gammas(gamma)
        for n in N_ref:
            load = n.mean*n.nhours
            # for gamma < 1: remember that part of neg. mismatch 
            # is still covered conventionally (only "excess deficit" matters)
            dfc = sum(get_positive(-n.mismatch))/load - (1.-gamma_vs_year[year-1990])
            deficit[year-1990,n.id] = dfc

    lc = linecap
    if (lc.find('balred') == -1):
        lc = 'linecap_'+lc
    i=0
    for iso in ISO:
        outfile = 'balancing_and_deficit_%s_step_%u_%s.npy' % (lc,step,iso)
        data=np.array([np.array(years),np.array(balancing[:,i]),np.array(deficit[:,i])])
        np.save(path+outfile,data)
        i += 1
    return


def get_curtailment_and_excess(path='./results/',datpath='./data/',step=2,linecap='copper'):
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']
        
    years = arange(1990,2050+1,1)
    curtailment = np.zeros((len(years),len(ISO)))
    excess = np.zeros((len(years),len(ISO)))
    N_ref = Nodes()
    for year in years:
        fname = 'logfit_gamma_year_%u' % year
        if ('balred' in linecap): fname += '_'+linecap
        else: fname += '_linecap_'+linecap
        fname += '_step_%u' % step
        if ('balred' in linecap):
            fname += '_capped_investment'
        fname += '_nodes.npz'
        N = Nodes(load_filename=fname)
        gamma = np.zeros(len(N))
        alpha = np.zeros(len(N))
        for n in N:
            load = n.mean*n.nhours
            curt = sum(n.curtailment)/load
            curtailment[year-1990,n.id] = curt
            alpha[n.id] = n.alpha
            gamma[n.id] = n.gamma
        del N
        N_ref.set_alphas(alpha)
        N_ref.set_gammas(gamma)
        for n in N_ref:
            load = n.mean*n.nhours
            exc = sum(get_positive(n.mismatch))/load
            excess[year-1990,n.id] = exc

    lc = linecap
    if (lc.find('balred') == -1):
        lc = 'linecap_'+lc
    i=0
    for iso in ISO:
        outfile = 'curtailment_and_excess_%s_step_%u_%s.npy' % (lc,step,iso)
        data=np.array([np.array(years),np.array(curtailment[:,i]),np.array(excess[:,i])])
        np.save(path+outfile,data)
        i += 1
    return


def get_balancing_quantiles_vs_year(path='./results/',datpath='./data/',step=2,linecap='copper'):
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']
        
    years = arange(1990,2050+1,1)
    quant=[0.9,0.99,1.]
    bal_quant = np.zeros((len(years),len(ISO),len(quant)))
    for year in years:
        fname = 'logfit_gamma_year_%u' % year
        if ('balred' in linecap): fname += '_'+linecap
        else: fname += '_linecap_'+linecap
        fname += '_step_%u' % step
        if ('balred' in linecap):
            fname += '_capped_investment'
        fname += '_nodes.npz'
        N = Nodes(load_filename=fname)
        for n in N:
            load=n.mean*n.nhours
            bal = n.balancing/load
            hst = hist(bal,cumulative=True,bins=500,normed=True)
            for i in range(len(quant)):
                for j in range(len(hst[0])):
                    if hst[0][j]>=quant[i]:
                        bal_quant[year-1990,n.id,i]=hst[1][j]
                        break
                if j == len(hst[0]) - 1:
                    bal_quant[year-1990,n.id,i]=max(hst[1][:])
            clf()
        del N

    lc = linecap
    if (not ('balred' in lc)):
        lc = 'linecap_'+lc
    i=0
    for iso in ISO:
        outfile = 'balancing_quantiles_%s_step_%u_%s.npy' % (lc,step,iso)
        np.save(path+outfile,bal_quant[:,i,:])
        i += 1
    return


def get_excess_quantiles_vs_year(path='./results/',datpath='./data/',step=2,linecap='copper'):
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']
        
    years = arange(1990,2050+1,1)
    quant=[0.9,0.99,1.]
    excess_quant = np.zeros((len(years),len(ISO),len(quant)))
    N = Nodes()
    for year in years:
        gammas=array(get_basepath_gamma(year,step=step))
        alphas=array(get_basepath_alpha(year,step=step))
        N.set_alphas(alphas)
        N.set_gammas(gammas)
        for n in N:
            load=n.mean*n.nhours
            excess = get_positive(n.mismatch)/load
            hst = hist(excess,cumulative=True,bins=500,normed=True)
            for i in range(len(quant)):
                for j in range(len(hst[0])):
                    if hst[0][j]>=quant[i]:
                        excess_quant[year-1990,n.id,i]=hst[1][j]
                        break
                if j == len(hst[0]) - 1:
                    excess_quant[year-1990,n.id,i]=max(hst[1][:])
            clf()
    del N

    lc = linecap
    if (not ('balred' in lc)):
        lc = 'linecap_'+lc
    i=0
    for iso in ISO:
        outfile = 'excess_quantiles_%s_step_%u_%s.npy' % (lc,step,iso)
        np.save(path+outfile,excess_quant[:,i,:])
        i += 1
    return


def get_deficit_quantiles_vs_year(path='./results/',datpath='./data/',step=2,linecap='copper'):
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']
        
    years = arange(1990,2050+1,1)
    quant=[0.9,0.99,1.]
    deficit_quant = np.zeros((len(years),len(ISO),len(quant)))
    N = Nodes()
    for year in years:
        gammas=array(get_basepath_gamma(year,step=step))
        alphas=array(get_basepath_alpha(year,step=step))
        N.set_alphas(alphas)
        N.set_gammas(gammas)
        for n in N:
            load=n.mean*n.nhours
            deficit = get_positive(-n.mismatch)/load
            hst = hist(deficit,cumulative=True,bins=500,normed=True)
            for i in range(len(quant)):
                for j in range(len(hst[0])):
                    if hst[0][j]>=quant[i]:
                        deficit_quant[year-1990,n.id,i]=hst[1][j]
                        break
                if j == len(hst[0]) - 1:
                    deficit_quant[year-1990,n.id,i]=max(hst[1][:])
            clf()
    del N

    lc = linecap
    if (not ('balred' in lc)):
        lc = 'linecap_'+lc
    i=0
    for iso in ISO:
        outfile = 'deficit_quantiles_%s_step_%u_%s.npy' % (lc,step,iso)
        np.save(path+outfile,deficit_quant[:,i,:])
        i += 1
    return


######################################################
########### Plot the results #########################
######################################################

# some nice color palettes
# cl = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] # Ocean Five (COLOURlovers) *
# cl = ['#CB9FBE','#3D6BA5','#99CD26','#270323','#9A8224'] # my rainbow :-) (CL)
# cl = ['#8378B0','#FFC352','#B94F11','#001119','#B49382'] # black & tan (CL)
# cl = ['#1BA9D0','#8F5A42','#FFAD7A','#36291F','#E6304B'] # Click (CL) *
# cl = ['#8F07BC','#4007BC','#F48709','#0C9D4C','#9D0C2B'] # lynnn (CL) *
# cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL) *

def plot_balancing_vs_gamma(path='./results/',prefix='homogenous_gamma',linecaps=['0.40Q','today','balred_0.50','balred_0.90','copper'],label=['no transmission','line capacities as of today',r'50$\,$% bal. reduction quantile line capacities',r'90$\,$% bal. reduction quantile line capacities','copper plate'],title_=r"homogenous increase in $\gamma\,$; $ \alpha_{\rm W}=0.7 $",picname='balancing_vs_homogenous_gamma.png'):
    fig=figure(1); clf()
    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    pp = []
    i=0
    for cap in linecaps:
        print 'Processing '+cap
        data = get_balancing_vs_gamma(linecap=cap)
        pp_x = array(data[:,0])
        pp_y = array(data[:,1]) # balancing
        pp_y = pp_y - (1-pp_x)  # excess balancing
        pp_ = plot(pp_x,pp_y,lw=1.5,color=cl[i])
        pp.extend(pp_)
        i += 1

    fig.suptitle(title_,fontsize=14,y=0.96)
    axis(xmin=0,xmax=1,ymin=0,ymax=0.3)
    xlabel('share $\gamma$ of VRES in total electricity production')
    ylabel(r'excess balancing/av.h.l. (bal.$ - (1-\gamma)$)')

    pp_label = label
    leg = legend(pp,pp_label,loc='upper left');
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    save_figure(picname)


def plot_flows_vs_gamma(path='./results/',prefix='homogenous_gamma',linecaps=['0.40Q','today','balred_0.50','balred_0.90','copper'],title_=r"homogenous increase in $\gamma\,$; $ \alpha_{\rm W}=0.7 $",label=['no transmission','line capacities as of today',r'50$\,$% bal. reduction quantile line capacities',r'90$\,$% bal. reduction quantile line capacities','copper plate'],picname='flows_vs_homogenous_gamma.png'):
    fig=figure(1); clf()
    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    pp = []
    i=0
    for cap in linecaps:
        fname = 'Fvsg_'+prefix
        if cap.find('balred')>=0:
            fname += '_'+cap+'.npy'
        else:
            fname += '_linecap_'+cap+'.npy'
        print 'Processing '+fname
        if not os.path.exists(path+fname):
            data = get_flows_vs_gamma(linecap=cap)
        else:
            data = np.load(path+fname)
        pp_x = array(data[:,0])
        pp_y = array(data[:,1])/1.e9
        pp_ = plot(pp_x,pp_y,lw=1.5,color=cl[i])
        pp.extend(pp_)
        i += 1

    fig.suptitle(title_,fontsize=14,y=0.96)
    axis(xmin=0,xmax=1,ymin=-0.05,ymax=5)
    xlabel('share $\gamma$ of VRES in total electricity production')
    ylabel(r'sum of absolute flows/$10^9\,$MW')

    pp_label = label
    leg = legend(pp,pp_label,loc='upper left');
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    save_figure(picname)
    

def plot_linecaps_vs_gamma(path='./results/',prefix='homogenous_gamma_balred_quantiles_refined',title_=r"homogenous increase in $\gamma\,$; $ \alpha_{\rm W}=0.7 $",label=['needed for 50 % balancing reduction','needed for 90 % balancing reduction']):
    linecaps=get_linecaps_vs_gamma(path=path,prefix=prefix)
    fig=figure(1); clf()
    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    pp = []
    pp_x=arange(0.,1.01,0.01)
    for i in range(2):
        pp_y=linecaps[i,:]*1e-5
        print shape(pp_y)
        pp_=plot(pp_x,pp_y,lw=1.5,color=cl[i+2])
        pp.extend(pp_)
    pp_label=label
    leg=legend(pp,pp_label,loc='upper left')
    axis(xmin=0,xmax=1,ymin=0,ymax=10.)
    xlabel('gamma')
    ylabel(r'total new line capacities/$10^5\,$MW')

    fig.suptitle(title_,fontsize=14,y=0.96)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    picname = 'linecaps_vs_homogenous_gamma' + '.png'
    save_figure(picname)


def plot_investment_vs_gamma(path='./results/',prefix='homogenous_gamma_balred_quantiles_refined',title_=r"homogenous increase in $\gamma\,$; $ \alpha_{\rm W}=0.7 $",label=['needed for 50 % balancing reduction','needed for 90 % balancing reduction'],title_leg=r"2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $"):
    linecaps=get_linecaps_vs_gamma(path=path,prefix=prefix)
    fig=figure(1); clf()
    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    pp = []
    pp_x=arange(0.,1.01,0.01)
    for i in range(2):
        caps=linecaps[i,:]*1e-5
        # get positive part of derivative
        pp_y=np.zeros(len(caps))
        for j in range(len(pp_y)):
            if j==0: continue
            pp_y[j] = get_positive(caps[j]-max(caps[:j]))
        pp_=plot(pp_x,pp_y,lw=1.5,color=cl[i+2])
        pp.extend(pp_)
    pp_label=label
    leg=legend(pp,pp_label,title=title_leg,loc='upper left')

    fig.suptitle(title_,fontsize=14,y=0.96)
    axis(xmin=0,xmax=1,ymin=0,ymax=2.)
    xlabel('gamma')
    ylabel(r'necessary total new line capacities/$10^5\,$MW per gamma')

    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    picname = 'investment_vs_homogenous_gamma' + '.png'
    save_figure(picname)


def plot_balancing_vs_year(path='./results/',prefix='logfit_gamma',linecaps = ['0.40Q','today','balred_0.50','balred_0.90','copper'],label=['no transmission','line capacities as of today',r'50$\,$% bal. reduction quantile line capacities',r'90$\,$% bal. reduction quantile line capacities','copper plate'],picname='balancing_vs_year',step=2,capped_inv=False):

    fig=figure(1); clf()
    #calculate average gamma as weighted mean of all countries for all years
    gamma_vs_year= array(get_gamma_vs_year(step=step))
    gamma_hom=arange(0.,1.01,0.01)
    fig.subplots_adjust(top=0.8) # leave more top margin for the extra ticks
    host = host_subplot(111)

    par = host.twin()
    axis(xmin=1990,xmax=2050,ymin=0.,ymax=0.3)

    host.set_xlabel('year')
    host.set_ylabel(r'excess balancing/av.h.l. (bal.$ - (1-\gamma)$)')
    par.set_xlabel(r'effective $\gamma$')
    gamma_ticks = arange(1990,2050+1,5)
    gamma_labels = ['%.3f' % gamma_vs_year[i-1990] for i in gamma_ticks]
    # for i in range(len(gamma_ticks)):
    #     print gamma_ticks[i], gamma_labels[i]
    par.set_xticks(gamma_ticks)
    par.set_xticklabels(gamma_labels,rotation=30)#,shift=3)

    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    pp = []
    i=0
    for linecap in linecaps:
        if (linecap.find('balred')>=0):
            name=prefix+'_'+linecap
        else:
            name = prefix+'_linecap_'+linecap
        name += ('_step_%u' % step)
        if (linecap.find('balred')>=0):
            if (capped_inv):
                name += '_capped_investment'
        fname = 'Bvsg_'+name+'.npy'
        print 'Processing '+fname
        if not os.path.exists(path+fname):
            data = get_balancing_vs_year(prefix=prefix,step=step,linecap=linecap,capped_inv=capped_inv)
        else:
            data = np.load(path+fname)
        pp_x = array(data[:,0])
        pp_y = array(data[:,1]) # balancing
        # excess balancing is not so easy here
        pp_y = pp_y - (1.-gamma_vs_year)
        pp_ = host.plot(pp_x,pp_y,lw=1.5,color=cl[i])
        pp.extend(pp_)
        # plot homogenous balancing for comparison
        bal_hom=get_balancing_vs_gamma(linecap=linecap)
        bal_hom=bal_hom[:,1]
        bal_hom_eff=np.interp(gamma_vs_year,gamma_hom,bal_hom)
        pp_y = bal_hom_eff - (1.-gamma_vs_year)
        pp_ = host.plot(pp_x,pp_y,lw=1.5,color=cl[i],ls='--')
        i += 1

    title_ = ''
    if capped_inv: title_ = 'Capped transmission investment'
    else: title_ = 'Full transmission investment'
    fig.suptitle(title_,fontsize=14,y=0.96)

    title_leg=''
    if (step == 2): title_leg=r"2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $"
    if (step == 3): title_leg=r"2050 target:  76% VRES, $\alpha_{\rm W}=0.7 $"
    pp_label = label
    leg = legend(pp,pp_label,loc='upper left',title=title_leg);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()
    text(1991,0.03,'dashed lines indicate corresponding homogenous',fontsize='small',fontstyle='italic')
    text(1991,0.015,'gamma-growth balancing',fontsize='small',fontstyle='italic')

    picname = picname +'_'+ prefix + ('_step_%u' %step)
    if capped_inv:
        picname += '_capped_investment'
    picname += '.png'
    save_figure(picname)


def plot_flows_vs_year(path='./results/',prefix='logfit_gamma',linecaps = ['0.40Q','today','balred_0.50','balred_0.90','copper'],label=['no transmission','line capacities as of today',r'50$\,$% bal. reduction quantile line capacities',r'90$\,$% bal. reduction quantile line capacities','copper plate'],picname='flows_vs_year',step=2,capped_inv=False):
    fig=figure(1); clf()
    #calculate average gamma as weighted mean of all countries for all years
    gamma_vs_year= array(get_gamma_vs_year(step=step))
    gamma_hom=arange(0.,1.01,0.01)
    fig.subplots_adjust(top=0.8) # leave more top margin for the extra ticks
    host = host_subplot(111)

    par = host.twin()
    axis(xmin=1990,xmax=2050,ymin=0.,ymax=0.3)

    host.set_xlabel('year')
    host.set_ylabel(r'excess balancing/av.h.l. (bal.$ - (1-\gamma)$)')
    par.set_xlabel(r'effective $\gamma$')
    gamma_ticks = arange(1990,2050+1,5)
    gamma_labels = ['%.3f' % gamma_vs_year[i-1990] for i in gamma_ticks]
    # for i in range(len(gamma_ticks)):
    #     print gamma_ticks[i], gamma_labels[i]
    par.set_xticks(gamma_ticks)
    par.set_xticklabels(gamma_labels,rotation=30)#,shift=3)
    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    pp = []
    i=0
    #calculate average gamma as weighted mean of all countries for all years
    for linecap in linecaps:
        if (linecap.find('balred')>=0):
            name=prefix+'_'+linecap
        else:
            name = prefix+'_linecap_'+linecap
        name += ('_step_%u' % step)
        if (linecap.find('balred')>=0):
            if (capped_inv):
                name += '_capped_investment'
        fname = 'Fvsg_'+name+'.npy'
        print 'Processing '+fname
        if not os.path.exists(path+fname):
            data = get_flows_vs_year(prefix=prefix,step=step,linecap=linecap)
        else:
            data = np.load(path+fname)
        pp_x = array(data[:,0])
        pp_y = array(data[:,1])/1.e9
        pp_ = plot(pp_x,pp_y,lw=1.5,color=cl[i])
        pp.extend(pp_)
        fl_hom=get_flows_vs_gamma(linecap=linecap)
        fl_hom=fl_hom[:,1]
        fl_hom_eff=np.interp(gamma_vs_year,gamma_hom,fl_hom)
        pp_y = fl_hom_eff/1.e9
        pp_ = host.plot(pp_x,pp_y,lw=1.5,color=cl[i],ls='--')
        i += 1

    title_ = ''
    if capped_inv: title_ = 'Capped transmission investment'
    else: title_ = 'Full transmission investment'
    fig.suptitle(title_,fontsize=14,y=0.96)

    title_leg=''
    if (step == 2): title_leg=r"2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $"
    if (step == 3): title_leg=r"2050 target:  76% VRES, $\alpha_{\rm W}=0.7 $"
    pp_label = label
    leg = legend(pp,pp_label,loc='upper left',title=title_leg);
    axis(xmin=1990,xmax=2050,ymin=-0.05,ymax=5.)
    xlabel('year')
    ylabel(r'sum of absolute flows/$10^9\,$MW')

    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()
    text(1991,0.4,'dashed lines indicate corresponding homogenous',fontsize='small',fontstyle='italic')
    text(1991,0.15,'gamma-growth flows',fontsize='small',fontstyle='italic')

    picname = picname +'_'+ prefix + ('_step_%u' %step)
    if capped_inv: picname += '_capped_investment'
    picname += '.png'
    save_figure(picname)


def plot_linecaps_vs_year(path='./results/',prefix='logistic_gamma_balred_quantiles_refined',step=2,label=['needed for 50 % balancing reduction','needed for 90 % balancing reduction'],capped_inv=False):
    linecaps=get_linecaps_vs_year(path=path,prefix=prefix,step=step,capped_inv=capped_inv)
    fig=figure(1); clf()
    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    pp = []
    pp_x=arange(1990,2050+1,1)
    for i in range(2):
        pp_y=linecaps[i,:]*1e-5
        pp_=plot(pp_x,pp_y,lw=1.5,color=cl[i+2])
        pp.extend(pp_)
    pp_label=label
    title_leg=''
    if (step == 2): title_leg=r"2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $"
    if (step == 3): title_leg=r"2050 target:  76% VRES, $\alpha_{\rm W}=0.7 $"
    leg=legend(pp,pp_label,title=title_leg,loc='upper left')
    axis(xmin=1990,xmax=2050,ymin=0,ymax=10.)
    xlabel('year')
    ylabel(r'necessary total new line capacities/$10^5\,$MW')

    title_ = ''
    if capped_inv: title_ = 'Capped transmission investment'
    else: title_ = 'Full transmission investment'
    fig.suptitle(title_,fontsize=14,y=0.92)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    picname = 'linecaps_vs_year' + ('_step_%u' %step)
    if capped_inv: picname += '_capped_investment'
    picname += '.png'
    save_figure(picname)


def plot_investment_vs_year(path='./results/',prefix='logistic_gamma_balred_quantiles_refined',step=2,label=['needed for 50 % balancing reduction','needed for 90 % balancing reduction'],capped_inv=False):
    linecaps=get_linecaps_vs_year(path=path,prefix=prefix,step=step,capped_inv=capped_inv)
    fig=figure(1); clf()
    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    pp = []
    pp_x=arange(1990,2050+1,1)
    for i in range(2):
        caps=linecaps[i,:]*1e-5
        # get positive part of derivative
        pp_y=np.zeros(len(caps))
        for j in range(len(pp_y)):
            if j==0: continue
            pp_y[j] = get_positive(caps[j]-max(caps[:j]))
        pp_=plot(pp_x,pp_y,lw=1.5,color=cl[i+2])
        pp.extend(pp_)
    pp_label=label
    title_leg=''
    if (step == 2): title_leg=r"2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $"
    if (step == 3): title_leg=r"2050 target:  76% VRES, $\alpha_{\rm W}=0.7 $"
    leg=legend(pp,pp_label,title=title_leg,loc='upper left')
    axis(xmin=1990,xmax=2050,ymin=0,ymax=2.)
    xlabel('year')
    ylabel(r'necessary total new line capacities/$10^5\,$MW per year')

    title_ = ''
    if capped_inv: title_ = 'Capped transmission investment'
    else: title_ = 'Full transmission investment'
    fig.suptitle(title_,fontsize=14,y=0.92)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    picname = 'investment_vs_year' + ('_step_%u' %step)
    if capped_inv: picname += '_capped_investment'
    picname += '.png'
    save_figure(picname)


def plot_export_and_curtailment(path='./results/',step=2,title_='Export opportunities for '):
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']

    linecap = ['copper','balred_0.90','balred_0.50','today']
    lc = []
    for lica in linecap:
        if ('balred' in lica): lc.append(lica)
        else: lc.append('linecap_'+lica)

    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    for iso in ISO:
        fig=figure(1); clf()
        ax=fig.add_subplot(1,1,1)
        i=0
        pp_label = ['overproduction','export with copper plate transmission','export with 90% balancing reduction transmission','export with 50% balancing reduction transmission','export with transmission lines as of today']
        for lica in lc:
            datfile = 'curtailment_and_excess_%s_step_%u_%s.npy' % (lica,step,iso)
            data = np.load(path+datfile)
            years = data[0,:]
            curtailment = data[1,:]
            excess = data[2,:]
            pp_x=years
            if (i == 0):
                pp_y=excess
                ax.bar(pp_x-0.35,pp_y,color=cl[0],label=pp_label[0])
            pp_y=excess-curtailment # export
            ax.bar(pp_x-0.35,pp_y,color=cl[i+1],label=pp_label[i+1])
            i += 1

        if (step == 2): title_leg = r'2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $'
        if (step == 3): title_leg = r'2050 target:  76% VRES, $\alpha_{\rm W}=0.7 $'
        leg=ax.legend(loc='upper left',title=title_leg)
        ltext  = leg.get_texts();
        setp(ltext, fontsize='small')    # the legend text fontsize

        axis(xmin=1990,xmax=2050.5,ymin=0,ymax=0.3)
        xlabel('year')
        ylabel(r'production/av.l.h.')
        tlte = title_
        if (iso=='NL' or iso=='CZ'): tlte += 'the ' # get the grammar right
        tlte += ISO2name(ISO=iso)
        fig.suptitle(tlte,fontsize=14,y=0.96)

        picname = 'excess_and_export_vs_year'+ ('_step_%u' %step) + ('_%s.png' % iso)
        save_figure(picname)


def plot_import_and_deficit(path='./results/',step=2,title_='Import opportunities for '):
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']

    linecap = ['copper','balred_0.90','balred_0.50','today']
    lc = []
    for lica in linecap:
        if ('balred' in lica): lc.append(lica)
        else: lc.append('linecap_'+lica)

    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    for iso in ISO:
        fig=figure(1); clf()
        ax=fig.add_subplot(1,1,1)
        i=0
        pp_label = ['deficit','import with copper plate transmission','import with 90% balancing reduction transmission','import with 50% balancing reduction transmission','import with transmission lines as of today']
        for lica in lc:
            datfile = 'balancing_and_deficit_%s_step_%u_%s.npy' % (lica,step,iso)
            data = np.load(path+datfile)
            years = data[0,:]
            balancing = data[1,:]
            deficit = data[2,:]
            pp_x=years
            if (i == 0):
                pp_y=deficit
                ax.bar(pp_x-0.35,pp_y,color=cl[0],label=pp_label[0])
            pp_y=deficit-balancing # import
            ax.bar(pp_x-0.35,pp_y,color=cl[i+1],label=pp_label[i+1])
            i += 1

        if (step == 2): title_leg = r'2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $'
        if (step == 3): title_leg = r'2050 target:  76% VRES, $\alpha_{\rm W}=0.7 $'
        leg=ax.legend(loc='upper left',title=title_leg)
        ltext  = leg.get_texts();
        setp(ltext, fontsize='small')    # the legend text fontsize

        axis(xmin=1990,xmax=2050.5,ymin=0,ymax=.3)
        xlabel('year')
        ylabel(r'deficit/av.l.h.')
        tlte = title_
        if (iso=='NL' or iso=='CZ'): tlte += 'the ' # get the grammar right
        tlte += ISO2name(ISO=iso)
        fig.suptitle(tlte,fontsize=14,y=0.96)

        picname = 'deficit_and_import_vs_year'+ ('_step_%u' %step) + ('_%s.png' % iso)
        save_figure(picname)


def plot_quantiles_vs_year(path='./results/',qtype='balancing',quant=0.9,step=2):
    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']

    years=arange(1990,2050+1,1)
    linecap=['linecap_0.40Q','linecap_today','balred_0.50','balred_0.90','linecap_copper']
    lbl=['no transmission','line capacities as of today',r'50$\,$% bal. reduction quantile line capacities',r'90$\,$% bal. reduction quantile line capacities','copper plate']

    nhours = 70128
    for iso in ISO:
        fig=figure(1); clf()
        ax=subplot(1,1,1)
        pp_x=years
        i=0
        for lc in linecap:
            filename='%s_quantiles_%s_step_%u_%s.npy' % (qtype,lc,step,iso)
            quantiles=np.load(path+filename) # (years,quant) array
            if (quant==0.9):
                quantile=nhours*quantiles[:,0]
            elif (quant==0.99):
                quantile=nhours*quantiles[:,1]
            elif (quant==1.):
                quantile=nhours*quantiles[:,2]
            else:
                print ('quant %.2f not available!' % quant)
                return
            ax.bar(pp_x-0.35,quantile,color=cl[i],label=lbl[i])
            i += 1
        if (step == 2): title_leg = r'2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $'
        if (step == 3): title_leg = r'2050 target:  76% VRES, $\alpha_{\rm W}=0.7 $'
        leg=ax.legend(loc='upper left',title=title_leg)
        ltext  = leg.get_texts();
        setp(ltext, fontsize='small')    # the legend text fontsize

        axis(xmin=1990,xmax=2050.5,ymin=0,ymax=2.5)
        xlabel('year')
        ylabel(('%s quantiles/av.l.h.' % qtype))
        title_ = ('%u%% %s quantiles for ' % (100.*quant,qtype))+ ISO2name(ISO=iso)
        fig.suptitle(title_,fontsize=14,y=0.96)

        picname = ('%u_percent_%s_quantiles_vs_year' %(100*quant,qtype))+ ('_step_%u' %step)+ ('_%s.png' % iso)
        save_figure(picname)


def plot_cumulative_quantiles_vs_year(path='./results/',qtype='balancing',quant=0.9,step=2):
    cl = ['#490A3D','#BD1550','#E97F02','#F8CA00','#8A9B0F'] # by sugar (CL)
    ISO = ['AT', 'FI', 'NL', 'BA', 'FR', 'NO', 'BE', 'GB', 'PL', 'BG', 'GR', 'PT',
           'CH', 'HR', 'RO', 'CZ', 'HU', 'RS', 'DE', 'IE', 'SE', 'DK', 'IT', 'SI',
           'ES', 'LU', 'SK']

    years=arange(1990,2050+1,1)

    linecap=['linecap_0.40Q','linecap_today','balred_0.50','balred_0.90','linecap_copper']
    lbl=['no transmission','line capacities as of today',r'50$\,$% bal. reduction quantile line capacities',r'90$\,$% bal. reduction quantile line capacities','copper plate']

    N=Nodes()
    weight=[]
    nhours=0
    for n in N:
        weight.append(n.mean*n.nhours)
    nhours=N[0].nhours
    del N

    fig=figure(1); clf()
    ax=subplot(1,1,1)
    pp_x=years
    i=0
    for lc in linecap:
        quantile=np.zeros(len(years))
        j=0
        for iso in ISO:
            filename='%s_quantiles_%s_step_%u_%s.npy' % (qtype,lc,step,iso)
            quantiles=np.load(path+filename) # (years,quant) array
            if (quant==0.9):
                quantile += weight[j]*nhours*np.array(quantiles[:,0])
            elif (quant==0.99):
                quantile += weight[j]*nhours*np.array(quantiles[:,1])
            elif (quant==1.):
                quantile += weight[j]*nhours*np.array(quantiles[:,2])
            else:
                print ('quant %.2f not available!' % quant)
                return
            j += 1
        quantile /= sum(weight)
        ax.bar(pp_x-0.35,quantile,color=cl[i],label=lbl[i])
        i += 1
    if (step == 2): title_leg = r'2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $'
    if (step == 3): title_leg = r'2050 target:  76% VRES, $\alpha_{\rm W}=0.7 $'
    leg=ax.legend(loc='upper left',title=title_leg)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    axis(xmin=1990,xmax=2050.5,ymin=0,ymax=2.5)
    xlabel('year')
    ylabel(('%s quantiles/av.l.h.' % qtype))
    title_ = ('%u%% %s quantiles for Europe' % (100.*quant,qtype))
    fig.suptitle(title_,fontsize=14,y=0.96)

    picname = ('%u_percent_cumulative_%s_quantiles_vs_year' %(100*quant,qtype))+ ('_step_%u' %step)+ '.png'
    save_figure(picname)
        

def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
	
    figure(fignumber)
    savefig(path + figname, dpi=dpi)
    print 'Saved figure:',path + figname
    sys.stdout.flush()


