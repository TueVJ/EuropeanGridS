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


def gamma_homogenous(linecap='copper',start=None,stop=None):
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
    N.set_alphas(0.7)
    gvals= arange(skip,skip_end,1)
    for gval in gvals:
        gamma = gval*0.01
        print "Now calculating for gamma = ",gamma
        N.set_gammas(gamma)
        F = sdcpf(N,copper=copper,h0=h0)
        name = ('homogenous_gamma_%.2f' % gamma)+'_linecap_'+linecap
        print name
        N.save_nodes_small(name+'_nodes')
        np.save('./results/'+name+'_flows',F)


def gamma_homogenous_balred(red=[0.50,0.90],guess=[0.885,0.98],start=None,stop=None):
    if start != None:
        skip = start
    else:
        skip = 0
    if stop != None:
        skip_end = stop
    else:
        skip_end=101

    gvals= arange(skip,skip_end,1)
    last_quant=guess
    quantiles=[]
    for gval in gvals:
        gamma = gval*0.01
        print "Now calculating for gamma = ",gamma
        save_filename = []
        for i in range(len(red)):
            save_filename.append(('homogenous_gamma_%.2f_balred_%.2f' % (gamma,red[i])))
        f_copper='homogenous_gamma_%.2f_linecap_copper_nodes.npz' % gamma
        f_notrans='homogenous_gamma_%.2f_linecap_0.40Q_nodes.npz' % gamma
        if (gamma >= 0.32 and gamma <= 0.38):
            last_quant=guess
        last_quant=find_balancing_reduction_quantiles(reduction=red,eps=1.e-3,guess=last_quant,stepsize=0.0025,file_copper=f_copper,file_notrans=f_notrans,gamma=gamma,alpha=None,save_filename=save_filename)
        quantiles.append(np.array(last_quant).copy())
    qsave_file='homogenous_gamma'
    if (start != None):
        qsave_file += '_from_%.2f' % min(gvals)
    if (stop != None):
        qsave_file += '_to_%.2f' % max(gvals)
    qsave_file += '_balred_quantiles'
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


def gamma_logfit_balred(red=[0.50,0.90],guess=[0.80,0.92],step=2,start=None,stop=None):
        
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
    last_quant=guess
    for year in years:
        gammas=array(get_basepath_gamma(year,step=step))
        alphas=array(get_basepath_alpha(year,step=step))
        print "Now calculating for year = ",year
        save_filename = []
        for i in range(len(red)):
            save_filename.append(('logfit_gamma_year_%u_balred_%.2f_step_%u' % (year,red[i],step)))
        f_copper='logfit_gamma_year_%u_linecap_copper_step_%u_nodes.npz' % (year,step)
        f_notrans='logfit_gamma_year_%u_linecap_0.40Q_step_%u_nodes.npz' % (year,step)
        if (year == 2018): # first year where bal. becomes significant
            last_quant=guess
        last_quant=find_balancing_reduction_quantiles(reduction=red,eps=1.e-3,guess=last_quant,stepsize=0.005,file_copper=f_copper,file_notrans=f_notrans,gamma=gammas,alpha=alphas,save_filename=save_filename)
        quantiles.append(np.array(last_quant).copy())
        print quantiles
    qsave_file='logistic_gamma'
    if (start != None):
        qsave_file += '_from_%u' % min(years)
    if (stop != None):
        qsave_file += '_to_%u' % max(years)
    qsave_file += '_balred_quantiles'
    qsave_file += '_step_%u' % step
    np.save('./results/'+qsave_file,quantiles)


######################################################
########### Process the results ######################
######################################################

def get_balancing_vs_gamma(path='./results/',prefix='homogenous_gamma',linecap='copper'):
    filename='Bvsg_'+prefix+('_linecap_%s.npy' % linecap)
    if os.path.exists(path+filename):
        Fvsg=np.load(filename)
        return Fvsg
    gammas=arange(0,1.01,0.01)
    Bvsg=np.zeros((len(gammas),2))
    j=0
    for gamma in gammas:
        Bvsg[j,0]=gamma
        load_filename=(prefix+'_%.2f_linecap_%s_nodes.npz' % (gamma,linecap))
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
    filename=path+'Fvsg_'+prefix+('_linecap_%s' % linecap)
    if os.path.exists(filename+'.npy'):
        Fvsg=np.load(filename)
        return Fvsg
    gammas=arange(0,1.01,0.01)
    Fvsg=np.zeros((len(gammas),3))
    j=0
    for gamma in gammas:
        Fvsg[j,0]=gamma
        load_filename=(prefix+'_%.2f_linecap_%s_flows.npy' % (gamma,linecap))
        print load_filename
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
    np.save(filename,Fvsg)
    return Fvsg


def get_balancing_vs_year(prefix='logfit_gamma',linecap='copper',step=2):
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
        load_filename += ('_step_%u' % step)+'_nodes.npz'
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
    filename = prefix
    if (linecap.find('balred')>=0):
        filename += '_'+linecap
    else:
        filename += '_linecap_'+linecap
    filename +=('_step_%u' % step)
    print filename
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


def get_flows_vs_year(prefix='logfit_gamma',path='./results/',linecap='copper',step=2):
    filename=path+'Fvsg_'+prefix
    if (linecap.find('balred')>=0):
        filename += '_'+linecap
    else:
        filename += '_linecap_'+linecap
    filename += '_step_%u' % step
    if os.path.exists(filename+'.npy'):
        Fvsg=np.load(filename)
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
        load_filename += ('_step_%u' % step)+'_flows.npy'
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


######################################################
########### Plot the results #########################
######################################################

def plot_balancing_vs_gamma(path='./results/',prefix='homogenous_gamma',linecaps=['0.40Q','today','balred_0.50','balred_0.90','copper'],title_=r"homogenous increase in $\gamma\,$; $ \alpha_{\rm W}=0.7 $",label=['no transmission','line capacities as of today',r'90$\,$% quantile line capacities','copper plate'],picname='balancing_vs_homogenous_gamma.png'):
    figure(1); clf()
    cl = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers
    pp = []
    i=0
    for cap in linecaps:
        fname = 'Bvsg_'+prefix
        if cap.find('balred')>=0:
            fname += '_'+cap+'.npy'
        else:
            fname += '_linecap_'+cap+'.npy'
        print 'Processing '+fname
        if not os.path.exists(path+fname):
            data = get_balancing_vs_gamma(linecap=cap)
        else:
            data = np.load(path+fname)
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


def plot_flows_vs_gamma(path='./results/',prefix='homogenous_gamma',linecaps=['0.40Q','today','0.90Q','copper'],title_=r"homogenous increase in $\gamma\,$; $ \alpha_{\rm W}=0.7 $",label=['no transmission','line capacities as of today',r'90$\,$% quantile line capacities','copper plate'],picname='flows_vs_homogenous_gamma.png'):
    figure(1); clf()
    cl = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers
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

    title(title_)
    axis(xmin=0,xmax=1,ymin=-0.05,ymax=5)
    xlabel('share $\gamma$ of VRES in total electricity production')
    ylabel(r'sum of absolute flows/$10^9\,$MW')

    pp_label = label
    leg = legend(pp,pp_label,loc='upper left');
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    save_figure(picname)
    

def plot_balancing_vs_year(path='./results/',prefix='logfit_gamma',linecaps = ['0.40Q','today','balred_0.50','balred_0.90','copper'],title_=r"Logistic growth of wind and solar generation",title_leg=r"2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $",label=['no transmission','line capacities as of today',r'50$\,$% bal. reduction quantile line capacities',r'90$\,$% bal. reduction quantile line capacities','copper plate'],picname='balancing_vs_year',step=2):
    fig=figure(1); clf()
    #calculate average gamma as weighted mean of all countries for all years
    gamma_vs_year= array(get_gamma_vs_year(step=step))
    fig.subplots_adjust(top=0.8) # leave more top margin for the extra ticks
    host = host_subplot(111)

    par = host.twin()
    axis(xmin=1990,xmax=2050,ymin=0.,ymax=0.3)

    host.set_xlabel('year')
    host.set_ylabel(r'excess balancing/av.h.l. (bal.$ - (1-\gamma)$)')
    par.set_xlabel(r'effective $\gamma$')
    gamma_ticks = arange(1990,2050+1,5)
    gamma_labels = ['%.3f' % gamma_vs_year[i-1990] for i in gamma_ticks]
    for i in range(len(gamma_ticks)):
        print gamma_ticks[i], gamma_labels[i]
    par.set_xticks(gamma_ticks)
    par.set_xticklabels(gamma_labels,rotation=30)#,shift=3)

    cl = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers
    pp = []
    i=0
    for linecap in linecaps:
        if (linecap.find('balred')>=0):
            print 'balred'
            name=prefix+'_'+linecap
        else:
            name = prefix+'_linecap_'+linecap
        name += ('_step_%u' % step)
        fname = 'Bvsg_'+name+'.npy'
        print 'Processing '+fname
        if not os.path.exists(path+fname):
            data = get_balancing_vs_year(prefix=prefix,step=step,linecap=linecap)
        else:
            data = np.load(path+fname)
        pp_x = array(data[:,0])
        pp_y = array(data[:,1]) # balancing
        # excess balancing is not so easy here
        pp_y = pp_y - (1.-gamma_vs_year)
        pp_ = host.plot(pp_x,pp_y,lw=1.5,color=cl[i])
        pp.extend(pp_)
        i += 1

    fig.suptitle(title_,fontsize=14,y=0.96)

    pp_label = label
    leg = legend(pp,pp_label,loc='upper left',title=title_leg);
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    picname = picname +'_'+ prefix + ('_step_%u' %step) + '.png'
    save_figure(picname)


def plot_flows_vs_year(path='./results/',prefix='logfit_gamma',linecaps = ['0.40Q','today','balred_0.50','balred_0.90','copper'],title_=r"Logistic growth of wind and solar generation",title_leg=r"2050 target: 100% VRES, $\alpha_{\rm W}=0.7 $",label=['no transmission','line capacities as of today',r'50$\,$% bal. reduction quantile line capacities',r'90$\,$% bal. reduction quantile line capacities','copper plate'],picname='flows_vs_year',step=2):
    fig=figure(1); clf()
    #calculate average gamma as weighted mean of all countries for all years
    gamma_vs_year= array(get_gamma_vs_year(step=step))
    fig.subplots_adjust(top=0.8) # leave more top margin for the extra ticks
    host = host_subplot(111)

    par = host.twin()
    axis(xmin=1990,xmax=2050,ymin=0.,ymax=0.3)

    host.set_xlabel('year')
    host.set_ylabel(r'excess balancing/av.h.l. (bal.$ - (1-\gamma)$)')
    par.set_xlabel(r'effective $\gamma$')
    gamma_ticks = arange(1990,2050+1,5)
    gamma_labels = ['%.3f' % gamma_vs_year[i-1990] for i in gamma_ticks]
    for i in range(len(gamma_ticks)):
        print gamma_ticks[i], gamma_labels[i]
    par.set_xticks(gamma_ticks)
    par.set_xticklabels(gamma_labels,rotation=30)#,shift=3)
    cl = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers
    pp = []
    i=0
    #calculate average gamma as weighted mean of all countries for all years
    for linecap in linecaps:
        if (linecap.find('balred')>=0):
            print 'balred'
            name=prefix+'_'+linecap
        else:
            name = prefix+'_linecap_'+linecap
        name += ('_step_%u' % step)
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
        i += 1

    fig.suptitle(title_,fontsize=14,y=0.96)

    pp_label = label
    leg = legend(pp,pp_label,loc='upper left',title=title_leg);
    axis(xmin=1990,xmax=2050,ymin=-0.05,ymax=5.)
    xlabel('year')
    ylabel(r'sum of absolute flows/$10^9\,$MW')

    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    picname = picname +'_'+ prefix + ('_step_%u' %step) + '.png'
    save_figure(picname)


def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
	
    figure(fignumber)
    savefig(path + figname, dpi=dpi)
    print 'Saved figure:',path + figname
    sys.stdout.flush()


