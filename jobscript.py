from zdcpf import *
import scipy.optimize as optimize
from logistic_gammas import *

def copper_flow():
    N = Nodes()
    N.set_alphas(0.7)
    N.set_gammas(1.0)
    N,F = sdcpf(N,copper=1)
    N.save_nodes('copper_nodes')
    save('./results/'+'copper_flows',F)

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
        N,F = sdcpf(N,copper=copper,h0=h0)
        name = ('homogenous_gamma_%.2f' % gamma)+'_linecap_'+linecap
        print name
        N.save_nodes(name+'_nodes')
        save('./results/'+name+'_flows',F)

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
        N,F = sdcpf(N,h0=h0,copper=copper)
        name = 'logfit_gamma_year_%u_linecap_%s_step_%u' % (year,linecap,step)
        print name
        N.save_nodes(name+'_nodes')
        save('./results/'+name+'_flows',F)
        del N


def get_balancing_vs_gamma(prefix='homogenous_gamma',linecap='copper'):
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
    filename=prefix+('_linecap_%s' % linecap)
    save('./results/Bvsg_'+filename,Bvsg)
    return Bvsg

def plot_balancing_vs_gamma(filenames=['quant_0.40_gamma','today_linecap_gamma','quant_0.90_gamma','copper_gamma'],title_=r"homogenous increase in $\gamma\,$; $ \alpha_{\rm W}=0.7 $",label=['no transmission','line capacities as of today',r'90$\,$% quantile line capacities','copper plate'],picname='balancing_vs_gamma.png'):
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
    save(filename,Fvsg)
    return Fvsg
    
def get_balancing_vs_year(prefix='logfit_gamma',linecap='copper',step=2):
    years=arange(1990,2050+1,1)
    Bvsg=np.zeros((len(years),2))
    j=0
    for year in years:
        Bvsg[j,0]=year
        load_filename = prefix+('_year_%u' %year)
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
    filename = prefix+'_linecap_'+linecap+('_step_%u' % step)
    print filename
    save('./results/Bvsg_'+filename,Bvsg)
    return Bvsg

def get_gamma_vs_year(path='./results/',prefix='logfit_gamma',step=2):
    filename = prefix+'_vs_year_step_%u.npy' % step
    if os.path.exists(path+filename):
        gamma_vs_year = np.load(filename)
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
        save(filename,gamma_vs_year)
    return gamma_vs_year
    

def plot_balancing_vs_year(prefix='gamma_logfit',title_=r"logarithmic growth of $\gamma\,$; final $\alpha_{\rm W}=0.7 $",label=['no transmission','line capacities as of today',r'90$\,$% quantile line capacities','copper plate'],picname='balancing_vs_gamma_logfit',step=2):
    figure(1); clf()
    cl = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers
    pp = []
    i=0
    #calculate average gamma as weighted mean of all countries for all years
    linecaps = ['0.40Q','today','0.90Q','copper']
    for linecap in linecaps:
        name = prefix+'_linecap_'+linecap+('_step_%u' % step)
        fname = './results/Bvsg_'+name+'.npy'
        print 'Processing '+fname
        if not os.path.exists(fname):
            data = get_balancing_vs_year(prefix=prefix,step=step,linecap=linecap)
        else:
            data = np.load(fname)
        pp_x = array(data[:,0])
        pp_y = array(data[:,1]) # balancing
        # excess balancing is not so easy here
        gamma_vs_year= array(get_gamma_vs_year(step=step))
        pp_y = pp_y - (1.-gamma_vs_year)
        pp_ = plot(pp_x,pp_y,lw=1.5,color=cl[i])
        pp.extend(pp_)
        i += 1

    title(title_)
    axis(xmin=1990,xmax=2050,ymin=0.,ymax=0.3)
    xlabel('year')
    ylabel(r'excess balancing/av.h.l. (bal.$ - (1-\gamma)$)')

    pp_label = label
    leg = legend(pp,pp_label,loc='upper left');
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    legend()

    picname = picname + ('_step_%u' %step) + '.png'
    save_figure(picname)

def get_flows_vs_year(prefix='logfit_gamma',path='./results/',linecap='copper',step=2):
    filename=path+'Fvsg_'+prefix+('_linecap_%s_step_%u' % (linecap,step))
    if os.path.exists(filename+'.npy'):
        Fvsg=np.load(filename)
        return Fvsg
    years=arange(1990,2050+1,1)
    Fvsg=np.zeros((len(years),3))
    j=0
    for year in years:
        Fvsg[j,0]=year
        load_filename=(prefix+'_year_%u_linecap_%s_step_%u_flows.npy' % (year,linecap,step))
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
    save(filename,Fvsg)
    return Fvsg

'''def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
	
    figure(fignumber)
    savefig(path + figname, dpi=dpi)
    print 'Saved figure:',path + figname
    sys.stdout.flush()

'''
