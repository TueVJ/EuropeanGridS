from zdcpf import *
from scipy import optimize as optimize

def sort_to_node_order(arr):
    '''For some reason, Rolando has a strange node order that does not
    match the alphabetical I am using. So, I have to convert.'''
    sortarr = []
    sortarr.append(arr[0])
    sortarr.append(arr[10])
    sortarr.append(arr[19])
    sortarr.append(arr[3])
    sortarr.append(arr[9])
    sortarr.append(arr[18])
    sortarr.append(arr[1])
    sortarr.append(arr[11])
    sortarr.append(arr[21])
    sortarr.append(arr[2])
    sortarr.append(arr[12])
    sortarr.append(arr[20])
    sortarr.append(arr[5])
    sortarr.append(arr[16])
    sortarr.append(arr[22])
    sortarr.append(arr[4])
    sortarr.append(arr[13])
    sortarr.append(arr[26])
    sortarr.append(arr[6])
    sortarr.append(arr[15])
    sortarr.append(arr[23])
    sortarr.append(arr[7])
    sortarr.append(arr[14])
    sortarr.append(arr[25])
    sortarr.append(arr[8])
    sortarr.append(arr[17])
    sortarr.append(arr[24])

    return sortarr

def get_basepath_gamma(year,filename='./results/basepath_gamma.npy'):

    print "Loading: {0}. Warning columns not pre-labeled!!".format(filename)
    data = np.load(filename)

    Gamma = zeros(len(data)-1)
    for i in arange(len(Gamma)):
        Gamma[i] = interp(year,data[0],data[i+1])
    gamma = sort_to_node_order(Gamma)
    
    return gamma

def get_basepath_alpha(year,filename='./results/basepath_alpha_w.npy'):

    print "Loading: {0}. Warning columns not pre-labeled!!".format(filename)
    data = np.load(filename)

    Alpha = zeros(len(data)-1)
    for i in arange(len(Alpha)):
        Alpha[i] = interp(year,data[0],data[i+1])
    alpha = sort_to_node_order(Alpha)

    return alpha

def generate_basepath_gamma_alpha(txtfile='../DataAndPredictionsGammaAlpha/gamma.csv',year0=1980,year_hist=2010,plot_on=True,step=2):

    print "Loading: {0}. Warning columns not pre-labeled!!".format(txtfile)
    txttitles = ['Year', 'Austria (wind)', 'Austria (solar)', 'Belgium (wind)', 'Belgium (solar)','Bulgaria (wind)','Bulgaria (solar)','Bosnia and Herzegovina (wind)','Bosnia and Herzegovina (solar)', 'Czech Republic (wind)', 'Czech Republic (solar)', 'Switzerland (wind)', 'Switzerland (solar)', 'Germany (wind)', 'Germany (solar)', 'Denmark (wind)', 'Denmark (solar)', 'Spain (wind)', 'Spain (solar)', 'France (wind)', 'France (solar)', 'Finland (wind)', 'Finland (solar)', 'Great Britain (wind)', 'Great Britain (solar)', 'Greece (wind)', 'Greece (solar)', 'Hungary (wind)', 'Hungary (solar)', 'Italy (wind)', 'Italy (solar)', 'Ireland (wind)', 'Ireland (solar)', 'Croatia (wind)', 'Croatia (solar)', 'Luxembourg (wind)', 'Luxembourg (solar)', 'Norway (wind)', 'Norway (solar)', 'Netherlands (wind)', 'Netherlands (solar)', 'Portugal (wind)', 'Portugal (solar)', 'Poland (wind)', 'Poland (solar)', 'Romania (wind)', 'Romania (solar)', 'Sweden (wind)', 'Sweden (solar)', 'Slovakia (wind)', 'Slovakia (solar)', 'Slovenia (wind)', 'Slovenia (solar)', 'Serbia (wind)', 'Serbia (solar)']
    txtlabels = ['year', 'AT-wind', 'AT-solar', 'BE-wind', 'BE-solar',
    'BG-wind', 'BG-solar', 'BA-wind', 'BA-solar', 'CZ-wind',
    'CZ-solar', 'CH-wind', 'CH-solar', 'DE-wind', 'DE-solar',
    'DK-wind', 'DK-solar', 'ES-wind', 'ES-solar', 'FR-wind',
    'FR-solar', 'FI-wind', 'FI-solar', 'GB-wind', 'GB-solar',
    'GR-wind', 'GR-solar', 'HU-wind', 'HU-solar', 'IT-wind',
    'IT-solar', 'IE-wind', 'IE-solar', 'HR-wind', 'HR-solar',
    'LU-wind', 'LU-solar', 'NO-wind', 'NO-solar', 'NL-wind',
    'NL-solar', 'PT-wind', 'PT-solar', 'PL-wind', 'PL-solar',
    'RO-wind', 'RO-solar', 'SE-wind', 'SE-solar', 'SK-wind',
    'SK-solar', 'SI-wind', 'SI-solar', 'RS-wind', 'RS-solar']
    data = np.genfromtxt(txtfile,delimiter=',',skip_header=0)
    
    p_year = array(data[0][2:-1])
    # print p_year
    year = arange(amin(p_year),amax(p_year),1)
    
    if plot_on==True:
        p_historical = p_year<=year_hist
    else:
        p_historical = None
    
    gamma, alpha_w = [], []
    for i in arange(1,data.shape[0]-1,2):
        wind = array(data[i][2:-1])
        solar = array(data[i+1][2:-1])
        if (step == 3):
            wind[-1] = data[i][-1] # take alternative values for 2050 in step 3
            solar[-1] = data[i+1][-1]

        if ~all(isnan(wind)):
            i_data = find(~isnan(wind))
            gamma_wind = get_logistic_fit(p_year[i_data],wind[i_data],year0=year0,year=year,plot_on=plot_on,p_historical=p_historical[i_data],txtlabel=txtlabels[i],txttitle=txttitles[i])
        else:
            gamma_wind = zeros(wind.shape)
            
        if ~all(isnan(solar)):
            i_data = find(~isnan(solar))
            gamma_solar = get_logistic_fit(p_year[i_data],solar[i_data],year=year,plot_on=plot_on,p_historical=p_historical[i_data],txtlabel=txtlabels[i+1],txttitle=txttitles[i+1])
        else:
            gamma_solar = zeros(wind.shape)

        gamma.append(gamma_wind+gamma_solar)
        alpha_w.append(gamma_wind/(gamma_wind+gamma_solar))

    if plot_on==True:
        weight = array([13672.325571167074,16629.057925672601,2312.6684067162068,1630.6221299198942,18730.927464722623])
        #plot_basepath_gamma_alpha(year,array(gamma),array(alpha_w),weight)

    np.save('./results/basepath_gamma',concatenate([array(year,ndmin=2),array(gamma)]))
    print 'Saved file: basepath_gamma.npy'
    np.save('./results/basepath_alpha_w',concatenate([array(year,ndmin=2),array(alpha_w)]))
    print 'Saved file: basepath_alpha_w.npy'


def plot_basepath_gamma_alpha(year,gamma,alpha_w,weight,txtlabels=None):

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    figure(1); clf()
    
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])    
    
    pp = []
    for i in arange(len(gamma)):
        pp_ = plot(year,gamma[i],'-',lw=1.5,color=colors_countries[i])
        pp.extend(pp_)
    
    gamma_mean = sum(gamma*kron(array(weight,ndmin=2).transpose(),ones(gamma.shape[1]))/sum(weight),axis=0)
    gamma_mean_DK = sum(gamma[2:4]*kron(array(weight[2:4],ndmin=2).transpose(),ones(gamma[2:4].shape[1]))/sum(weight[2:4]),axis=0)
    
    pp_mean = plot(year,gamma_mean,'k-',lw=2)
    pp.extend(pp_mean)
    
    pp_mean_DK = plot(year,gamma_mean_DK,'k--',lw=2)
    pp.extend(pp_mean_DK)
    
    
    axis(xmin=amin(year),xmax=2053,ymin=0,ymax=1.3)
    xlabel('Reference year')
    ylabel(r'Share of electricity demand ($\gamma_X$)')
    
    pp_text = ['Norway','Sweden','West Denmark','East Denmark','North Germany','Region (mean)','Denmark (mean)']
    
    leg = legend(pp,pp_text,loc='upper left')
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    #tight_layout(pad=.2)
    save_figure('plot_basepath_gamma_vs_year.png')
    
    figure(1); clf()
    
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])    
    
    pp = []
    for i in arange(len(gamma)):
        pp_ = plot(gamma_mean_DK,gamma[i],'-',lw=1.5,color=colors_countries[i])
        pp.extend(pp_)
    
    pp_mean = plot(gamma_mean_DK,gamma_mean,'k-',lw=2)
    pp.extend(pp_mean)
    
    pp_mean_DK = plot(gamma_mean_DK,gamma_mean_DK,'k--',lw=2)
    pp.extend(pp_mean_DK)
    
    
    axis(xmin=0,xmax=1.025,ymin=0,ymax=1.3)
    xlabel(r'Danish share of electricity demand ($\gamma_{\rm DK}$)')
    ylabel(r'Share of electricity demand ($\gamma_{\rm X}$)')
    
    pp_text = ['Norway','Sweden','West Denmark','East Denmark','North Germany','Region (mean)','Denmark (mean)']
    
    leg = legend(pp,pp_text,loc='upper left')
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    
    #tight_layout(pad=.2)
    save_figure('plot_basepath_gamma_vs_gamma_DK.png')

def get_logistic_fit(p_year,p_gamma,year0=1980,year=None,plot_on=False,p_historical=None,txtlabel=None,txttitle=None):
    
    if year==None:
        year = arange(amin(p_year),amax(p_year),1)
    
    fitfunc = lambda p, x: p[0]*abs(p[3])*exp(p[1]*(x-p[2]))/(abs(p[3])+p[0]*(exp(p[1]*(x-p[2]))-1.))
    # parameters:
    # fitfunc(p[2])=p[0]
    # p[1]: "slope"
    # |p[3]|=\lim_{x\to\infty} fitfunc(x)
    errfunc = lambda p, x, y, weight: (fitfunc(p, x) - y)/weight # Distance to the target function

    # p = [P_0, r, year0, K]
    # p_0 = [amin(p_gamma), .01, year0, amax(p_gamma)] # Initial guess for the parameters
    p_0 = [amin(p_gamma), .1, year0, p_gamma[-1]] # Initial guess for the parameters
    p_weight = ones(p_year.shape)
    p_fit, success = optimize.leastsq(errfunc, p_0[:], args=(p_year, p_gamma, p_weight))

    if plot_on==True:
        plot_logistic_fit(year,fitfunc(p_fit,year),p_year,p_gamma,p_historical,txtlabel,txttitle)

    return fitfunc(p_fit,year)


def plot_logistic_fit(year,gamma_fit,p_year,p_gamma,p_historical=None,txtlabel=None,txttitle=None):

    if p_historical==None:
        p_historical = ones(year.shape)

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    figure(1); clf()
    
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    
    pp_hist = plot(array(p_year[find(p_historical)]),array(p_gamma[find(p_historical)]),'go')
    pp_target = plot(p_year[find(~p_historical)],p_gamma[find(~p_historical)],'ro')
    pp_fit = plot(year,gamma_fit,'k--',lw=1.5)
    
    axis(xmin=amin(year),xmax=2053,ymin=0,ymax=1.3)
    
    xlabel('Reference year')
    ylabel(r'Share of total electricity demand ($\gamma_{\rm '+txtlabel+'}$)')


    pp = concatenate([pp_hist,pp_target,pp_fit])
    pp_text = ['Historical values','Target values','Logistic fit']
    leg = legend(pp,pp_text,loc='upper left',title=txttitle)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    #tight_layout(pad=.2)
    save_figure('plot_logistic_fit_' + txtlabel + '.png')

def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
	
    figure(fignumber)
    savefig(path + figname, dpi=dpi)
    print 'Saved figure:',path + figname
    sys.stdout.flush()
