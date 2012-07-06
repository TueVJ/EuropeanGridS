from zdcpf import *
from scipy import optimize as optimize
from Database_v2 import *


######################################################
########### Access the results #######################
######################################################

def get_basepath_gamma(year,filename='./results/basepath_gamma',step=2,combifit=False):
    '''Return basepath gamma for year year and the basepath specified by
    filename, step, and combifit.'''
    filename += '_step_%u' % step
    if (combifit):
        filename += '_combined_fit.npy'
    else:
        filename += '_separate_fit.npy'

    print "Loading: {0}. Warning columns not pre-labeled!!".format(filename)
    data = np.load(filename)

    Gamma = zeros(len(data)-1)
    for i in arange(len(Gamma)):
        Gamma[i] = interp(year,data[0],data[i+1])
    gamma = sort_to_node_order(Gamma)
    
    return gamma


def get_basepath_alpha(year,filename='./results/basepath_alpha_w',step=2,combifit=False):
    '''Return basepath alpha for year year and the basepath specified by
    filename, step, and combifit.'''
    filename += '_step_%u' % step
    if (combifit):
        filename += '_combined_fit.npy'
    else:
        filename += '_separate_fit.npy'

    print "Loading: {0}. Warning columns not pre-labeled!!".format(filename)
    data = np.load(filename)

    Alpha = zeros(len(data)-1)
    for i in arange(len(Alpha)):
        Alpha[i] = interp(year,data[0],data[i+1])
    alpha = sort_to_node_order(Alpha)

    return alpha


######################################################
########### Generate the results #####################
######################################################

def generate_basepath_gamma_alpha(txtfile='../DataAndPredictionsGammaAlpha/gamma.csv',year0=1980,year_hist=2010,plot_on=True,combifit=False,step=2):

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
    
    p_year = array(data[0][2:25])
    # print p_year
    # print data[0][23:]
    # print data[1][23:]
    year = arange(amin(p_year),amax(p_year)+1,1)
    
    if plot_on==True:
        p_historical = p_year<=year_hist
    else:
        p_historical = None
    #print 'p_historical ',p_historical
    #print len(arange(1,data.shape[0]-1,2))
    
    gamma, alpha_w = [], []
    for i in arange(1,data.shape[0]-1,2): # loop over countries
        wind = array(data[i][2:25])
        solar = array(data[i+1][2:25])
        for j in arange(21,26,1):
            if (step == j):
                wind[-1] = data[i][j+4] # take alternative values for 2050
                solar[-1] = data[i+1][j+4]
        if (step == 3):
            wind[-1] = data[i][30] # take alternative values for 2050
            solar[-1] = data[i+1][30]
        for j in arange(31,36,1):
            if (step == j):
                wind[-1] = data[i][j] # take alternative values for 2050
                solar[-1] = data[i+1][j]
        if (step == 4):
            wind[-1] = data[i][-1] # take alternative values for 2050 in step 4
            solar[-1] = data[i+1][-1]
        # print 'country: ',txttitles[i]
        # print 'wind: ',wind[-1]
        # print 'solar: ',solar[-1]

        if (combifit and txtlabels[i] != 'ES-wind'): # ES is messed up in combifit, so resort to separate fit in this case
            i_data_w = find(~isnan(wind))
            i_data_s = find(~isnan(solar))
            i_data_ws = []
            for j in range(len(i_data_w)):
                if i_data_w[j] in i_data_s:
                    i_data_ws.append(i_data_w[j])
            i_data_ws=array(i_data_ws)

            wind = wind[i_data_ws]
            solar = solar[i_data_ws]
            p_gamma=(wind+solar)
            for j in range(len(p_gamma)): # avoid division by zero
                if p_gamma[j]<1.e-6:
                    p_gamma[j]=1.e-6
            p_alpha_w=zeros(len(p_gamma))
            for j in range(len(p_gamma)):
                p_alpha_w[j]=wind[j]/p_gamma[j]
            year, f_gamma, f_alpha_w = get_wind_solar_logistic_fit(p_year[i_data_ws],p_gamma,p_alpha_w,year0=year0,year=year)            

            if plot_on==True:
                plot_logistic_fit(year,f_gamma*f_alpha_w,p_year[i_data_ws],p_gamma*p_alpha_w,p_historical[i_data_ws],txtlabel=txtlabels[i],txttitle=txttitles[i],step=step,combifit=combifit)
                plot_logistic_fit(year,f_gamma*(1-f_alpha_w),p_year[i_data_ws],p_gamma*(1-p_alpha_w),p_historical[i_data_ws],txtlabel=txtlabels[i+1],txttitle=txttitles[i+1],step=step,combifit=combifit)

            gamma.append(f_gamma)
            alpha_w.append(f_alpha_w)

        else:
            if ~all(isnan(wind)):
                i_data = find(~isnan(wind))
                year, gamma_wind = get_logistic_fit(p_year[i_data],wind[i_data],year0=year0,year=year)
                if plot_on==True:
                    plot_logistic_fit(year,gamma_wind,p_year[i_data],wind[i_data],p_historical[i_data],txtlabel=txtlabels[i],txttitle=txttitles[i],step=step,combifit=combifit)
            else:
                gamma_wind = zeros(wind.shape)
        
            if ~all(isnan(solar)):
                i_data = find(~isnan(solar))
                year, gamma_solar = get_logistic_fit(p_year[i_data],solar[i_data],year0=year0,year=year)
                if plot_on:
                    plot_logistic_fit(year,gamma_solar,p_year[i_data],solar[i_data],p_historical[i_data],txtlabel=txtlabels[i+1],txttitle=txttitles[i+1],step=step,combifit=combifit)
                    
            else:
                gamma_solar = zeros(wind.shape)

            gamma.append(gamma_wind+gamma_solar)
            alpha_w.append(gamma_wind/(gamma_wind+gamma_solar))

    save_filename='basepath_gamma_step_%u' % step
    if combifit:
        save_filename += '_combined_fit'
    else:
        save_filename += '_separate_fit'
    np.save('./results/'+save_filename,concatenate([array(year,ndmin=2),array(gamma)]))
    print 'Saved file: '+save_filename
    save_filename = 'basepath_alpha_w_step_%u' %step
    if combifit:
        save_filename += '_combined_fit'
    else:
        save_filename += '_separate_fit'
    np.save('./results/'+save_filename,concatenate([array(year,ndmin=2),array(alpha_w)]))
    print 'Saved file: '+save_filename


def generate_basepath_gamma_optimal_alpha(txtfile='../DataAndPredictionsGammaAlpha/gamma.csv',year0=1980,year_hist=2010,CS=None,rel_tol=1.e-2,plot_on=True,combifit=False,step=22):

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
    
    p_year = array(data[0][2:-4])
    # print p_year
    year = arange(amin(p_year),amax(p_year)+1,1)
    
    if plot_on==True:
        p_historical = p_year<=year_hist
    else:
        p_historical = None
    #print 'p_historical ',p_historical
    #print len(arange(1,data.shape[0]-1,2))
    
    gamma, alpha_w = [], []
    for i in arange(1,data.shape[0]-1,2): # loop over countries arange(1,7,2):
        wind = array(data[i][2:-4])
        solar = array(data[i+1][2:-4])
        if (step == 3):
            wind[-1] = data[i][-4] # take alternative values for 2050 in step 3
            solar[-1] = data[i+1][-4]
        elif (step == 4):
            wind[-1] = data[i][-3] # take alternative values for 2050 in step 4
            solar[-1] = data[i+1][-3]
        elif (step == 22):
            wind[-1] = data[i][-2] # take alternative values for 2050 in step 22
            solar[-1] = data[i+1][-2]
        elif (step == 32):
            wind[-1] = data[i][-1] # take alternative values for 2050 in step 32
            solar[-1] = data[i+1][-1]
        # print 'country: ',txttitles[i]
        # print 'wind: ',wind
        # print 'solar: ',solar

        p_gamma = wind+solar
        p_alpha = zeros_like(p_gamma)
        for j in range(len(p_gamma)): # avoid division by zero
            if p_gamma[j]<1.e-6:
                p_gamma[j]=1.e-6
            p_alpha[j] = wind[j]/p_gamma[j]
        i_data = find(~isnan(p_gamma))
        ISO = txtlabels[i][0:2]
        yr, f_gamma, f_alpha_w = get_optimal_path_logistic_fit(p_year[i_data], p_gamma[i_data], ISO=ISO, year0=1980., year=year, CS=CS, rel_tol=rel_tol)
        if plot_on==True:
            #plot_logistic_fit(year,f_gamma*f_alpha_w,p_year[i_data],wind[i_data],p_historical[i_data],txtlabel=txtlabels[i],txttitle=txttitles[i],step=step,combifit='optimal_alpha')
            #plot_logistic_fit(year,f_gamma*(1-f_alpha_w),p_year[i_data],solar[i_data],p_historical[i_data],txtlabel=txtlabels[i+1],txttitle=txttitles[i+1],step=step,combifit='optimal_alpha')
            plot_optimal_alpha_logistic_fit(year,f_gamma,f_alpha_w,p_year[i_data],p_gamma[i_data],p_alpha[i_data],p_historical[i_data],txtlabel=ISO,txttitle=txttitles[i][:-7],step=step)


        gamma.append(f_gamma)
        alpha_w.append(f_alpha_w)

    save_filename='basepath_gamma_step_%u' % step
    save_filename += '_optimal_alpha'
    np.save('./results/'+save_filename,concatenate([array(year,ndmin=2),array(gamma)]))
    print 'Saved file: '+save_filename
    save_filename = 'basepath_alpha_w_step_%u' %step
    save_filename += '_optimal_alpha'
    np.save('./results/'+save_filename,concatenate([array(year,ndmin=2),array(alpha_w)]))
    print 'Saved file: '+save_filename

    return year,gamma,alpha_w

######################################################
########### Various fit functions ####################
######################################################

def get_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w, year0=1980, year=None):
    """Combine targets for both gamma and alpha_w in one fit. Both wind and solar are required to follow a logistic growth."""

    p_year, p_gamma, p_alpha_w = array(p_year,ndmin=1), array(p_gamma,ndmin=1), array(p_alpha_w,ndmin=1)

    if year==None:
        year = arange(amin(p_year),amax(p_year)+1,1)

    f_logistic = lambda p, x, x0: abs(p[0])*abs(p[2])*exp(abs(p[1])*(x-x0))/(abs(p[2])+abs(p[0])*(exp(abs(p[1])*(x-x0))-1.)) # p is a length 3 array

    f_gamma_w = lambda pp, year: f_logistic(pp[0:3],year,year0)
    f_gamma_s = lambda pp, year: f_logistic(pp[3:6],year,year0)
    f_gamma = lambda pp, year: (f_gamma_w(pp,year) + f_gamma_s(pp,year)) if (f_gamma_w(pp,year)[0] + f_gamma_s(pp,year)[0] >= 1.e-6) else 1.e-6
    f_alpha_w = lambda pp, year: f_gamma_w(pp,year)/f_gamma(pp,year)
    
    #errfunc = lambda pp, year, gamma, alpha_w: concatenate([(f_gamma(pp,year)-gamma),(f_alpha_w(pp,year)-alpha_w)]) # pp is a length 6 array.
    
    #This error function treats wind and solar in the same way. In the one above, emphasis can be placed on either mix or on gamma.
    errfunc = lambda pp, year, gamma, alpha_w: concatenate([(f_alpha_w(pp,year)*f_gamma(pp,year)-alpha_w*gamma),((1-f_alpha_w(pp,year))*f_gamma(pp,year)-(1-alpha_w)*gamma)]) # pp is a length 6 array.

    pp_0 = [amin(p_alpha_w*p_gamma), .01, amax(p_alpha_w*p_gamma), amin((1-p_alpha_w)*p_gamma), .01, amax((1-p_alpha_w)*p_gamma)] # Initial guess for the parameters
    pp_fit, success = optimize.leastsq(errfunc, pp_0[:], args=(p_year, p_gamma, p_alpha_w))

    return year, f_gamma(pp_fit,year), f_alpha_w(pp_fit,year)


def get_logistic_fit(p_year,p_gamma,year0=1980,year=None):
    
    p_year, p_gamma = array(p_year,ndmin=1), array(p_gamma,ndmin=1)

    if year==None:
        year = arange(amin(p_year),amax(p_year)+1,1)
    
    fitfunc = lambda p, x: p[0]*abs(p[3])*exp(p[1]*(x-p[2]))/(abs(p[3])+p[0]*(exp(p[1]*(x-p[2]))-1.))
    # parameters:
    # fitfunc(p[2])=p[0]
    # p[1]: "slope"
    # |p[3]|=\lim_{x\to\infty} fitfunc(x)
    errfunc = lambda p, x, y: (fitfunc(p, x) - y) # Distance to the target function

    p_0 = [amin(p_gamma), .1, year0, p_gamma[-1]] # Initial guess for the parameters
    p_fit, success = optimize.leastsq(errfunc, p_0[:], args=(p_year, p_gamma))

    return year, fitfunc(p_fit,year)


#
# get_optimal_path_logistic_fit(p_year=[1980,2000,2010,2020,2050], p_gamma=[0.01,0.1,.25,0.5,1.0]);
#
def get_optimal_path_logistic_fit(p_year, p_gamma, ISO='DK', year0=1980., year=None, CS=None, rel_tol=1e-2):
    """Combines targets for gamma with the optimal path to get an optimal build-up of wind and solar. Both follow a logistic growth."""

    p_year, p_gamma = array(p_year,ndmin=1), array(p_gamma,ndmin=1)
    
    if year==None:
        year = arange(amin(p_year),amax(p_year)+1,1)
        
    t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO)
    
    ## Initial estimate of gamma vs year
    p_alpha_w_opt = get_optimal_path_balancing(L,Gw,Gs,p_gamma,CS=CS)
    #print p_year,p_gamma,p_alpha_w_opt
    year, gamma, alpha_w = get_wind_solar_logistic_fit(p_year, p_gamma, p_alpha_w_opt, year0, year)
    #print year, gamma, alpha_w
    for i in range(len(gamma)):
        if gamma[i] < 1.e-6:
            gamma[i] =1.e-6
    
    ## Iterate to get better match with the optimal mix. The gamma vs year estimate is used to provide a better and continous match to the optimal path. Otherwise few targets can result in large deviations from the optimal path.
    rel_err = 1000; rel_goodness_old = 1000; i=0
    while (rel_err>rel_tol or i<3):
        i+=1
        gamma_w_old = gamma*alpha_w
        
        alpha_w_opt = get_optimal_path_balancing(L,Gw,Gs,gamma,CS=CS)
        year, gamma, alpha_w = get_wind_solar_logistic_fit(year, gamma, alpha_w_opt, year0, year)
        # for i in range(len(gamma)):
        #     if gamma[i] < 1.e-6:
        #         gamma[i] = 1.e-6

        rel_goodness = amax(abs(gamma_w_old - gamma*alpha_w))
        rel_err = abs(rel_goodness - rel_goodness_old)
        rel_goodness_old = rel_goodness    
    
    return year, gamma, alpha_w


def get_optimal_path_balancing(L, GW, GS, gamma=linspace(0,1,5), p_interval=0.01, CS=None, returnall=False, normalized=True):
    """Wraper for get_optimal_mix_balancing(). This function allows gamma to be an array."""

    gamma = array(gamma,ndmin=1)
    
    if returnall==True:
        alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p = get_optimal_mix_balancing(L, GW, GS, gamma[0], p_interval, CS, returnall, normalized)
        alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p = expand2array(alpha_w_opt,gamma), expand2array(alpha_w_opt_1p_interval,gamma), expand2array(res_load_sum_opt,gamma), expand2array(mismatch_opt,gamma), expand2array(res_load_sum_1p,gamma)
        
        for i in arange(1,len(gamma)):
            alpha_w_opt[i], alpha_w_opt_1p_interval[i], res_load_sum_opt[i], mismatch_opt[i], res_load_sum_1p[i] = get_optimal_mix_balancing(L, GW, GS, gamma[i], p_interval, CS, returnall, normalized)
        
        return alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p

    else:
        alpha_w_opt = zeros(gamma.shape)
        for i in arange(len(gamma)):
            alpha_w_opt[i] = get_optimal_mix_balancing(L, GW, GS, gamma[i], p_interval, CS, returnall, normalized)
        
        return alpha_w_opt
        

def get_optimal_mix_balancing(L, GW, GS, gamma=1., p_interval=0.01, CS=None, returnall=False, normalized=True, DefaultWind=True):

    L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l

    if CS==None:
        res_load_sum = lambda alpha_w: sum(get_positive(-mismatch(alpha_w))) - alpha_w*0.001*sign(DefaultWind-.5)
    else:
        res_load_sum = lambda alpha_w: sum(get_positive(-get_policy_2_storage(mismatch(alpha_w),storage_capacity = CS)[0])) - alpha_w*0.001*sign(DefaultWind-.5)
    
    alpha_w_opt = optimize.fmin(res_load_sum,0.5,disp=False)

    if alpha_w_opt>1.:
        alpha_w_opt = 1.
    elif alpha_w_opt<0.:
        alpha_w_opt = 0.

    if normalized:
        if CS==None:
            mismatch_opt = mismatch(alpha_w_opt)
        else:
            mismatch_opt = get_policy_2_storage(mismatch(alpha_w_opt),storage_capacity = CS)[0]
    else:
        if CS==None:
            mismatch_opt = mismatch(alpha_w_opt)*mean(sum(L,axis=0))
        else:
            mismatch_opt = get_policy_2_storage(mismatch(alpha_w_opt),storage_capacity = CS)[0]*mean(sum(L,axis=0))
            
    res_load_sum_opt = res_load_sum(alpha_w_opt)

    if returnall:
        res_load_sum_1p_interval = lambda alpha_w: res_load_sum(alpha_w)-(res_load_sum(alpha_w_opt)+p_interval*sum(l))
        
        if sign(res_load_sum_1p_interval(0))!=sign(res_load_sum_1p_interval(alpha_w_opt)):
            lower_bound = brentq(res_load_sum_1p_interval, 0, alpha_w_opt)
        else:
            lower_bound = 0
        
        if sign(res_load_sum_1p_interval(1))!=sign(res_load_sum_1p_interval(alpha_w_opt)):
            upper_bound = brentq(res_load_sum_1p_interval, alpha_w_opt, 1)
        else:
            upper_bound = 1
        
        alpha_w_opt_1p_interval = array([lower_bound,upper_bound])
        res_load_sum_1p = amax([res_load_sum(lower_bound),res_load_sum(upper_bound)])
        
        #Returns: alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt
        return alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt, res_load_sum_1p
    else:
        return alpha_w_opt


######################################################
########### Plot functions ###########################
######################################################

def plot_logistic_fit(year,gamma_fit,p_year,p_gamma,p_historical=None,txtlabel=None,txttitle=None,step=2,combifit=False):

    if p_historical==None:
        p_historical = ones(year.shape)

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    figure(1); clf()
    
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    #print 'p_year', p_year
    #print 'p_year[find(p_historical)]', p_year[find(p_historical)]
    #print 'p_gamma[find(p_historical)]', p_gamma[find(p_historical)]
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
    figname = 'plot_logistic_fit_'
    if (combifit==True):
        figname += 'combined'
    elif (combifit==False):
        figname += 'separate'
    else:
        figname += combifit
    figname += ('_step_%u_' % step) + txtlabel + '.png'
    save_figure(figname)


def plot_optimal_alpha_logistic_fit(year,f_gamma,f_alpha,p_year,p_gamma,p_alpha,p_historical=None,txtlabel=None,txttitle=None,step=22):

    if p_historical==None:
        p_historical = ones(year.shape)

    #Set plot options	
    matplotlib.rcParams['font.size'] = 10

    figure(1); clf()
    
    gcf().set_dpi(300)
    gcf().set_size_inches([5.25,3.5])
    #print 'p_year', p_year
    #print 'p_year[find(p_historical)]', p_year[find(p_historical)]
    #print 'p_gamma[find(p_historical)]', p_gamma[find(p_historical)]
    pp_hist = plot(p_year[find(p_historical)],p_gamma[find(p_historical)],'go')
    pp_target = plot(p_year[find(~p_historical)],p_gamma[find(~p_historical)],'ro')
    pp_alpha_h = plot(p_year[find(p_historical)],p_gamma[find(p_historical)]*p_alpha[find(p_historical)],color='black',marker='.',ls='')
    pp_alpha_t = plot(p_year[find(~p_historical)],p_gamma[find(~p_historical)]*p_alpha[find(~p_historical)],color='gray',marker='.',ls='')
    pp_fit = plot(year,f_gamma,'k--',lw=1.5)
    pp_fit_a = plot(year,f_gamma*f_alpha,color='gray',ls=':',lw=1.5)
    
    axis(xmin=amin(year),xmax=2053,ymin=0,ymax=1.3)
    
    xlabel('Reference year')
    ylabel(r'Share of total electricity demand ($\gamma_{\rm '+txtlabel+'}$)')

    pp = concatenate([pp_hist,pp_target,pp_alpha_h,pp_alpha_t,pp_fit,pp_fit_a])
    pp_text = [r'Historical total share $\gamma$','Targeted total share $\gamma$',r'Historical wind share $\gamma_{\rm W}$',r'Targeted wind share $\gamma_{\rm W}$',r'Combined logistic fit to $\gamma$',r'Optimal $\gamma_{\rm W}$ path']
    leg = legend(pp,pp_text,loc='upper left',title=txttitle)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize

    #tight_layout(pad=.2)
    figname = 'plot_optimal_alpha_logistic_fit'
    figname += ('_step_%u_' % step) + txtlabel + '.png'
    save_figure(figname)


######################################################
########### Little helpers ###########################
######################################################

def sort_to_node_order(arr):
    '''Rolando has his own node order that does not
    match the one I am using here. So, I have to convert.'''
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


# def save_figure(figname='TestFigure.png', fignumber=gcf().number, path='./figures/', dpi=300):
	
#     figure(fignumber)
#     savefig(path + figname, dpi=dpi)
#     print 'Saved figure:',path + figname
#     sys.stdout.flush()

