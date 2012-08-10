from jobscript import *

linecap=['copper','today','0.40Q']
step=[2,21,22]#,3,31,32]
step=[3,31,32]

for l in linecap:
    for s in step:
        gamma_logfit(linecap=l,step=s)

# these files are needed as input to the balred calculations!
for i in linecap:
    # get_balancing_vs_gamma(linecap=i)
    for j in step:
        get_balancing_vs_year(linecap=i,step=j)

for s in step:
    gamma_logfit_balred(step=s)

linecap=['balred_0.70','balred_0.90']
for i in linecap:
    for s in step:
        get_balancing_vs_year(linecap=i,step=s,capped_inv=False)

for s in step:
    gamma_logfit_balred_capped(step=s)

linecap=['balred_0.70','balred_0.90']
for i in linecap:
    for s in step:
        get_balancing_vs_year(linecap=i,step=s,capped_inv=True)

linecap=['copper','today','0.40Q']
for l in linecap:
    for s in step:
        get_flows_vs_year(linecap=l,step=s)

linecap=['balred_0.70','balred_0.90']
for l in linecap:
    for s in step:
        get_flows_vs_year(linecap=l,step=s,capped_inv=False)
        get_flows_vs_year(linecap=l,step=s,capped_inv=True)

linecap=['copper','today','0.40Q']
for l in linecap:
    for s in step:
        get_balancing_quantiles_vs_year(step=s,linecap=l)

linecap=['balred_0.70','balred_0.90']
for l in linecap:
    for s in step:
        get_balancing_quantiles_vs_year(step=s,linecap=l,capped_inv=False)
        get_balancing_quantiles_vs_year(step=s,linecap=l,capped_inv=True)

linecap=['balred_0.70','balred_0.90']
for s in step:
    for l in linecap:
        get_import_and_deficit(step=s,linecap=l)
        get_curtailment_and_excess(step=s,linecap=l)

# plot_balancing_vs_year(step=2,capped_inv=False)
# plot_balancing_vs_year(step=2,capped_inv=True)
# plot_flows_vs_year(step=2,capped_inv=False)
# plot_flows_vs_year(step=2,capped_inv=True)
# plot_linecaps_vs_year(step=2,capped_inv=False)
# plot_linecaps_vs_year(step=2,capped_inv=True)
# plot_investment_vs_year(step=2,capped_inv=False)
# plot_investment_vs_year(step=2,capped_inv=True)
    
