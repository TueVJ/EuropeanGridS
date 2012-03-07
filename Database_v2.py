#
#  Database_v2.py
#  Fetch data from the pepsi sql database.
#
#  Created by Gorm Bruun Andresen on 16/12/2010.
#  Modified by Sarah Becker on 01/03/2012.
#  Copyright (c) 2010 University of Aarhus. All rights reserved.
#

import numpy
from pylab import *

from sqlalchemy.engine import reflection
from sqlalchemy import create_engine
from sqlalchemy.schema import (MetaData,Table)

#cp: engine, inspector, metadata = connect2database()
def connect2database(url='postgresql://dheide:050877-4125@localhost/iset_data',localhost=True,echo=False):

	if localhost:
		## Works from anywhere else if you setup an ssh tunnel first: ssh -L5432:localhost:5432 username@pepsi.imf.au.dk
		engine = create_engine(url, echo=echo)
	else:
		## Works when you are on the IMF intranet.
		engine = create_engine(url, echo=echo)

	inspector = reflection.Inspector.from_engine(engine)

	metadata = MetaData(engine)
	
	return engine, inspector, metadata

def get_all_table_names(inspector=None):

	if inspector==None:
		engine, inspector, metadata = connect2database(localhost=True);

	for schema_name in inspector.get_schema_names():
		print 'Tables in schema:',schema_name
		print '---------------------'
		for table_name in inspector.get_table_names(schema=schema_name):
			print table_name
		print
	
def get_datatype(table='europe_load_dt',schema='agg_avg_1hour_pdata_caps_eu2020', inspector=None):

	if inspector==None:
		engine, inspector, metadata = connect2database(localhost=True);
	
	datatype = inspector.get_columns(table,schema=schema)
	
	return datatype
	
def get_data(table='europe_load_dt',schema='agg_avg_1hour_pdata_caps_eu2020', metadata=None):
	
	if metadata==None:
		engine, inspector, metadata = connect2database(localhost=True);
	
	data = Table(table,metadata,schema=schema,autoload=True)
	s = data.select()
	rs = s.execute() 
	
	data = list(empty(rs.rowcount))
	for i in arange(rs.rowcount):
		data[i] = list(rs.fetchone())
	
	return_data = transposed(data)
	#data = rs.fetchall()
	
	return return_data

def transposed(lists):
   if not lists: return []
   return map(lambda *row: list(row), *lists)
	
def get_data_countries(schema='agg_avg_1hour_pdata_caps_eu2020', metadata=None, inspector=None, localhost=True):

	if metadata==None or inspector==None:
		engine, inspector, metadata = connect2database(localhost=True);
	
	#get labels
	datatype = get_datatype('countries_load_dt',schema, inspector)
	datalabels = list(empty(size(datatype)-1));
	for i in arange(0,size(datalabels)):
		datalabels[i] = datatype[i+1]['name']
		
	#get data
	load = get_data('countries_load_dt',schema,metadata)
	time = load[0];
	datetime_offset = date2num(time[0]);
	time = date2num(time) - datetime_offset #Use num2date(time+datetime_offset) to convert back to date-type objects.
	
	load = array(load[1:]);
	wind = array(get_data('countries_wind',schema,metadata)[1:])
	solar = array(get_data('countries_pv',schema,metadata)[1:])
	
	return time, load, wind, solar, datetime_offset, datalabels
	
def get_data_regions(schema='agg_avg_1hour_pdata_caps_eu2020', metadata=None, inspector=None, localhost=True):

	if metadata==None or inspector==None:
		engine, inspector, metadata = connect2database(localhost=True);
	
	#get labels
	datatype = get_datatype('regions_load_dt',schema, inspector)
	datalabels = list(empty(size(datatype)-1));
	for i in arange(0,size(datalabels)):
		datalabels[i] = datatype[i+1]['name']
		
	#get data
	load = get_data('regions_load_dt',schema,metadata)
	time = load[0];
	datetime_offset = date2num(time[0]);
	time = date2num(time) - datetime_offset #Use num2date(time+datetime_offset) to convert back to date-type objects.
	
	load = array(load[1:]);
	wind = array(get_data('regions_wind',schema,metadata)[1:])
	solar = array(get_data('regions_pv',schema,metadata)[1:])
	
	return time, load, wind, solar, datetime_offset, datalabels

	
def plot_load_wind_solar_countries(Ndays=365, localhost=True):
	"""Test plot load, wind, and solar data for all countries."""
	
	
	time, load, wind, solar, datetime_offset, datalabels = get_data_countries(localhost=localhost);
	
	figure(1); clf();
	for i in range(load.shape[0]):
		plot(time[:Ndays*24],load[i][:Ndays*24]/sum(load[i]),label=datalabels[i])	

	xlabel('Time [hours]')
	ylabel('Load [norm.]')
	
	legend()
	
	savefig('CountriesLoad.png')
	
	figure(1); clf();
	for i in range(load.shape[0]):
		plot(time[:Ndays*24],wind[i][:Ndays*24]/sum(wind[i]),label=datalabels[i])	

	xlabel('Time [hours]')
	ylabel('Wind power [norm.]')
	
	legend()
	
	savefig('CountriesWind.png')
	
	figure(1); clf();
	for i in range(load.shape[0]):
		plot(time[:Ndays*24],solar[i][:Ndays*24]/sum(solar[i]),label=datalabels[i])	
	
	xlabel('Time [hours]')
	ylabel('Solar power [norm.]')
	
	legend()
	
	savefig('CountriesSolar.png')



######
# Convenient access to ISET country or region data.
# Main functions:  get_ISET_country_data(), get_ISET_region_data(),
# get_ISET_data()
###

def get_ISET_data(name='AT',path='./data/'):
    """
    Wrapper around get_ISET_country_data() and get_ISET_region_data()
    Determines which one of those functions to call depending on the 
    name of the region/country
    """

    if valid_ISO(name):
        t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data(ISO=name,path=path)
    elif valid_REG(name):
        t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_region_data(REG=name,path=path)
    else:
        sys.exit("Error (45nlksd): No such country or region ({0}). For a list of countries, use get_ISET_country_names(), for regions use get_ISET_region_names().".format(name))

    return t, L, Gw, Gs, datetime_offset, datalabel

def get_ISET_country_data(ISO='DK',path='./data/'):
    """
    Returns data for a specific country in the ISET data set. 
    
    The country is specified using ISO two letter names. For a list of names use get_ISET_country_names().
    
    The function retrives/stores data from/in ./data as default. If the data is not available the function
     attempts to download it from pepsi. This requires ssh access: 
     ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk 
     
    Filenames are: ISET_country_ISO.npz
    
    Returns
    -------
    t: array
        hours numbered from 0 to N-1
    L: array
        load in MW
    Gw: array
        normalized wind power generation
    Gs: array
        normalized solar power generation
    datetime_offset: scalar
        num2date(datetime_offset) yields the date and hour of t=0
    datalabel: string
        country ISO
        
    Example: Returns the Danish load, wind and solar power generation
        t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_country_data('DK')
    
    """
    
    if not valid_ISO(ISO):
        sys.exit("Error (43nlksd): No such country ISO ({0}). For a list of names use get_ISET_country_names().".format(ISO))
    
    filename = 'ISET_country_' + ISO + '.npz'
    
    try:
        #Load the data file if it exists:
        npzfile = np.load(path + filename)
        print 'Loaded file: ', path + filename
        sys.stdout.flush()
        
    except IOError:
        print 'Datafile does not exist:', path + filename
        print 'Trying to download data from pepsi...'
        sys.stdout.flush()
        try: 
            #Normalized data: t, l, Gw, Gs, datetime_offset, datalabels = get_data_countries(schema='norm_agg_avg_1hour_pdata_caps_eu2020',localhost=True);
            t, L, GW, GS, datetime_offset, datalabels = get_data_countries(localhost=True);
        except:
            sys.exit("Error (sdf3dz1): Could not connect to pepsi. Setup ssh access first: ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk")
            
        #Save all country files:
        for i in arange(len(datalabels)):
            ISO_ = ISET2ISO_country_codes(datalabels[i])
            filename_ =  'ISET_country_' + ISO_ + '.npz'
            np.savez(path + filename_,t=t, L=L[i], Gw=GW[i]/mean(GW[i]), Gs=GS[i]/mean(GS[i]), datetime_offset=datetime_offset, datalabel=ISO_)
            print 'Saved file: ', path + filename_
            sys.stdout.flush()
         
        #Load the relevant file now that it has been created:       
        npzfile = np.load(path + filename)
        print 'Loaded file: ', path + filename
        sys.stdout.flush()
        
    return npzfile['t'], npzfile['L'], npzfile['Gw'], npzfile['Gs'], npzfile['datetime_offset'], npzfile['datalabel']

def valid_ISO(ISO='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    return (ISO in table['ISO'])
    
def ISO2ISET_country_codes(ISO='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    ISET = table['ISET'][find(table['ISO']==ISO)][0]
    
    return ISET

def ISET2ISO_country_codes(ISET='DK',filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    ISO = table['ISO'][find(table['ISET']==ISET)][0]
    
    return ISO

def get_ISET_country_names(filename='ISET2ISO_country_codes.npy',path='./settings/'):

    table = np.load(path+filename)

    print 'Index\tISO\tISET\tName'
    print '====================================='
    for i in arange(len(table)):
        print '{0}\t{1}\t{2}\t{3}'.format(i, table[i]['ISO'], table[i]['ISET'], table[i]['name'])
    sys.stdout.flush()

def save_ISET_country_codes():

	names = ['Austria','Belgium','Bulgaria','Bosnia and Herzegovina','Czech Republic','Switzerland','Germany','Denmark','Spain','France','Finland','Great Britain','Greece','Hungary','Italy','Ireland','Croatia','Luxembourg','Norway','Netherlands','Portugal','Poland','Romania','Sweden','Slovakia','Slovenia','Serbia']
	ISO_codes  = ['AT','BE','BG','BA','CZ','CH','DE','DK','ES','FR','FI','GB','GR','HU','IT','IE','HR','LU','NO','NL','PT','PL','RO','SE','SK','SI','RS']
	ISET_codes  = ['A','B','BG','BH','CZ','Ch','D','DK','ES','F','FIN','GB','GR','H','I','IRL','Kro','Lux','N','NL','P','PL','Ro','S','SK','SLO','SRB']
	
	table = array([empty(len(names))],dtype={'names': ('name', 'ISO', 'ISET'),'formats': ('S30', 'S8', 'S8')})
	table['name'] = names
	table['ISO'] = ISO_codes
	table['ISET'] = ISET_codes
 
	np.save('./settings/ISET2ISO_country_codes2.npy',table[0])


##############
### same again for regions

def get_ISET_region_data(REG='AT',path='./data/'):
    """
    Returns data for a specific region in the ISET data set. 
    
    The function retrives/stores data from/in ./data as default. If
     the data is not available the function 
     attempts to download it from pepsi. This requires ssh access: 
     ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk 
     
    Filenames are: ISET_region_<datalabel>.npz
    
    Returns
    -------
    t: array
        hours numbered from 0 to N-1
    L: array
        load in MW
    Gw: array
        normalized wind power generation
    Gs: array
        normalized solar power generation
    datetime_offset: scalar
        num2date(datetime_offset) yields the date and hour of t=0
    datalabel: string
        region ISO
        
    Example: Returns the West Denmark load, wind and solar power generation
        t, L, Gw, Gs, datetime_offset, datalabel = get_ISET_region_data('DK_W')
    
    """
    
    if not valid_REG(REG):
        sys.exit("Error (44nlksd): No such REG ({0}). For a list of names use get_ISET_region_names().".format(REG))
    
    filename = 'ISET_region_' + REG + '.npz'
    
    try:
        #Load the data file if it exists:
        npzfile = np.load(path + filename)
        print 'Loaded file: ', path + filename
        sys.stdout.flush()
        
    except IOError:
        print 'Datafile does not exist:', path + filename
        print 'Trying to download data from pepsi...'
        sys.stdout.flush()
        try: 
            t, L, GW, GS, datetime_offset, datalabels = get_data_regions(localhost=True);
        except:
            sys.exit("Error (sdf4dz1): Could not connect to pepsi. Setup ssh access first: ssh -L5432:localhost:5432 USERNAME@pepsi.imf.au.dk")
            
        #Save all region files:
        for i in arange(len(datalabels)):
            ISO_ = ISET2ISO_region_codes(datalabels[i])
            filename_ =  'ISET_region_' + ISO_ + '.npz'
            if (ISO_.endswith('off')):
                np.savez(path + filename_,t=t,L=None, Gw=GW[i]/mean(GW[i]), Gs=None, datetime_offset=datetime_offset, datalabel=ISO_)
            else:
                np.savez(path + filename_,t=t, L=L[i], Gw=GW[i]/mean(GW[i]), Gs=GS[i]/mean(GS[i]), datetime_offset=datetime_offset, datalabel=ISO_)
            print 'Saved file: ', path + filename_
            sys.stdout.flush()
         
        #Load the relevant file now that it has been created:       
        npzfile = np.load(path + filename)
        print 'Loaded file: ', path + filename
        sys.stdout.flush()
        
    return npzfile['t'], npzfile['L'], npzfile['Gw'], npzfile['Gs'], npzfile['datetime_offset'], npzfile['datalabel']

def valid_REG(REG='AT',filename='ISET2ISO_region_codes.npy',path='./settings/'):
    table = np.load(path+filename)

    return (REG in table['ISO'])
    
def ISO2ISET_region_codes(REG='AT',filename='ISET2ISO_region_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    ISET = table['ISET'][find(table['ISO']==REG)][0]
    
    return ISET

def ISET2ISO_region_codes(ISET='AT',filename='ISET2ISO_region_codes.npy',path='./settings/'):

    table = np.load(path+filename)
    
    ISO = table['ISO'][find(table['ISET']==ISET)][0]
    
    return ISO

def get_ISET_region_names(filename='ISET2ISO_region_codes.npy',path='./settings/'):

    table = np.load(path+filename)

    print 'Index\tISO\t\tISET\t\tName'
    print '=========================================================================='
    for i in arange(len(table)):
        print '{0}\t{1:10s}\t{2:10s}\t{3}'.format(i, table[i]['ISO'], table[i]['ISET'], table[i]['name'])
    sys.stdout.flush()

def save_ISET_region_codes():

    ISET_codes = ['A', 'B', 'BG', 'BG_off', 'BH_Co', 'BH_Co_off',
     'B_off', 'CZ', 'Ch', 'DK_O', 'DK_O_off', 'DK_W',
     'DK_W_off', 'D_EON_M', 'D_EON_N', 'D_EON_S', 'D_EnBW', 'D_N_off',
     'D_O_off', 'D_RWE', 'D_VET', 'ES_MM_off', 'ES_NAK_off', 'ES_NW',
     'ES_O', 'ES_S', 'ES_SAK_off', 'FIN_N', 'FIN_N_off', 'FIN_S',
     'FIN_S_off', 'F_AEK_off', 'F_AK_off', 'F_MM_off', 'F_NO', 'F_NW',
     'F_SO', 'F_SW', 'GB_N', 'GB_O_off', 'GB_S', 'GB_Sc_off',
     'GB_W_off', 'GR', 'GR_off', 'H', 'IRL', 'IRL_NI',
     'IRL_NI_off', 'I_N', 'I_N_off', 'I_S', 'I_S_off', 'I_Sar',
     'I_Siz', 'Kro', 'Kro_off', 'Lux', 'NL', 'NL_off',
     'N_M', 'N_N', 'N_S', 'N_S_off', 'PL', 'PL_off',
     'P_M', 'P_M_off', 'P_N', 'P_N_off', 'P_S', 'P_S_off',
     'Ro', 'Ro_off', 'SK', 'SLO', 'SRB', 'S_M',
     'S_M_off', 'S_N', 'S_N_off', 'S_S', 'S_S_off']
    ISO_codes = ['AT', 'BE', 'BG', 'BG_off', 'BA', 'BA_off',
     'BE_off', 'CZ', 'CH', 'DK_E', 'DK_E_off', 'DK_W',
     'DK_W_off', 'DE_M', 'DE_NW', 'DE_SE', 'DE_SW', 'DE_NW_off',
     'DE_NE_off', 'DE_W', 'DE_NE', 'ES_SE_off', 'ES_NW_off', 'ES_NW',
     'ES_SE', 'ES_SW', 'ES_SW_off', 'FI_N', 'FI_N_off', 'FI_S',
     'FI_S_off', 'FR_NW_off', 'FR_SW_off', 'FR_SE_off', 'FR_NE', 'FR_NW',
     'FR_SE', 'FR_SW', 'GB_N', 'GB_E_off', 'GB_S', 'GB_N_off',
     'GB_W_off', 'GR', 'GR_off', 'HU', 'IE', 'IE_N',
     'IE_off', 'IT_N', 'IT_N_off', 'IT_S', 'IT_S_off', 'IT_Sar',
     'IT_Sic', 'HR', 'HR_off', 'LU', 'NL', 'NL_off',
     'NO_M', 'NO_N', 'NO_S', 'NO_S_off', 'PL', 'PL_off',
     'PT_M', 'PT_M_off', 'PT_N', 'PT_N_off', 'PT_S', 'PT_S_off',
     'RO', 'RO_off', 'SK', 'SI', 'RS', 'SE_M',
     'SE_M_off', 'SE_N', 'SE_N_off', 'SE_S', 'SE_S_off']
    names = ['Austria', 'Belgium', 'Bulgaria', 'BG offshore',
     'Bosnia and Herzegovina, Montenegro and Albania', 'BA offshore',
     'Belgium offshore', 'Czech Republic', 'Switzerland',
     'Denmark East', 'Denmark East offshore', 'Denmark West',
     'Denmark West offshore', 'Germany Middle', 'Germany North West',
     'Germany South East', 'Germany South West', 'Germany North West offshore',
     'Germany North East offshore', 'Germany West', 'Germany North East',
     'Spain South East offshore', 'Spain North West offshore',
     'Spain North West',
     'Spain South East', 'Spain South West', 'Spain South West offshore',
     'Finland North', 'Finland North offshore', 'Finland South',
     'Finland South offshore', 'France North West offshore',
     'France South West offshore', 'France South East offshore',
     'France North East', 'France North West',
     'France South East', 'France South West', 'Great Britain North',
     'Great Britain East offshore', 'Great Britain South',
     'Great Britain North offshore',
     'Great Britain West offshore', 'Greece', 'Greece offshore',
     'Hungary', 'Ireland', 'Northern Ireland',
     'Ireland offshore', 'Italy North', 'Italy North offshore',
     'Italy South', 'Italy South offshore', 'Italy Sardegna',
     'Italy Sicily', 'Croatia', 'Croatia offshore',
     'Luxembourg', 'Netherlands', 'Netherlands offshore',
     'Norway Middle', 'Norway North', 'Norway South',
     'Norway South offshore', 'Poland', 'Poland offshore',
     'Portugal Middle', 'Portugal Middle offshore', 'Portugal North',
     'Portugal North offshore', 'Portugal South', 'Portugal South offshore',
     'Romania', 'Romania offshore', 'Slovakia',
     'Slovenia', 'Serbia', 'Sweden Middle',
     'Sweden Middle offshore', 'Sweden North', 'Sweden North offshore',
     'Sweden South', 'Sweden South offshore']
	
    table = array([empty(len(names))],dtype={'names': ('name', 'ISO', 'ISET'),'formats': ('S50', 'S10', 'S10')})
    table['name'] = names
    table['ISO'] = ISO_codes
    table['ISET'] = ISET_codes
 
    np.save('./settings/ISET2ISO_region_codes.npy',table[0])
