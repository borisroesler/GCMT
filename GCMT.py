#!/usr/bin/env python

import numpy as np
import pandas as pd
import calendar
import datetime
import urllib
from bs4 import BeautifulSoup
import re
import obspy


#%% Functions
def updateCatalog(catalog,event):
    catalog = catalog.append({'date': str(obspy.UTCDateTime(event[2].split(': ')[1])).split('T')[0],
                              'time': '%02i:%02i:%05.2f' % (int(event[11][7:11]), int(event[11][12:16]), float(event[11][17:23])),
                              'event name': event[0].strip(),
                              'half duration': float(event[13].split(':')[1]),
                              'latitude': float(event[11][26:33]),
                              'longitude': float(event[11][35:42]),
                              'depth': float(event[11][43:49]),
                              'Mrr': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[18][6:13]),
                              'Mtt': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[18][15:21]),
                              'Mpp': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[18][22:29]),
                              'Mrt': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[18][30:37]),
                              'Mrp': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[18][38:45]),
                              'Mtp': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[18][46:54]),
                              'MrrError': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[19][6:13]),
                              'MttError': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[19][15:21]),
                              'MppError': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[19][22:29]),
                              'MrtError': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[19][30:37]),
                              'MrpError': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[19][38:45]),
                              'MtpError': 10**int(event[16].split('    ')[0].split(': ')[1])*float(event[19][46:54]),
                              }, ignore_index = True)
    return catalog


def makeMatrix(MT):
    MT = np.array([[MT[0], MT[3], MT[4]], [MT[3], MT[1], MT[5]], [MT[4], MT[5], MT[2]]])
    return MT


def unmakeMatrix(MT):
    MT = np.array([MT[0,0], MT[1,1], MT[2,2], MT[0,1], MT[0,2], MT[1,2]])
    return MT


def Moment(M):
    
    trace = np.trace(M)
    M_iso = np.diag(np.array([trace,trace,trace]))
    
    M_devi = M - M_iso
    eigenw, _ = np.linalg.eig(M_devi)
    
    M0 = np.sqrt(np.sum((eigenw**2)/2))
    
    return (M0)


def Magnitude(M0):
    return 2/3*(np.log10(M0)-16.1)


#%% Download
GCMT_catalog = pd.DataFrame(columns = ['date', 'time', 'event name', 'latitude', 'longitude', 'depth', 'Mrr', 'Mtt', 'Mpp', 'Mrt', 'Mrp', 'Mtp', 'MrrError', 'MttError', 'MppError', 'MrtError', 'MrpError', 'MtpError'])
years = np.arange(1976,2022,1).tolist()

counter = 0
for i,year in zip(range(len(years)),years):
    
    if calendar.isleap(year):
        days=366
    else:
        days=365
    
    for j in np.arange(0,days,1):
        start = datetime.datetime(year,1,1)+datetime.timedelta(days=int(j))
        print (str(start).split(' ')[0])
    
        results = 'https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr='+str(start.year)+'&mo='+str(start.month)+'&day='+str(start.day)+'&otype=nd&nday=1&list=5'
        
        try:
            text = urllib.request.urlopen(results)
            webtext = text.read()
            cleantext = BeautifulSoup(webtext, 'lxml').text.strip()
        except urllib.error.URLError:
            cleantext = ''
        
        # days with more than 9 events
        while 'More solutions' in str(webtext):
            link_regex = re.compile('((https?):((//)|(\\\\))+([\w\d:#@%/;$()~_?\+-=\\\.&amp;](#!)?)*)', re.DOTALL)
            link = re.findall(link_regex, str(webtext))[1][0]
        
            cleantext = cleantext[cleantext.find('Event name'):cleantext.find('More solutions')]
            
            number_of_events = cleantext.count('Event name')
            cleantext = cleantext.split('Event name: ')[1:]
            for i in range(number_of_events):
                event = cleantext[i].split('\n')[:-1]
                GCMT_catalog = updateCatalog(GCMT_catalog,event)
    
                print ('added event '+str(event[0].strip()))
            
            try:
                text = urllib.request.urlopen(link)
                webtext = text.read()
                cleantext = BeautifulSoup(webtext, 'lxml').text.strip()
            except urllib.error.URLError:
                cleantext = ''
            
            counter += number_of_events
            
            
        # days with less than 9 events
        cleantext = cleantext[cleantext.find('Event name'):cleantext.find('End of events')]
        
        number_of_events = cleantext.count('Event name')
        cleantext = cleantext.split('Event name: ')[1:]
        for i in range(number_of_events):
            event = cleantext[i].split('\n')[:-1]
            GCMT_catalog = updateCatalog(GCMT_catalog,event)
            
            print ('added event '+str(event[0].strip()))
        
        counter += number_of_events
        
        print ('events in GCMT catalog: ',counter)
        print()

#%% Additional Information
M0 = np.zeros((len(GCMT_catalog)))
Mw = np.zeros((len(GCMT_catalog)))
exponent = np.zeros((len(GCMT_catalog)),dtype=int)
for i in range(len(GCMT_catalog)):
    Mrr = GCMT_catalog.iloc[i]['Mrr']
    Mtt = GCMT_catalog.iloc[i]['Mtt']
    Mpp = GCMT_catalog.iloc[i]['Mpp']
    Mrt = GCMT_catalog.iloc[i]['Mrt']
    Mrp = GCMT_catalog.iloc[i]['Mrp']
    Mtp = GCMT_catalog.iloc[i]['Mtp']
    MT = np.array([Mrr,Mtt,Mpp,Mrt,Mrp,Mtp])
    M0[i] = '{:.4g}'.format(Moment(makeMatrix(MT)))
    Mw[i] = round(Magnitude(M0[i]),2)
    exponent[i] = int(np.log10(M0[i]))

GCMT_catalog['moment'] = M0
GCMT_catalog['magnitude'] = Mw
GCMT_catalog['exponent'] = exponent

#%%
GCMT_catalog = GCMT_catalog.sort_values(by=['date','time'])
GCMT_catalog.reset_index(drop=True, inplace=True)
GCMT_catalog.to_csv('Data/GCMT_catalog.txt', index=False)
