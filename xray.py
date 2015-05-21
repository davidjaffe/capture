#!/usr/bin/env python
'''
20150520 move to CAPTURE/ and modify
20140710 re-purpose for metal-loaded LS studies of photon attenuation length
20140605 investigate feasibility of x-raying TPC to validate xy recon for DarkSide 
'''
import math
import sys
import random


class xray():
    def __init__(self, Experiment='PROSPECT'):
        if Experiment=='DarkSide':
            self.pawdir = '/Users/djaffe/work/paw/DARKSIDE'
        
            self.prefix = '/Users/djaffe/Documents/DarkMatter/Experiments/DarkSide/RandDonTPC_cryo_recovery/G2_TPC/StainlessSteel/'
        
            self.xcom_filenames = ['304L_SS_highest_Fe_content.dat',\
                               '304L_SS_lowest_Fe_content.dat',\
                               'argon.dat']
        elif Experiment=='PROSPECT':
            self.prefix = 'XCOM/'
            self.pawdir = '/Users/djaffe/work/paw/PROSPECT'
            self.xcom_filenames = ['benzene.txt',\
                                   'carbon.txt',\
                                   'cu10percent_toluene.txt',\
                                   'gd5percent_toluene.txt',\
                                   'w5percent_toluene.txt',\
                                   'toluene.txt',\
                                   'water.txt',\
                                   'lead.txt' ]
        else:
            sys.exit('FAIL. Invalid Experiment ' + Experiment)


            
        self.col0 = 6 # first column to consider for meaningful data
        # densities toluene, benzene, Cu, Gd from wikipedia
        self.density = {'304L_SS_solid': 8.03,\
                        '316L_SS_solid': 7.99,\
                        'argon_gas' :1.784e-3, \
                        'argon_liquid': 1.40,\
                        'toluene_liquid': 0.87,\
                        'water_liquid': 1.00,\
                        'benzene_liquid': 0.8765,\
                        'cu10percent_toluene': 0.87*0.90+8.96*0.10,\
                        'gd5percent_toluene': 0.87*0.95+7.90*0.05,\
                        'w5percent_toluene': 0.87*0.95 + 19.25*0.05,\
                        'lead':11.35,
                        'lead_solid':11.35\
                        }

        self.ColumnNames = {'E'     :'PHOTON  ENERGY   (MeV)',\
                            'cohS'  : 'SCATTERING  COHERENT  (cm2/g)',\
                            'incS'  : 'SCATTERING  INCOHER. (cm2/g)',\
                            'peABS' : 'PHOTOELECTRIC ABSORPTION (cm2/g)',\
                            'pairN' : 'PAIR PRODUCTION IN NUCLEAR FIELD (cm2/g)',\
                            'pairE' : 'PAIR PRODUCTION IN NUCLEAR FIELD (cm2/g)',\
                            'attwCS': 'TOTAL ATTENUATION WITH COHERENT SCATT. (cm2/g)',\
                            'attwoCS':'TOTAL ATTENUATION WITHOUT COHERENT SCATT. (cm2/g)'\
                            }
        self.AttenTable = {}
        return
    def getAtten(self,material='lead',state='solid',Egamma=1.):
        '''
        return attenuation factor for gamma of energy Egamm incident on
        specified material.
        Fill table on first call for a new material+state.
        '''
        name = material + '_' + state
        if name not in self.AttenTable: 
            XCOM = self.getXCOM_data(material=material,state=state)
            icol = 6 # total atten with coherent scatt
            Energy = []
            for E in XCOM:
                Energy.append(E)
            Energy.sort()
            Lambda = []
            for i,E in enumerate(Energy):
                lam = XCOM[E][icol]
                Lambda.append(lam)
            self.AttenTable[name] = [Energy,Lambda]
        else:
            Energy = self.AttenTable[name][0]
            Lambda = self.AttenTable[name][1]

        #
        i1,i2 = self.getNeigbors(Egamma,Energy)
        #print 'i1,i2',i1,i2,'Energys',Energy[i1],Energy[i2],'Lambdas',Lambda[i1],Lambda[i2]
        lam = self.yint(Energy[i1],Lambda[i1],Energy[i2],Lambda[i2],Egamma)
        return lam
    def getNeigbors(self,a,X):
        '''
        return indices j,k of X such that X[j]<a<X[k] and k=j or k=j+1
        Assumes X is monotonically increasing
        '''
        j = None
        dist = 1.e20
        for i in range(len(X)):
            Q = abs(X[i]-a)
            if Q>dist:
                break
            if Q<dist:
                j = i
                dist = Q
        if a<X[j]:
            k = max(0,j-1)
        else:
            k = min(len(X)-1,j+1)
        return j,k
    def yint(self,x1,y1,x2,y2,xt):
        dx = x2-x1
        if dx==0. : return y1
        m = (y2-y1)/dx
        b = y2 - m*x2
        return m*xt+b   
    def getXCOM_data(self,material='argon',state='liquid',special=''):
        '''
        determine appropriate file, determine density (g/cc), read file and
        convert attenuation length from g/cm2 to cm
        output is dict with energy as key and a list corresponding to names in self.ColumnNames
        '''
        ####
        print 'getXCOM_data: material',material,'state',state,'special',special
        #### open file
        f = None
        sm = material
        if special!='': sm += '_' + special
        for fn in self.xcom_filenames:
            #print 'sm',sm,'fn',fn
        
            if sm == fn.split('.')[0] :
                f = open(self.prefix + fn,'r')
                print 'getXCOM_data: Opening',fn
                break
        if f is None:
            words = 'getXCOM_data: ERROR No file found. Invalid material '+sm
            sys.exit(words)
        # get density
        ms = material + '_' + state
        rho = None
        for name in self.density:
            if ms.lower()==name.lower():
                rho = self.density[name]
                print 'getXCOM_data: Using density of',rho,'g/cc for',ms
                break
        if rho is None:
            words = 'getXCOM_data: ERROR invalid material + state '+ ms
            sys.exit(words)

        # parse file
        XCOM = {}
        for line in f:
            if line[0:1]!='*':
                #print 'line',line[:-1]
                sl = line[self.col0:].split()
                #print 'split',sl
                if len(sl)>0:
                    data = []
                    for i,w in enumerate(sl):
                        r = 1.
                        if i>0: r = rho # 0th column is energy in MeV, >0th columns are mu/rho
                        data.append(float(w)*r) # convert to 1/cm
                    XCOM[data[0]] = data # fill dict
                    #print 'Energy',data[0],'data',data
        f.close()
        return XCOM
    def buildDet(self,model=0):
        '''
        build a detector using input model number
        Mainly look at effect of changing solid material to LAr
        '''
        print 'buildDet: model',model
        Det = []
        dx0 = 3.0  # cm
        if model==2 or model==3: dx0 = 0.5
        if model>3: dx0 = 10.0
        
        if model==0 or model==2:
            dx = dx0
            material = '304L_SS'
            state = 'solid'
            special = 'lowest_Fe_content'
            XCOM = self.getXCOM_data(material=material,state=state,special=special)
            Det.append( [dx, material, state, special, XCOM] )
        elif model==1 or model==3:
            dx = dx0
            material = 'argon'
            state = 'liquid'
            special = ''
            XCOM = self.getXCOM_data(material=material,state=state,special=special)
            Det.append( [dx, material, state, special, XCOM] )
        elif model==4 or model==5 or model==6 or model==7 or model==8:
            dx = dx0
            if model==4: material = 'gd5percent'
            if model==5: material = 'cu10percent'
            if model==7: material = 'w5percent'
            state = 'toluene'
            special = 'toluene'
            if model==6:
                material = 'toluene'
                state = 'liquid'
                special = ''
            if model==8: material,state,special = 'lead','solid',''
            XCOM = self.getXCOM_data(material=material, state=state,special=special)
            Det.append( [dx, material, state, special, XCOM] )
        else:
            words = 'buildDet: ERROR Invalid model ' + str(model)
            sys.exit(words)
        return Det
    def writeFlux(self, model=0, filename='',attlen_only=False):
        '''
        write file for paw of reduction in flux or attenuation length as function of energy for given model
        '''
        Det = self.buildDet(model=model)
        # use first material in file to get a list of energies, then sort the list
        XCOM = Det[0][4]
        Energies = []
        for E in XCOM:
            Energies.append(E)
        Energies.sort()

        # compute the reduction in flux at each energy for each layer
        flux = {}
        for E in Energies:
            for layer in Det:
                dx = layer[0]
                XCOM = layer[4]
                if E in XCOM:
                    lamb = XCOM[E][6]  # attenuation with coherent scattering
                    if attlen_only:
                        att = 1./lamb
                    else:
                        att = math.exp(-lamb*dx)
                    if E in flux:
                        if attlen_only:
                            print 'writeFlux: Found same energy',E,'twice???'
                        else:
                            flux[E] *= att
                    else:
                        flux[E] = att
                else:
                    print '\ncould not find energy',E,'for',layer,'\n'
                    
        # write to file
        if filename=='':
            fn = self.pawdir + '/model' + str(model).zfill(3) + '.vec'
        else:
            fn = self.pawdir + '/' + filename + '.vec'
        f = open(fn,'w')
        for E in Energies:
            s = str(E) + ' ' +  str(flux[E]) + ' \n'
            f.write(s)
        f.close()
        print 'writeFlux: Wrote',fn
        return 
if __name__ == '__main__':
    Experiment = 'PROSPECT'
    x = xray(Experiment=Experiment)
    if Experiment=='DarkSide':
        w = x.writeFlux(model=0,filename='SS_30mm')
        w = x.writeFlux(model=1,filename='LAr_30mm')
        w = x.writeFlux(model=2,filename='SS_5mm')
        w = x.writeFlux(model=3,filename='LAr_5mm')
    elif Experiment=='PROSPECT':
        #w = x.writeFlux(model=4,filename='gd5pc_toluene',attlen_only=True)
        #w = x.writeFlux(model=7,filename='w5pc_toluene',attlen_only=True)
        #w = x.writeFlux(model=5,filename='cu10pc_toluene',attlen_only=True)
        #w = x.writeFlux(model=6,filename='toluene_only',attlen_only=True)
        w = x.writeFlux(model=8,filename='lead',attlen_only=True)

    else:
      print 'xray.py: Invalid Experiment',Experiment  
    
