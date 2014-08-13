#!/usr/bin/env python
'''
20140711 process, use thermal neutron capture cross-section file
'''
import math
import sys
#import random


class capture():
    def __init__(self,Use_natLi = False):
        self.debug = 0
        
        self.xsecfile = '/Users/djaffe/PythonScripts/NEUTRONCAPTURE/thermal_cross_sections_lbl.txt'
        self.xseckey = {'thermal' : 'sigma(0)', 'absorption' : 'sigma(a)'}
        self.col0 = 29  # data after this value
        self.avonum = 6.022e23
        self.vthermal = 2200. * 100. # cm/s
        self.barn = 1.e-24
        # isotope abundances: NOTE THAT Lithium is pure Li6. natLi = 6Li(7.59%), 7Li(92.41%)
        self.abundance = {'H' : {1: 0.999885, 2: 0.000115}, \
                          'C' : {12: 0.9893, 13: 0.0107},\
                          'CU': {63:0.69, 65:0.31},\
                          'GD': {152:0.002, 154:0.022, 155:0.148, 156:0.205, 157:0.157, 158:0.248, 160:0.219},\
                          'W' : {180:0.0012,182:0.265, 183:0.143, 184:0.306, 186:0.284},\
                          'O' : { 16:0.99757, 17:0.00038, 18:0.00205},\
                          'P' : { 17:1.00}, \
                          'N' : { 14:0.99636, 15:0.00364}\
                          }
        if Use_natLi:
            self.abundance['LI'] = {6:0.0759, 7:0.9241}
            print '\n *** USING NATURAL LITHIUM *** \n'
        else:
            self.abundance['LI'] = {6:1.00}
            print '\n *** USING PURE 6LI *** \n'
        self.atomwt = {}
        self.captureProcess = {}
        # normalize abundances. also generate a table of capture processes for each isotope
        for element in self.abundance:
            a = 0.
            for iso in self.abundance[element]:
                a += self.abundance[element][iso]
                isotope = str(iso) + element
                process = 'thermal'
                if isotope=='6LI': process = 'absorption'
                if isotope=='17O': process = 'absorption'
                self.captureProcess[isotope] = process
            if self.debug>1: print 'element',element,'total abundance',a,'before normalizing'
            for iso in self.abundance[element]:
                self.abundance[element][iso] = self.abundance[element][iso]/a
            a = 0.
            aw = 0.
            for iso in self.abundance[element]:
                f = self.abundance[element][iso]
                a += f
                aw += float(iso)*f
            if self.debug>1: print 'element',element,'total abundance',a,'atomic weight',aw,'after normalizing'
            self.atomwt[element] = aw

        # compound[name] =  [ {elemental composition}, density in g/cm3 ]
        # toluene, 5%GdLS, 5%WLS, 10%CuLS, 0.1%GdLS, etc.
        self.complist = []
        self.compound = {'toluene' : [ {'H':8, 'C':7}, 0.87]}
        self.complist.append( 'toluene' )
        self.compound['pseudocumene'] = [ {'H':12, 'C':9}, 0.8761 ]
        self.complist.append('pseudocumene')
        for n in range(10,16+1):
            name = 'LAB_' + str(n)
            self.compound[name] = [ {'H': 5+2*n+1, 'C':6+n}, 0.86 ]
            self.complist.append(name)

        self.compound['LS'] = [ {'LAB_13':0.9999999999}, 0.86 ]
        self.complist.append('LS')

        # 20140813 add Ultima Gold AB 6013309 from PerkinElmer document
        # "LSC in Practice - LSC Cocktails - Elemental Composition"
        # density is from MSDS
        self.compound['UGAB'] = [ { 'C':0.763, 'H':0.097, 'N':0.0005, 'O':0.138, 'P':0.001}, 0.98 ]
        self.complist.append('UGAB')
        
        self.compound['Gd5LS'] = [ {'LS': 0.95, 'GD':0.05}, 1.22]
        self.complist.append('Gd5LS') 
        self.compound['W5LS']  = [ {'LS': 0.95, 'W' :0.05}, 1.79]
        self.complist.append('W5LS')  
        self.compound['Cu10LS']= [ {'LS': 0.90, 'CU':0.10}, 1.68]
        self.complist.append('Cu10LS')
        self.compound['Gd01LS']= [ {'LS': 0.999,'GD':0.001}, 0.87]
        self.complist.append('Gd01LS')

        self.compound['Li01LS']= [ {'LS': 0.999,'LI':0.001}, 0.87]
        self.complist.append('Li01LS')
        self.compound['Li03LS']= [ {'LS': 0.997,'LI':0.003}, 0.87]
        self.complist.append('Li03LS')
        self.compound['Li04LS']= [ {'LS': 0.996,'LI':0.004}, 0.87]
        self.complist.append('Li04LS')

        for fLi in [0.00000001, 0.0006, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]:
            name = 'Li' + '{0:03d}'.format(int(10000.*fLi+.5)) + 'UGAB'
            self.compound[name] = [ {'UGAB':1.0-fLi, 'LI':fLi}, 0.98 ]
            self.complist.append(name)

        for fLi in [0.001, 0.003, 0.004]:
            name = 'Li0' + str(int(1000.*fLi)) + 'Cu10LS'
            self.compound[name] = [ {'LS' : 0.90-fLi, 'CU':0.10, 'LI':fLi}, 1.68]
            self.complist.append(name)

        # identify main isotope of interest for captures
        self.mainIsotope = {}
        for stuff in self.complist:
            if 'Li' in stuff:
                self.mainIsotope[stuff] = '6LI'
            elif 'Gd' in stuff:
                self.mainIsotope[stuff] = 'GD'
            else:
                self.mainIsotope[stuff] = '1H'
        

        return
    def getNumDens(self, name='toluene'):
        '''
        return a dictionary of number densities for each isotope
        '''
        isoNumDens = {}
        if name in self.compound:
            d = self.compound[name][0]    # dictionary 
            rho = self.compound[name][1]   # density g/cm3

            for element in d:
                # if f<1. then f = fraction by weight
                # else f = number of atoms in molecule
                f = float( d[element] )
                if self.debug>0: print 'capture.getNumDens: name',name,'element',element,'f',f
                if f<1. :
                    # mixture: use mass fractions to calculate molecule
                    if element in self.compound:
                        partialIsoNumDens = self.getNumDens(name=element)
                        for isotope in partialIsoNumDens:
                            isoNumDens[isotope] = f * partialIsoNumDens[isotope]
                    else:
                        isodict = self.abundance[element]
                        elementNumDens = f*rho*self.avonum/self.atomwt[element]
                        for iso in isodict:
                            abund = isodict[iso]
                            isotope = str(iso) + element
                            isoNumDens[isotope] = abund * elementNumDens
                            
                else:
                    # molecule: compute molecular weight
                    isodict = self.abundance[element]
                    mw = 0.
                    for e in d: mw += d[e]*self.atomwt[e]
                    if self.debug>0: print 'molecule',name,'molecular wt',mw
                    # calculate element number density, the use isotopic abundances to calculate isotopic
                    # number densities
                    elementNumDens = f*rho*self.avonum/mw
                    for iso in isodict:
                        abund = isodict[iso]
                        isotope = str(iso) + element
                        isoNumDens[isotope] = abund * elementNumDens
                if self.debug>0: print element,isoNumDens
        return isoNumDens
    def getXSEC(self, isotope='1H', xsecType='thermal'):
        '''
        search cross-section file and find and return appropriate cross-section in barns
        no entry will return zero cross-section
        '''
        xsec = 0.
        if xsecType in self.xseckey:
            keyword = self.xseckey[xsecType]
        else:
            return xsec
        
        # append blank character to isotope name to make search unique
        iso = isotope
        if isotope[-1]!=' ' : iso += ' '
        
        f = open(self.xsecfile,'r')
        pastHeader = False
        Found = False
        for line in f:
            if pastHeader:
                if keyword in line and iso in line:
                    if self.debug>0: print 'capture.getXSEC: Parsing line:',line[:-1]
                    word =  line[self.col0:].split()[0]
                    if Found : # deal with case of multiple entries for an isotope
                        if self.debug>0:
                            print 'capture.getXSEC: Already have xsec',xsec,'for isotope',isotope,'xsecType',xsecType
                            print 'capture.getXSEC: Found line:',line[:-1]
                    if '<' in word:
                        xsec = 0
                    else:
                        xsec = float(word)
                        Found = True 
                        
            else:
                if 'Type' in line: pastHeader = True
        f.close()
        return xsec
    def moleFormula(self, name='toluene'):
        '''
        return string with molecular formula
        '''
        formula = ''
        if name in self.compound:
            d = self.compound[name][0]
            for element in d:
                n = d[element]
                if n>=1:
                    formula += element + str(n)
                else:
                    formula += '{0}({1:5.4f}) '.format(element,n)
        else:
            formula = 'UNKNOWN'
        return formula
    def calcCaptureTime(self, name='toluene', mainIsotope=None):
        '''
        calculate capture time given input material
        '''
        isoNumDens = self.getNumDens(name=name)
        width = 0.
        partialWidth = 0.
        for isotope in isoNumDens:
            xs = self.getXSEC(isotope=isotope, xsecType=self.captureProcess[isotope])
            w = isoNumDens[isotope] * xs * self.barn * self.vthermal
            if self.debug>0: print 'calcCaptureTime',name,isotope,w,isoNumDens[isotope]/1e23,xs
            width += w
            if (mainIsotope is not None) and mainIsotope in isotope: partialWidth += w
        captime = 1./width
        return captime, partialWidth/width
    def reportXS(self):
        '''
        prints cross-section table
        '''
        print '{0:8}& {1:}& {2:} \\\\'.format('isotope','cross-section (barns)','capture type')
        il = []
        for isotope in self.captureProcess: il.append(isotope)
        il.sort()
        for isotope in il:
            xs = self.getXSEC(isotope=isotope, xsecType=self.captureProcess[isotope])
            print '{0:8} & {1:8.1f} & {2} \\\\'.format(isotope,xs,self.captureProcess[isotope])
        return
    def reportCaptime(self,pawFile=None,Keys=['Li','UGAB']):
        '''
        print a table of capture times, formula, proton content
        write file for paw if requested
        '''
        HinLSNumDens = self.getNumDens(name='LS')['1H'] # for comparison
        if pawFile is not None: fpaw = open(pawFile,'w')
        
        print '{0:>15} & {2:} &{1:10} & {3} & {4} & {5}\\\\'.format('Compound','Formula','Capture time(microseconds)','proton density(1/cm3/1e23)','protons/protons(LS)','fraction (main isotope)')
        for stuff in self.complist:
            captime,fraction = self.calcCaptureTime(name=stuff, mainIsotope=self.mainIsotope[stuff])
            formula = self.moleFormula(name=stuff)
            HNumDens = self.getNumDens(name=stuff)['1H']
            Hnorm = HNumDens/HinLSNumDens
            print '{0:>15} & {2:6.2f} & {1:>10} & {3:5.3f} & {4:5.3f} & {5:5.3f} ({6})\\\\'.format(stuff,formula,captime*1.e6,HNumDens/1e23, Hnorm, fraction, self.mainIsotope[stuff])
            if pawFile is not None:
                Good = True
                for k in Keys:
                    if k not in stuff: Good = False
                if Good:
                    f = -1.
                    for s in formula.split(' '):
                        if 'LI' in s: f = float(s.split('(')[1].replace(')',''))
                    line = str(f) + ' ' + str(captime*1e6) + ' \n'
                    fpaw.write(line)
        if pawFile is not None:
            print 'capture.reportCaptime: close file',pawFile
            fpaw.close()
        return
if __name__ == '__main__':
    pawDir = '/Users/djaffe/work/paw/PROSPECT/'
    C = capture()
    print ' '
    C.reportXS()
    print ' ' 
    C.reportCaptime(pawFile=pawDir+'6Li_UGAB.vec')

    Cnat = capture(Use_natLi=True)
    print ' '
    Cnat.reportXS()
    print ' '
    Cnat.reportCaptime(pawFile=pawDir+'natLi_UGAB.vec')
