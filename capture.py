#!/usr/bin/env python
'''
20140711 process, use thermal neutron capture cross-section file
'''
import math
import sys
from operator import xor
import os
import xray 
#import random


class capture():
    def __init__(self,Use_natLi = False, xsecOrigin='NIST'):
        self.debug = 0

        if xsecOrigin=='LBL':
            self.xsecfile = '/Users/djaffe/work/GIT/CAPTURE/NEUTRONCAPTURE/thermal_cross_sections_lbl.txt'
            self.col0 = 29  # data after this value in lbl file
            self.xseckey = {'thermal' : ['sigma(0)'], 'absorption' : ['sigma(a)', 'sigma(A)'] } # for lbl file
        if xsecOrigin=='NIST':
            self.xsecfile = '/Users/djaffe/work/GIT/CAPTURE/NEUTRONCAPTURE/thermal_cross_sections_nist.txt'
        self.avonum = 6.022e23
        self.vthermal = 2200. * 100. # cm/s
        self.barn = 1.e-24
        # isotope abundances:
        self.abundance = {'H' : {1: 0.999885, 2: 0.000115}, \
                          'B' : {10: 0.199, 11: 0.801}, \
                          'C' : {12: 0.9893, 13: 0.0107},\
                          'CU': {63:0.69, 65:0.31},\
                          'GD': {152:0.002, 154:0.022, 155:0.148, 156:0.205, 157:0.157, 158:0.248, 160:0.219},\
                          'W' : {180:0.0012,182:0.265, 183:0.143, 184:0.306, 186:0.284},\
                          'O' : { 16:0.99757, 17:0.00038, 18:0.00205},\
                          'P' : { 17:1.00}, \
                          'AL': { 27:1.00}, \
                          'FE': { 54:0.05845, 56:0.91754, 57:0.02119, 58:0.00282},\
                          'CR': { 50:0.04345, 52:0.83789, 53:0.09501, 54:0.02365},\
                          'NI': { 58:0.68077, 60:0.26223, 61:0.011399,62:0.036346},\
                          'N' : { 14:0.99636, 15:0.00364}\
                          }
        self.alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        self.Header = 'Using '
        # either use 6Li-enriched or natural lithium, natLi = 6Li(7.59%), 7Li(92.41%)
        if Use_natLi:
            self.abundance['LI'] = {6:0.0759, 7:0.9241}
            print '\n *** USING NATURAL LITHIUM *** \n'
            self.Header += 'Natural Lithium '
        else:
            self.abundance['LI'] = {6:1.00}
            print '\n *** USING PURE 6LI *** \n'
            self.Header += 'Pure 6Li '
        f = open(self.xsecfile)
        l = f.readline()
        r = l.replace('#',' ')
        self.Header += 'and cross-sections' + r[:-1]
        f.close()
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
                if isotope=='10B': process = 'absorption'
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
        self.compound['HDPE'] = [ {'H':4, 'C':2}, 0.97 ]
        self.complist.append('HDPE')
        self.compound['SS304'] = [ {'CR': 0.19, 'NI': 0.10, 'FE':0.71}, 8.03]
        self.complist.append('SS304')
        for n in range(10,16+1):
            name = 'LAB-' + str(n)
            self.compound[name] = [ {'H': 5+2*n+1, 'C':6+n}, 0.86 ]
            self.complist.append(name)

        self.compound['LS'] = [ {'LAB-13':0.9999999999}, 0.86 ]
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

        self.compound['Al'] = [ {'AL': 1.0}, 2.70]
        self.complist.append('Al')
        self.compound['Iron'] = [ {'FE': 1.0}, 7.87]
        self.complist.append('Iron')

        self.compound['BPE']  = [ {'B': 0.05, 'HDPE':0.95}, 1.01]
        self.complist.append('BPE')


        for nscrew in [6., 30., 60., 120.]:
            screwOD = 3./8.
            sheetL = 96.
            sheetW = 48.
            f = nscrew * math.pi * math.pow(screwOD/2.,2)/sheetL/sheetW
            #print 'SS fraction in BPE sheet',f,'for',nscrew,'bolts with OD',screwOD,'inch'
            suffix = str(int(nscrew)) + 'b'
            self.makeMixture(name1='BPE',name2='SS304',frac1=1.-f,suffix=suffix)
            self.makeMixture(name1='HDPE',name2='SS304',frac1=1.-f,suffix=suffix)
            self.makeMixture(name1='BPE',name2='Al',frac1=1.-f,suffix=suffix)
            self.makeMixture(name1='HDPE',name2='Al',frac1=1.-f,suffix=suffix)

        for fLi in [0.001]: # [0.00000001, 0.0006, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.0073, 0.008, 0.009, 0.010]:
            name = 'Li' + '{0:03d}'.format(int(10000.*fLi+.5)) + 'UGAB'
            self.compound[name] = [ {'UGAB':1.0-fLi, 'LI':fLi}, 0.98 ]
            self.complist.append(name)

        for fLi in [0.001]: #[0.001, 0.003, 0.004]:
            name = 'Li0' + str(int(1000.*fLi)) + 'Cu10LS'
            self.compound[name] = [ {'LS' : 0.90-fLi, 'CU':0.10, 'LI':fLi}, 1.68]
            self.complist.append(name)

        self.complist.sort()

        # identify main isotope of interest for captures
        self.mainIsotope = {}
        for stuff in self.complist:
            if 'Li' in stuff:
                self.mainIsotope[stuff] = '6LI'
            elif 'Gd' in stuff:
                self.mainIsotope[stuff] = 'GD'
            elif 'SS304' in stuff:
                self.mainIsotope[stuff] = 'FE'
            elif 'Al' in stuff:
                self.mainIsotope[stuff] = 'AL'
            elif 'BPE' in stuff:
                self.mainIsotope[stuff] = '10B'
            else:
                self.mainIsotope[stuff] = '1H'

        # for gamma ray attenuation
        self.XX = xray.xray()
        self.XXfilename = 'NEUTRONCAPTURE/XXXX_thermal_neutron_capture_gammas.txt'

        return
    def makeMixture(self, name1=None, name2=None, frac1=0., suffix=None):
        '''
        create a mixture of two compounds given the names and the fraction of the
        first compound.
        name of new compound is the concantenation of the two names and an optional suffix
        '''
        if (name1 is None) or (name2 is None):
            print 'makeMixture: ERROR must specify names of two compounds'
            return
        if frac1<0. or frac1>1.:
            print 'makeMixture: ERROR Invalid fraction',frac1,'of',name1,'must be [0,1]'
            return
        newName = name1 + '-' + name2
        if suffix is not None: newName += '-' + suffix
        frac2 = 1.-frac1
        rho1 = self.compound[name1][1]
        rho2 = self.compound[name2][1]
        rho = frac1*rho1 + frac2*rho2
        self.compound[newName] = [ {name1:frac1, name2:frac2}, rho]
        self.complist.append(newName)
        #print 'makeMixture:',newName,{name1:frac1, name2:frac2},'density',rho
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
                if self.debug>0: print 'capture.getNumDens: name',name,'element',element,'frac by wt',f
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
                if self.debug>0: print 'capture.getNumDens:element',element,'isoNumDens',isoNumDens
        if self.debug>0: print 'capture.getNumDens: name',name,'final isoNumDens',isoNumDens
        return isoNumDens
    def getXSEC(self, inputIsotope='1H', xsecType='thermal'):
        '''
        given inputIsotope search cross-section file and find and return appropriate cross-section in barns.
        For NIST cross-section list, if isotope has 100% abundance, use the element name only in search.
        no entry will return zero cross-section
        '''
        NIST = 'nist' in self.xsecfile
        LBL  = 'lbl' in self.xsecfile
        if not xor(NIST,LBL):
            sys.exit('getXSEC: ERROR. Cannot process xsecfile ' + self.xsecfile)
            
        isotope = inputIsotope
        if NIST:
            element = self.getElement(inputIsotope)
            if len(self.abundance[element])==1:
                #print 'getXSEC: inputIsotope',inputIsotope,'has 100% abundance. Use',element,'in search for cross-section'
                isotope = element
        
        xsec = 0.

        if NIST:
            f = open(self.xsecfile,'r')
            for line in f:
                if line[0]!='#':
                    s = line[:-1].split('\t')
                    #print s
                    if isotope==s[0].upper():
                        cxsec = s[7]
                        if '(' in cxsec: cxsec = cxsec.split('(')[0]
                        xsec = float(cxsec)
                        #print '********************* found',isotope,'cxsec',cxsec,'xsec',xsec,'+++++++++++'
                        break
            f.close()
            return xsec
            
        if LBL:
            if xsecType in self.xseckey:
                keywords = self.xseckey[xsecType]
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
                    foundKey = False
                    for keyword in keywords:
                        if keyword in line: foundKey = True
                    if foundKey and iso in line:
                        if self.debug>0: print 'capture.getXSEC: Parsing line:',line[:-1]
                        word =  line[self.col0:].split()[0]
                        if Found : # deal with case of multiple entries for an isotope
                            if self.debug>0:
                                print 'capture.getXSEC: Already have xsec',xsec,'for isotope',isotope,'xsecType',xsecType
                                print 'capture.getXSEC: Found line:',line[:-1]
                        if '<' in word:
                            xsec += 0
                        else:
                            xsec += float(word)
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
                    formula += '{0}({1:4.3f}) '.format(element,n)
        else:
            formula = 'UNKNOWN'
        return formula
    def calcCaptureTime(self, name='toluene', mainIsotope=None):
        '''
        calculate capture time given input material.
        Also return list of isotopes and normalized partial capture widths
        '''
        isoNumDens = self.getNumDens(name=name)
        width = 0.
        partialWidth = 0.
        isotopes,pW = [],[] # isotopes and their partial widths
        for isotope in isoNumDens:
            xs = self.getXSEC(inputIsotope=isotope, xsecType=self.captureProcess[isotope])
            w = isoNumDens[isotope] * xs * self.barn * self.vthermal
            if self.debug>0: print 'calcCaptureTime:name',name,'isotope',isotope,'partialwidth',w,'isotopeNumDens/1e23',isoNumDens[isotope]/1e23,'xsec',xs
            width += w
            isotopes.append(isotope)
            pW.append(w)
            if (mainIsotope is not None) and mainIsotope in isotope: partialWidth += w
        captime = 1./width
        if self.debug>0: print 'calcCaptureTime:name',name,'width',width,'partialWidth',partialWidth,'(to',mainIsotope,'),captime',captime
        for i,w in enumerate(pW): pW[i]=w/width
        
        if self.debug>0:
            print 'isotope,pW',
            for i,w in zip(isotopes,pW): print i,w,
            print ''
        return captime, partialWidth/width, isotopes,pW
    def getAN(self,el):
        '''
        return integer corresponding to atomic number for input element of style NNXX
        '''
        AN = -1
        Done = False
        i = -1
        while not Done:
            try:
                AN = int(el[:i])
            except ValueError:
                i -= 1
            else:
                Done = True
        return AN
    def getElement(self,el):
        '''
        return element XX corresponding to input element of style NNXX
        '''
        element = ''
        for i in range(len(el)):
            if el[i] in self.alphabet: element += el[i]
        return element
    def reportXS(self):
        '''
        prints cross-section table
        '''
        words = 'The mean capture gamma energy and the energy of the strongest capture gamma, obtained from http://www.nndc.bnl.gov/capgam/, are also listed. A single energy of 478 keV for 10B is entered by hand.'
        print '\caption{Thermal cross-sections for isotopes considered in this note. '+words+self.Header+'}'
        print '\\begin{tabular}{r|r|l|r|r}'
        print '{0:8}& {1:}& {2:} & {3} & {4} \\\\ \hline'.format('isotope','$\sigma$ (barns)','capture type', 'mean $E^{\\rm cap}_{\gamma}$(keV)','strongest $E_\gamma$(keV)')
        il = []
        for isotope in self.captureProcess: il.append(isotope)
        il.sort(key=lambda ql: self.getAN(ql))
        for isotope in il:
            xs = self.getXSEC(inputIsotope=isotope, xsecType=self.captureProcess[isotope])
            niceIN = self.niceIsotopeName(isotope)
            iso,Emean,Estrongest = self.processCapGam(filename=self.XXfilename.replace('XXXX',niceIN),units='keV')
            print '{0:8} & {1:8.1f} & {2} & {3:5.1f} & {4:5.1f} \\\\ [-1.5ex]'.format(isotope,xs,self.captureProcess[isotope],Emean,Estrongest)
        print '\\end{tabular} \n'
        return
    def niceIsotopeName(self,isotope):
        '''
        nice character string for isotope name.e.g., 27AL is rendered as 27Al
        '''
        an = str(self.getAN(isotope))
        el = self.getElement(isotope)
        if len(el)==2: el = el[0]+el[1].lower()
        return an+el
    def calcRelativeGammaFlux(self,isotopes=['1H'],fractions=[1.],tlead=3.,option='meanEnergy'):
        '''
        calculate relative capture gamma flux given list of isotopes and fractions of gammas from each isotope
        assuming a thickness of lead tlead in cm
        '''
        flux = 0.
        for isotope,f in zip(isotopes,fractions):
            niceIN = self.niceIsotopeName(isotope)
            filename = self.XXfilename.replace('XXXX',niceIN)
            iso,Emean,Estrongest = self.processCapGam(filename=filename,units='MeV')
            Egamma = Estrongest
            if option=='meanEnergy': Egamma = Emean
            a = self.XX.getAtten(material='lead',state='solid',Egamma=Egamma)
            alpha = math.exp(-tlead*a)
            #print 'calcRelativeGammaFlux',an+el,'Eg',Emean,'lambda',a,'atten',alpha,'fraction',f
            flux += f*alpha
        return flux
    def reportCaptime(self,pawFile=None,Keys=['Li','UGAB']):
        '''
        print a table of capture times, formula, proton content
        write file for paw if requested
        '''
        tlead = 3. # cm
        
        HinLSNumDens = self.getNumDens(name='LS')['1H'] # for comparison
        c,f,isotopes,pW = self.calcCaptureTime(name='BPE',mainIsotope=self.mainIsotope['BPE'])
        BPErelGF_Emean = self.calcRelativeGammaFlux(isotopes=isotopes,fractions=pW,tlead=tlead,option='meanEnergy')
        BPErelGF_Estrong=self.calcRelativeGammaFlux(isotopes=isotopes,fractions=pW,tlead=tlead,option='strongestEnergy')
        if pawFile is not None: fpaw = open(pawFile,'w')

        words = 'fluxA uses the mean energy of capture gammas and fluxB uses the strongest gamma capture energy.'
        print '\caption{Calculated thermal neutron capture times, formula, proton density(${\\rm cm}^{-3} 10^{-23}$), proton fraction compared to LS and capture fraction on a selected isotope for each compound. The capture gamma flux relative to BPE is also given. '+words+self.Header+'}'
        print '\\begin{tabular}{l|r|c|r|r|r|r|r}'
        print '{0:>15} & {2:} &{1:10} & {3} & {4} & {5} & {6} & {7} \\\\ [-1.5ex]'.format(' ',' ','$\\tau_{\\rm Capture}$','',' ','fraction','Relative','Relative')
        print '{0:>15} & {2:} &{1:10} & {3} & {4} & {5} & {6} & {7}\\\\ \hline'.format('Compound','Formula','($\\mu$s)','$\\rho(p)$','$p/p({\\rm LS)}$','(main isotope)','$\\gamma$ fluxA','$\\gamma$ fluxB')
        for stuff in self.complist:
            captime,fraction,isotopes,pW = self.calcCaptureTime(name=stuff, mainIsotope=self.mainIsotope[stuff])
            relGamFlux_Emean = self.calcRelativeGammaFlux(isotopes=isotopes,fractions=pW,tlead=tlead,option='meanEnergy')
            relGamFlux_Estrong=self.calcRelativeGammaFlux(isotopes=isotopes,fractions=pW,tlead=tlead,option='strongestEnergy')
            #print stuff,'relativeGammaFlux',relGamFlux,'Relative to BPE',relGamFlux/BPErelGF
            formula = self.moleFormula(name=stuff)
            HNumDens = 0.
            if '1H' in self.getNumDens(name=stuff): HNumDens = self.getNumDens(name=stuff)['1H']
            Hnorm = HNumDens/HinLSNumDens
            GnormA = relGamFlux_Emean/BPErelGF_Emean
            GnormB = relGamFlux_Estrong/BPErelGF_Estrong
            print '{0:>15} & {2:6.2f} & {1:>10} & {3:5.3g} & {4:5.3f} & {5:0.2g} ({6}) & {7:6.4f}&{8:6.4f}\\\\ [-1.5ex]'.format(stuff,formula,captime*1.e6,HNumDens/1e23, Hnorm, fraction, self.mainIsotope[stuff], GnormA, GnormB)
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
        print '\\end{tabular}'
        if pawFile is not None:
            print 'capture.reportCaptime: close file',pawFile
            fpaw.close()
        return
    def processCapGam(self,filename='NEUTRONCAPTURE/27Al_thermal_neutron_capture_gammas.txt',units='MeV'):
        '''
        return target nucleus name,
        weighted mean capture gamma energy and most intense capture gamma energy.
        Note that energies in input file are in keV.
        '''
        
        maxStrength = -1.
        EmaxStrength= 0.
        wtave = 0.
        wts   = 0.
        tgtNuc = None
        if os.path.isfile(filename):
            f = open(filename,'r')
            for line in f:
                if 'Target Nucleus=' in line:
                    tgtNuc = line.split('=')[1]
                    tgtNuc = tgtNuc.replace('\n','')
                s = line.split('\t')
                #print s[0]
                try:
                    Egam = float(s[0])
                    #print 'accepted',s
                    strength = 0.
                    if len(s)>4: strength = float(s[4])
                    if strength>maxStrength:
                        maxStrength = strength
                        EmaxStrength= Egam
                    #print 'Egam',Egam,'strength',strength,'maxStrength',maxStrength,'Emax',EmaxStrength

                    wts += strength
                    wtave += Egam*strength
                except ValueError:
                    continue
                    #print 'rejected',line[:-1]
            if wts>0.: wtave = wtave/wts
            if tgtNuc not in filename:
                print 'processCapGam: ERROR tgt nucleus',tgtNuc,'inconsistent with filename',filename
                sys.exit('processCapGam: INCONSISTENT FILENAME AND TARGET NUCLEUS')
            f.close()
        if units.lower()=='mev': # change energies to MeV?
            wtave /= 1000.
            EmaxStrength /= 1000.
        return tgtNuc,wtave,EmaxStrength
if __name__ == '__main__':
    ### TESTING ####
    if 0: # check lead attenuation length estimate
        X = xray.xray(Experiment='PROSPECT')
        tlead = 3. # cm
        f10B = X.getAtten(material='lead',state='solid',Egamma= 0.478)
        alpha_10B = math.exp(-tlead*f10B)
        for Eg in [0.478, 5.362, 7.631]:
            att = X.getAtten(material='lead',state='solid',Egamma=Eg)
            alpha = math.exp(-tlead*att)
            print 'Eg',Eg,'att',att,'1/att',1/att,'alpha',alpha,'alpha/alpha(10B)',alpha/alpha_10B

    #### TESTING ####
    if 0:   # check processing of all capture gamma files
        C = capture()
        dn = 'NEUTRONCAPTURE'
        fileList = []
        for name in os.listdir(dn):
            if 'thermal_neutron_capture_gammas.txt' in name:
                fileList.append(  dn + '/' + name )
        for fileName in fileList:
            tgt,Eave,Emax = C.processCapGam(fileName)
            print tgt,'Eave',Eave,'Emax',Emax
        sys.exit("END OF TESTING")

    
    # got to change next line if you want to write files to your machine
    pawDir = '/Users/djaffe/work/paw/PROSPECT/'
    ListOfXsec = ['NIST'] # 'LBL' is another
    for xsecOrigin in ListOfXsec:
        # first do calculations with enriched Li (pure 6Li)
        C = capture(xsecOrigin=xsecOrigin)
        print ' '
        C.reportXS()
        print ' ' 
        C.reportCaptime() #pawFile=pawDir+'6Li_UGAB.vec')

        # now do calculations with natural Li
        Cnat = capture(Use_natLi=True,xsecOrigin=xsecOrigin)
        print ' '
        Cnat.reportXS()
        print ' '
        Cnat.reportCaptime() #pawFile=pawDir+'natLi_UGAB.vec')
