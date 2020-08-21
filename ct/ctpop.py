#! /usr/bin/env python3

import os
import sys
import time
import array
import cantera as ct

class CanteraOutput(object):
    
    def __init__(self, datfile):
        self.datfile = datfile
        dirname = os.path.dirname(datfile)
        if dirname != "": dirname = dirname + "/"
        self.info = []
        if not os.path.exists(datfile):
            raise IOError('"%s" does not exists.' % (datfile,))
        
        self.info.append("Data file: %s" % datfile)
        self.info.append("Last modified: %s" % 
                         (time.strftime('%Y-%m-%d %H:%M:%S',
                                        time.localtime(os.path.getmtime(datfile))),))
        fp = open(datfile)
        ll = next(fp).split()
        typ = ll[0]
        outname = ll[1]
        idtag = ll[2]
        
        self.sens = False
        self.sens_hp = False
        if len(ll) > 3:
            if ll[3].lower() == "sens":
                self.sens = True
            elif ll[3].lower() == "sens_hp":
                self.sens = True
                self.sens_hp = True
                
        ctmlfile = dirname + outname + ".xml"
        if not os.path.exists(ctmlfile):
            raise IOError('"%s" does not exists.' % (ctmlfile,))
        self.info.append("CTML file: %s" % ctmlfile)
        self.info.append("id: %s" % idtag)
        self.info.append("Data type: %s" % typ)

        # read propaties
        self.gas = ct.Solution(ctmlfile)
        self.ne = self.gas.n_elements  #mm
        self.ns = self.gas.n_species  #kk
        self.nr = self.gas.n_reactions  #ii
        self.info.append("Number of elements: %s" % self.ne)
        self.info.append("Number of species: %s" % self.ns)
        self.info.append("Number of reactions: %s" % self.nr)
        self.elements = self.gas.element_names
        self.species = self.gas.species_names
        self.reactions = self.gas.reaction_equations(list(range(self.nr)))
        self.weights = [float(x) for x in self.gas.molecular_weights]
        self.compositions = []
        for x in self.species:
            self.compositions.append([self.gas.n_atoms(x, el) for el in self.elements])

        # header
        ll = next(fp).split()
        self.xvarname = ll[0]
        self.addvnames = []
###        if self.sens: self.sens_vars = []
        if self.sens:
            if self.sens_hp: self.sens_vars = ["mass", "T(K)"]
            else: self.sens_vars = ["mass", "V(m3)", "T(K)"]
        
        itemp, ipres, iaddv = None, None, []
        for i in range(len(ll) - 1 - self.ns):
            x = ll[1+i]
            if x == "T(K)": itemp = i+1
            elif x == "P(Pa)": ipres = i+1
            else:
                self.addvnames.append(x)
                iaddv.append(i+1)
###            if self.sens: self.sens_vars.append(x)
        self.naddv = len(iaddv)
        if self.sens: self.sens_vars = self.sens_vars + list(self.species)
        if itemp is None or ipres is None: raise ValueError("wrong header")

        # profiles
        self.xvars = []
        self.press = []
        self.temps = []
        self.fracs = [[] for i in range(self.ns)]
        self.addvs = [[] for i in range(self.naddv)]
        for l in fp:
            vals = list(map(float, l.split()))
            self.xvars.append(vals[0])
            self.press.append(vals[ipres])
            self.temps.append(vals[itemp])
            for i in range(self.naddv):
                self.addvs[i].append(vals[iaddv[i]])
            for i in range(self.ns):
                self.fracs[i].append(vals[i+1+2+self.naddv])
            
        self.ndata = len(self.xvars)
        self.info.append("Number of dataset: %s" % (self.ndata))

        if self.sens:
            sensfile = dirname + outname + ".sens"
            if os.path.exists(sensfile):
                self.sfp = open(sensfile, "rb")
                self.info.append("Sensitivity analysis: %s" % sensfile)
            else:
                self.info.append("Sensitivity analysis: %s [File not found]" % sensfile)
                self.sens = False

    def ckqxp(self, p, temp, x):
        """Returns the rates of progress (mole/cm3-s) for the reactions given pressure (Pa),
        temperature (K) and mole fractions
        """
        self.gas.TPX = temp, p, x
        rop = self.gas.net_rates_of_progress
        for i in range(len(rop)):
            rop[i] *= 1e-3 # kmole/cm3-s => mole/cm3-s
        return rop

    def ckcont(self, k, q):
        """Returns the contributions (mole/cm3-s) of the reactions to the molar
        production rate of a k-th species.
        q is rates of progress for the reactions (mole/cm3-s)
        """
        cont = []
        for i in range(self.nr):
            rate = q[i]
            r, p = self.gas.reactant_stoich_coeff(k, i), self.gas.product_stoich_coeff(k, i)
            cont.append(- r*rate + p*rate)
        return cont

    def ckconth(self, q):
        """Returns the contributions (J/cm3-s) of the reactions to enthalpy
        formation rate.
        q is rates of progress for the reactions (mole/cm3-s)
        """
        cont = []
        for i in range(self.nr):
            cont.append(-self.gas.delta_enthalpy[i]*q[i]*1.0e-3)
        return cont

    def ckcontu(self, q, temp):
        """Returns the contributions (J/cm3-s) of the reactions to internal
        energy formation rate.
        q is rates of progress for the reactions (mole/cm3-s)
        """
        cont = []
        RT = ct.gas_constant*temp
        for i in range(self.nr):
            dn = 0
            #for x in self.gas.Reaction(i):
            #    for species in x.reactants:
            #        dn += -x.reactants[species]
            #    for species in x.products:
            #        dn += x.products[species]
            for j in range(self.ns):
                dn += -self.gas.reactant_stoich_coeff(j,i) + self.gas.product_stoich_coeff(j,i)
            cont.append(-(self.gas.delta_enthalpy[i]-dn*RT)*q[i]*1.0e-3)
        return cont
  
    def ckinu(self, i):
        """Returns the species indices in a reaction, and their indices
        and stoichiometric coefficients
        retruns reactants{}, products{} (key: species index, val: stoichiometric coeffs)
        """
        rd, pd = {}, {}
        for k in range(self.ns):
            r, p = self.gas.reactant_stoich_coeff(k, i), self.gas.product_stoich_coeff(k, i)
            if abs(r) > 1e-32: rd[k] = r
            if abs(p) > 1e-32: pd[k] = p
        return rd, pd

    # util methods
    def xvarAt(self, index):
        return self.xvars[index]
    
    def indexAt(self, xvar):
        """return index i with which xvar[i] <= xvar < xvar[i+1]"""
        l = self.xvars
        # binary search
        imin, imax = 0, len(l) - 1
        while imin < imax:
            imid = int((imin + imax) / 2)
            if xvar < l[imid]: imax = imid
            else: imin = imid + 1
        if imin > 0: imin -= 1
        return imin

    # profiles

    def getProfileAt(self, index, conc_unit=None):
        """return xvar, press, temp, frac[], addv[] at i-th solution"""
        xvar, press, temp = self.xvars[index], self.press[index], self.temps[index]
        frac = [self.fracs[j][index] for j in range(self.ns)]
        addv = [self.addvs[j][index] for j in range(self.naddv)]
        
        if conc_unit is None: conc_unit = "molefrac"
        if conc_unit == 'molefrac': pass
        elif conc_unit == 'massfrac':
            sumymw = sum([frac[j] * self.weights[j] for j in range(self.ns)])
            for j in range(self.ns):
                frac[j] = frac[j] * self.weights[j] / sumymw
        elif conc_unit == 'molecules/cm3':
            conv = press / (temp * 8.314510e6) * 6.022137e23
            for j in range(self.ns):
                frac[j] = frac[j] * conv
        elif conc_unit == 'mole/cm3':
            conv = press / (temp * 8.314510e6)
            for j in range(self.ns):
                frac[j] = frac[j] * conv
        else:
            raise ValueError('Unsupported concentration unit: %s' % (conc_unit,))

        return xvar, press, temp, frac, addv
        
    def getXvarList(self):
        return self.xvars
        
    def getPressList(self):
        return self.press

    def getTempList(self):
        return self.temps

    def getAddvsList(self):
        return self.addvs

    def getAddvList(self, addv):
        ind = self.addvnames.index(addv)
        return self.addvs[ind]

    def getConcLists(self, species=None, conc_unit=None):
        """species: None for all"""
        if conc_unit is None: conc_unit = "molefrac"
        
        clists = [[] for x in range(self.ns)]
        for i in range(len(self.xvars)):
            if conc_unit == 'molefrac':
                for j in range(self.ns):
                    clists[j].append(self.fracs[j][i])
            elif conc_unit == 'massfrac':
                sumymw = sum([self.fracs[j][i] * self.weights[j] for j in range(self.ns)])
                for j in range(self.ns):
                    clists[j].append(self.fracs[j][i] * self.weights[j] / sumymw)
            elif conc_unit == 'molecules/cm3':
                press, temp = self.press[i], self.temps[i]
                conv = press / (temp * 8.314510e6) * 6.022137e23
                for j in range(self.ns):
                    clists[j].append(self.fracs[j][i] * conv)
            elif conc_unit == 'mole/cm3':
                press, temp = self.press[i], self.temps[i]
                conv = press / (temp * 8.314510e6)
                for j in range(self.ns):
                    clists[j].append(self.fracs[j][i] * conv)
            else:
                raise ValueError('Unsupported concentration unit: %s' % (conc_unit,))
        if species is None:
            return clists
        else:
            return [clists[self.species.index(name)] for name in species]


    def exportProfiles(self, fp, xvar=True, press=True, temp=True, addv=True, species=None,
                       conc_unit=None, precision=6, sep=' ', linesep='\n', header='# '):
        """export xvar profiles.
        species: None for all, False for none
        header: give string to add header with prefix, None for not to add
        """
        headers = []
        exports = []
        fmtstr = '%%.%de' % precision
        if xvar: 
            headers.append(self.xvarname)
            exports.append(self.getXvarList())
        if press:
            headers.append('P(Pa)')
            exports.append(self.getPressList())
        if temp:
            headers.append('T(K)')
            exports.append(self.getTempList())
        if addv:
            headers.extend(self.addvnames)
            exports.extend(self.getAddvsList())
        if species is not False:
            sphdrs = self.species[:] if species is None else species[:]
            sphdrs = ['[%s]' % (s,) for s in sphdrs]
            sphdrs[0] = '%s(%s)' % (sphdrs[0], conc_unit or "molefrac")
            headers.extend(sphdrs)
            exports.extend(self.getConcLists(species, conc_unit))
        if header is not None:
            fp.write('%s%s%s' % (header, sep.join(headers),linesep))
        for x in zip(*exports):
            fp.write('%s%s' % (sep.join([fmtstr % (v,) for v in x]),linesep))

    # production rates (reaction contributions)

    def rateUnitConv(self, rates, unit=None):
        if unit is None: unit = 'molecules/cm3-s'
        if unit == 'molecules/cm3-s': return [x*6.022137e23 for x in rates]
        elif unit == 'mole/cm3-s': return rates[:]
        else: raise ValueError('Unsupported rate unit: %s' % (unit,))
    
    def getProductionRatesOfAt(self, species_name, index, unit=None):
        """contributions of the reactions at xvar[index]
        return totalrate, contributions[nr]"""
        xvar, p, temp, x, addv = self.getProfileAt(index)
        rates = self.ckqxp(p, temp, x)
        cntrbs = self.ckcont(self.species.index(species_name), rates)
        cntrbs =  self.rateUnitConv(cntrbs, unit)
        totrate = sum(cntrbs)
        return totrate, cntrbs
    
    def getHeatProductionRatesOfAt(self, heat_name, index):
        """contributions of the reactions at xvar[index]
        return totalrate, contributions[nr]"""
        xvar, p, temp, x, addv = self.getProfileAt(index)
        rates = self.ckqxp(p, temp, x)
        if heat_name is 'u':
            cntrbs = self.ckcontu(rates, temp)
        elif heat_name is 'h':
            cntrbs = self.ckconth(rates)
        totrate = sum(cntrbs)
        return totrate, cntrbs
    
    def getProductionRatesProfileOf(self, species_name, indices=None, unit=None):
        """contribution profiles of species_name for reactions specified by indices
        return totrates, cntrb_lists"""
        if indices is None:
            indices = list(range(self.nr))
        
        cntrb_lists = [[] for i in range(len(indices))]
        totrates = []
        k = self.species.index(species_name)
        for i in range(len(self.temps)):
            rates = self.ckqxp(self.press[i], self.temps[i],
                               [self.fracs[j][i] for j in range(self.ns)])
            cntrbs = self.rateUnitConv(self.ckcont(k, rates), unit)
            totrates.append(sum(cntrbs))
            for j, x in enumerate(indices):
                cntrb_lists[j].append(cntrbs[x])
        return totrates, cntrb_lists
    
    def getHeatProductionRatesProfileOf(self, heat_name, indices=None):
        """contribution profiles of heat for reactions specified by indices
        return totrates, cntrb_lists"""
        if indices is None:
            indices = list(range(self.nr))
        
        cntrb_lists = [[] for i in range(len(indices))]
        totrates = []
        if heat_name is 'u':
            for i in range(len(self.temps)):
                rates = self.ckqxp(self.press[i], self.temps[i],
                                   [self.fracs[j][i] for j in range(self.ns)])
                cntrbs = self.ckcontu(rates, self.temps[i])
                totrates.append(sum(cntrbs))
                for j, x in enumerate(indices):
                    cntrb_lists[j].append(cntrbs[x])
        elif heat_name is 'h':
            for i in range(len(self.temps)):
                rates = self.ckqxp(self.press[i], self.temps[i],
                                   [self.fracs[j][i] for j in range(self.ns)])
                cntrbs = self.ckcontu(rates, self.temps[i])
                totrates.append(sum(cntrbs))
                for j, x in enumerate(indices):
                    cntrb_lists[j].append(cntrbs[x])
        return totrates, cntrb_lists
    
    def exportProductionRatesOfAt(self, fp, species_name, index, rate_unit=None, 
                                  sort=False, precision=6):
        fmtstr = '%% .%de' % precision
        totrate, cntrbs = self.getProductionRatesOfAt(species_name, index, unit=rate_unit)
        crate = sum([x for x in cntrbs if x >= 0])
        drate = sum([x for x in cntrbs if x < 0])
        xvarstr = fmtstr % (self.xvarAt(index),)
        
        d = {}
        for i,x in enumerate(cntrbs): d[i] = x
        if sort: 
            items = sorted(d.items(), key=lambda x: abs(x[1]), reverse=True)
            fp.write('Production rates (sotred) ')
        else:
            items = d.items()
            fp.write('Production rates ')
        fp.write('of %s at %s sec (index=%d) in %s\n' %
                 (species_name, xvarstr, index, 
                  rate_unit or "molecules/cm3-s"))
        fp.write(' total = %s (creation = %s, destruction = %s)\n' %
                 (fmtstr % (totrate,), fmtstr % (crate,), fmtstr % (drate,)))
        fp.write('  i, rate, Reaction\n')
        for i,x in items:
            r = fmtstr % (x,)
            if x >= 0: s = '(%5.2f%% of creaction  )' % (100.0*x/crate,)
            else: s = '(%5.2f%% of destruction)' % (100.0*x/drate,)
            fp.write('%4d, %s, %s, %s\n' % (i+1, r, s, self.reactions[i]))
    
    def exportHeatProductionRatesOfAt(self, fp, heat_name, index, 
                                  sort=False, precision=6):
        fmtstr = '%% .%de' % precision
        totrate, cntrbs = self.getHeatProductionRatesOfAt(heat_name, index)
        crate = sum([x for x in cntrbs if x >= 0])
        drate = sum([x for x in cntrbs if x < 0])
        xvarstr = fmtstr % (self.xvarAt(index),)
        
        d = {}
        for i,x in enumerate(cntrbs): d[i] = x
        if sort: 
            items = sorted(d.items(), key=lambda x: abs(x[1]), reverse=True)
            fp.write('Heat (%s) production rates (sotred) ' % heat_name)
        else:
            items = d.items()
            fp.write('Heat (%s) production rates ' % heat_name)
        fp.write('of %s at %s sec (index=%d) in J/cm3-s\n' %
                 (heat_name, xvarstr, index,))
        fp.write(' total = %s (creation = %s, destruction = %s)\n' %
                 (fmtstr % (totrate,), fmtstr % (crate,), fmtstr % (drate,)))
        fp.write('  i, rate, Reaction\n')
        for i,x in items:
            r = fmtstr % (x,)
            if x >= 0: s = '(%5.2f%% of creaction  )' % (100.0*x/crate,)
            else: s = '(%5.2f%% of destruction)' % (100.0*x/drate,)
            fp.write('%4d, %s, %s, %s\n' % (i+1, r, s, self.reactions[i]))
    
    def exportProductionRatesProfileOf(self, fp, species_name, indices=None,
                                       rate_unit=None,
                                       precision=6, sep=' ', linesep='\n', header='# '):
        headers = []
        exports = []
        fmtstr = '%% .%de' % precision
        totrates, cntrb_lists = self.getProductionRatesProfileOf(species_name, indices,
                                                                 unit=rate_unit)
        headers.append(self.xvarname)
        exports.append(self.getXvarList())
        headers.append('total(%s)' % (rate_unit or "molecules/cm3-s",))
        exports.append(totrates)
        
        idhdrs = list(range(self.nr)) if indices is None else indices[:]
        idhdrs = ['R%d' % (i+1,) for i in idhdrs]
        headers.extend(idhdrs)
        exports.extend(cntrb_lists)
        if header is not None:
            fp.write('%s%s%s' % (header, sep.join(headers),linesep))
        for x in zip(*exports):
            fp.write('%s%s' % (sep.join([fmtstr % (v,) for v in x]),linesep))

    def exportHeatProductionRatesProfileOf(self, fp, heat_name, indices=None,
                                       precision=6, sep=' ', linesep='\n', header='# '):
        headers = []
        exports = []
        fmtstr = '%% .%de' % precision
        totrates, cntrb_lists = self.getHeatProductionRatesProfileOf(heat_name, indices)
        headers.append(self.xvarname)
        exports.append(self.getXvarList())
        headers.append('total(J/cm3-s)')
        exports.append(totrates)
        
        idhdrs = list(range(self.nr)) if indices is None else indices[:]
        idhdrs = ['R%d' % (i+1,) for i in idhdrs]
        headers.extend(idhdrs)
        exports.extend(cntrb_lists)
        if header is not None:
            fp.write('%s%s%s' % (header, sep.join(headers),linesep))
        for x in zip(*exports):
            fp.write('%s%s' % (sep.join([fmtstr % (v,) for v in x]),linesep))

    # reaction rates
        
    def getReactionRatesAt(self, index, unit=None):
        """net reaction rates (rates of progress) for the reactions at xvar[index]
        return q[nr]"""
        xvar, p, temp, x, addv = self.getProfileAt(index)
        rates = self.ckqxp(p, temp, x)
        return self.rateUnitConv(rates, unit)

    # sens
    
    def getSensCoeffsOfAt(self, ivar, index):
        """sens coeffs for ivar-th variable at xvar[index]"""
        dsize = array.array('d').itemsize
        nvar = len(self.sens_vars)
        self.sfp.seek(index * self.nr * nvar * dsize, 0)
        arr = array.array("d")
        for ir in range(self.nr):
            if ir == 0: self.sfp.seek(ivar * dsize, 1)
            else: self.sfp.seek((nvar - 1) * dsize, 1)
            arr.fromfile(self.sfp, 1)
        return arr.tolist()
    
    def getSensCoeffsProfileOf(self, ivar, indices=None):
        """sens coeffs profiles for ivar-th variable to reactions specified by indices"""
        dsize = array.array('d').itemsize
        nvar = len(self.sens_vars)
        if indices is None:
            indices = list(range(self.nr))
        sens_lists = [[] for i in range(len(indices))]
        for i in range(self.ndata):
            sens = self.getSensCoeffsOfAt(ivar, i)
            for j, x in enumerate(indices):
                sens_lists[j].append(sens[x])
        return sens_lists

    def exportSensCoeffsOfAt(self, fp, ivar, index, sort=False, precision=6):
        fmtstr = '%% .%de' % precision
        sens_list = self.getSensCoeffsOfAt(ivar, index)
        xvarstr = fmtstr % (self.xvarAt(index),)
        
        name = self.sens_vars[ivar]
        fp.write('Sensitivity coefficients ')
            
        d = {}
        for i,x in enumerate(sens_list): d[i] = x
        if sort: 
            items = sorted(d.items(), key=lambda x: abs(x[1]), reverse=True)
            fp.write('(sotred) ')
        else:
            items = d.items()
        fp.write('to variable "%s" at %s sec (index=%d)\n' % (name, xvarstr, index ))
        fp.write('  i, S, Reaction\n')
        for i,x in items:
            s = fmtstr % (x,)
            fp.write('%4d, %s, %s\n' % (i+1, s, self.reactions[i]))
        
    def exportSensCoeffsProfileOf(self, fp, ivar, indices=None,
                                  precision=6, sep=' ', linesep='\n', header='# '):
        headers = []
        exports = []
        fmtstr = '%% .%de' % precision
        sens_lists = self.getSensCoeffsProfileOf(ivar, indices)
        headers.append(self.xvarname)
        exports.append(self.getXvarList())
        
        idhdrs = list(range(self.nr)) if indices is None else indices[:]
        idhdrs = ['R%d' % (i+1,) for i in idhdrs]
        headers.extend(idhdrs)
        exports.extend(sens_lists)
        if header is not None:
            fp.write('%s%s%s' % (header, sep.join(headers),linesep))
        for x in zip(*exports):
            fp.write('%s%s' % (sep.join([fmtstr % (v,) for v in x]),linesep))

    # mass flux
    def getMassFractionAt(self, index):
        xvar, press, temp, frac, addv = self.getProfileAt(index)
        sumymw = sum([frac[j] * self.weights[j] for j in range(self.ns)])
        for j in range(self.ns):
            frac[j] = frac[j] * self.weights[j] / sumymw
        return frac
    
    def getNetMassFluxMatrixAt(self, index):
        """net mass flux matrix Akl at xvar[index] in g/cm3-s unit
        return A
        A[k][l]: net mass flux for species[k] -> species[l]
        """
        xvar, p, temp, x, addv = self.getProfileAt(index)
        molrates = self.ckqxp(p, temp, x)
        
        A = [[0.0 for k in range(self.ns)] for k in range(self.ns)]
        for i in range(self.nr):
            rdict, pdict = self.ckinu(i)
            wt = 0.
            for k in rdict: wt += rdict[k] * self.weights[k]
            for k in rdict:
                for l in pdict:
                    flux = molrates[i] * (self.weights[k] / wt) * self.weights[l]
                    A[k][l] += flux
                    A[l][k] -= flux
        #for k in range(self.ns):
        #    for l in range(self.ns): print self.species[k], self.species[l], A[k][l]
        return A
    
    def getNetMassFluxContrib(self, index, mat, connex):
        """net mass flux contributions
        return concontrib[(from,to)]: [(rind, frac)]
        """
        xvar, p, temp, x, addv = self.getProfileAt(index)
        molrates = self.ckqxp(p, temp, x)
        concontrib = {}
        for x in connex: concontrib[x] = []
        for i in range(self.nr):
            rdict, pdict = self.ckinu(i)
            wt = 0.
            for k in rdict: wt += rdict[k] * self.weights[k]
            for k in rdict:
                for l in pdict:
                    if (k,l) in connex:
                        flux = molrates[i] * (self.weights[k] / wt) * self.weights[l]
                        concontrib[(k,l)].append((i, flux / mat[k][l]))
                    elif (l,k) in connex:
                        flux = - molrates[i] * (self.weights[k] / wt) * self.weights[l]
                        concontrib[(l,k)].append((i, flux / mat[l][k]))
        for x in concontrib:
            concontrib[x].sort(key=lambda x: abs(x[1]), reverse=True)
        return concontrib

    # element flux
    
    def getElementFractionAt(self, index, element):
        ielem = self.elements.index(element)
        
        xvar, press, temp, frac, addv = self.getProfileAt(index)
        conv = 6.022137e23 * press / (temp * 8.314510e6) # molfrac -> molecule/cm3
        for k in range(self.ns): 
            frac[k] = frac[k] * conv * self.compositions[k][ielem]  # atoms/cm3
        sumatoms = sum(frac)
        if sumatoms > 1e-128:
            for k in range(self.ns): frac[k] = frac[k] / sumatoms
        return frac
        
    def getNetElementFluxMatrixAt(self, index, element):
        """net element flux matrix Akl at xvar[index] in atoms/cm3-s unit
        return A
        A[k][l]: net element flux for species[k] -> species[l]
        net element flux = rates[i] * (nk / N) * nl
        """
        xvar, p, temp, x, addv = self.getProfileAt(index)
        fluxes = [x*6.022137e23 for x in self.ckqxp(p, temp, x)]
        ielem = self.elements.index(element)
        A = [[0.0 for k in range(self.ns)] for k in range(self.ns)]
        for i in range(self.nr):
            rdict, pdict = self.ckinu(i)
            N = 0
            for k in rdict:
                N += self.compositions[k][ielem]
            if N == 0: continue
            for k in rdict:
                nk = self.compositions[k][ielem]
                if nk == 0: continue
                for l in pdict:
                    nl = self.compositions[l][ielem]
                    if nl == 0: continue
                    flux = fluxes[i] * (float(nk) / N) * nl
                    A[k][l] += flux
                    A[l][k] -= flux
        #for k in range(self.ns):
        #    for l in range(self.ns): print self.species[k], self.species[l], A[k][l]
        return A
    
    def getNetElementFluxContrib(self, index, element, mat, connex):
        """net element flux contributions
        return concontrib[(from,to)]: [(rind, frac)]
        """
        xvar, p, temp, x, addv = self.getProfileAt(index)
        fluxes = [x*6.022137e23 for x in self.ckqxp(p, temp, x)]
        ielem = self.elements.index(element)
        concontrib = {}
        for x in connex: concontrib[x] = []
        for i in range(self.nr):
            rdict, pdict = self.ckinu(i)
            N = 0
            for k in rdict:
                N += self.compositions[k][ielem]
            if N == 0: continue
            for k in rdict:
                nk = self.compositions[k][ielem]
                if nk == 0: continue
                for l in pdict:
                    nl = self.compositions[l][ielem]
                    if nl == 0: continue
                    if (k,l) in connex:
                        flux = fluxes[i] * (float(nk) / N) * nl
                        concontrib[(k,l)].append((i, flux / mat[k][l]))
                    elif (l,k) in connex:
                        flux = - fluxes[i] * (float(nk) / N) * nl
                        concontrib[(l,k)].append((i, flux / mat[l][k]))
        for x in concontrib:
            concontrib[x].sort(key=lambda x: abs(x[1]), reverse=True)
        return concontrib
        
    # flux analysis methods
    
    def filterFluxFromMatrix(self, mat, include=None, exclude=None, 
                             prior_included=False, included_only=False,
                             search_direction=0, remove_unconnected=True, rtol=1e-3):
        """filter flux matrix
          include: included species; None for all
          exclude: excluded species; None for none
          prior_included: included is piror or not
          included_only: no search
          search_direction: 0 for both, positive for consumption, negative for formation
          remove_unconnected: remove isolated species or not
          rtol: rel. tol. for flux (rel to maximum flux from included species)
        return spdict, connex
          spdict[id]: (sumfor, sumback, [(forind1, frac), ...], [(backind1, frac, ...)])
          connex[(from,to)]: (flux, consumption frac, formation frac)
        """
        if include is None: includes = list(range(self.ns))
        else: includes = [self.species.index(name) for name in include]
        
        if exclude is None: excludes = set()
        else: excludes = set([self.species.index(name) for name in exclude])
        if prior_included:
            for k in includes:
                if k in excludes: excludes.remove(k)
        else:
            for k in excludes:
                if k in includes: includes.remove(k)
        if included_only:
            excludes = set()
            for k in range(self.ns):
                if k not in includes: excludes.add(k)
        
        maxflux = 0.0
        for k in includes:
            for l in range(self.ns):
                val = mat[k][l]
                if (((l not in excludes) and (abs(val) > maxflux)) and
                    ((search_direction == 0) or
                     (search_direction > 0 and val > 0) or
                     (search_direction < 0 and val < 0))):
                    maxflux = abs(val)
        mintol = maxflux * rtol

        totfs = []
        for k in range(self.ns):
            sumfor, sumback = 0., 0.
            for l in range(self.ns):
                val = mat[k][l]
                if val > 0: sumfor += val
                elif val < 0: sumback += -val
            totfs.append((sumfor, sumback))

        que = includes[:]
        spinds = set()
        connex = {}
        i = 0
        while True:
            if i >= len(que): break
            k = que[i]
            i += 1
            if k in spinds: continue
            if k in excludes: continue
            spinds.add(k)
            for l in range(self.ns):
                val = mat[k][l]
                if (((l != k) and (l not in excludes) and (abs(val) > mintol)) and
                    ((search_direction == 0) or
                     (search_direction > 0 and val > 0) or
                     (search_direction < 0 and val < 0))):
                    if val > 0: con = (k,l)  # k -> l
                    else: con = (l,k)  # k <- l
                    if con not in connex:
                        # (flux, consumption frac, formation frac)
                        if val > 0:
                            flux, cfr, ffr = (val, val/totfs[k][0], val/totfs[l][1])
                        else:
                            flux, cfr, ffr = (-val, -val/totfs[l][0], -val/totfs[k][1])
                        connex[con] =(flux, cfr, ffr)
                        if l not in spinds and l not in que: que.append(l)
        if remove_unconnected:
            spinds2 = set()
            for k, l in connex:
                if k not in spinds2: spinds2.add(k)
                if l not in spinds2: spinds2.add(l)
            spinds = spinds2
        spdict = {}
        for ind in spinds:
            fl, bl = [], []
            for x in connex:
                k, l = x
                if k == ind: fl.append((l, connex[x][1]))
                elif l == ind: bl.append((k, connex[x][2]))
            # (sumfor, sumback, [(forind1, frac), ...], [(backind1, frac, ...)])
            fl.sort(key=lambda x: x[1], reverse=True)
            bl.sort(key=lambda x: x[1], reverse=True)
            spdict[ind] = (totfs[ind][0], totfs[ind][1], fl, bl)
        return spdict, connex
        

    def dumpFluxToTxt(self, fp, fracs, spdict, connex, concontrib, 
                      header="Flux", unit="g/cm3-s"):
        fp.write("%s\n" % (header,))
        fp.write("\n")
        fp.write("[Species]\n")
        for k in range(self.ns):
            if k not in spdict: continue
            fp.write('%s [id=%d] (%.3g %%):\n' \
                     % (self.species[k], k+1, fracs[k]*100))
            fp.write('  Consumption flux = %.4e %s\n' % (spdict[k][0], unit))
            for l in spdict[k][2]:
                fp.write('      -> %s [id=%d]  (%.3g %%)\n' %
                         (self.species[l[0]], l[0]+1, l[1]*100))
            fp.write('  Formation flux = %.4e %s\n' % (spdict[k][1], unit))
            for l in spdict[k][3]:
                fp.write('      <- %s [id=%d]  (%.3g %%)\n' %
                         (self.species[l[0]], l[0]+1, l[1]*100))
            fp.write('  Total flux = %.4e %s\n' % (spdict[k][1]-spdict[k][0], unit))
        fp.write("\n")
        fp.write("[Connections]\n")
        for k in range(self.ns):
            for l in range(self.ns):
                x = (k, l)
                if (k, l) not in connex: continue
                val, cfr, ffr = connex[x]
                contrib = concontrib[x]
                froms, tos = self.species[k], self.species[l]
                fp.write('%s [id=%d] -> %s [id=%d] (%.4e %s):\n' % \
                         (froms, k+1, tos, l+1, val, unit))
                fp.write('  Account for %.3g %% of %s consumption, ' % (cfr*100, froms))
                fp.write('%.3g %% of %s formation\n' % (ffr*100, tos))
                if len(contrib) > 0:
                    fp.write('  Contributions:\n')
                    for rind, rfrac in contrib:
                        fp.write('    %4d: %s  (%.3g %%)\n' % 
                                 (rind+1, self.reactions[rind], rfrac*100))
        return

    def dumpFluxToDot(self, fp, fracs, spdict, connex):
        fp.write("digraph g {\n")
        for k in spdict: fp.write('  "%s";\n' % (self.species[k],))
        
        maxflux = 0.0
        for x in connex: 
            if connex[x][0] > maxflux: maxflux = connex[x][0]
            
        for x in connex:
            flux, cfr, ffr = connex[x]
            k, l = x
            froms, tos = self.species[k], self.species[l]
            weight = max(1, int(100.*abs(flux/maxflux)))
            width = max(0.5, 10*abs(flux/maxflux)**0.5)
            fp.write('  "%s" -> "%s" [weight=%d, penwidth=%.2f];\n' %
                     (froms, tos, weight, width))
        fp.write("}\n")
                        



if __name__ == "__main__":
    ctout = CanteraOutput("test0d_1.dat")
    print(ctout.info)
    a = ctout.indexAt(1e-3)
    #ctout.exportProfiles(sys.stdout)
    print(ctout.getReactionRatesAt(a))
    frac = ctout.getMassFractionAt(a)
    mat = ctout.getNetMassFluxMatrixAt(a) 
    spdict, connex = ctout.filterFluxFromMatrix(mat, include=['H', "H2O"], rtol=1e-3)
    concontrib = ctout.getNetMassFluxContrib(a, mat, connex)
    ctout.dumpFluxToTxt(sys.stdout, frac, spdict, connex, concontrib)
    
