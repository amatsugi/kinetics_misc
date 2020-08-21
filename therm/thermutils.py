#! /usr/bin/env python3

import numpy as np
import sys
import os
import datetime 

import lsq


_RU = 8.3144621 # J / (K mole)
_cal2J = 4.184
_dHf_H = 218.0 # DHf of H  kJ/mol
_lstheader = " SPECIES    HF(298)  S(298)     CP300   CP400   CP500   CP600   CP800   CP1000  CP1500  DATE  REF           ELEMENTS "

def formcoeff(x):
    s = "% 16.9E" % x
    if s[-3] == "0": s = s[:-3] + s[-2:]
    else: s = s[:-5] + s[-4:]
    return s

def nasa7_cp(T, a):
    return (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4) * _RU

def nasa7_h(T, a):
    return (a[0]*T + a[1]/2.*T**2 + a[2]/3.*T**3 + a[3]/4.*T**4 + a[4]/5.*T**5 + a[5]) * _RU / 1000.

def nasa7_s(T, a):
    return (a[0]*np.log(T) + a[1]*T + a[2]/2.*T**2 + a[3]/3.*T**3 + a[4]/4.*T**4 + a[6]) * _RU

def nasa7_g(T, a):
    return nasa7_h(T, a) - nasa7_s(T, a) / 1000. * T
    
def readNasaChemfromstr(sl, def_lt=300., def_mt=1000., def_ht=5000.):
    l = sl[0]
    name = l[0:16].split()[0]
    comment = l[len(name):24].strip()
    atoms = [str2atom(l[24:29]), str2atom(l[29:34]),
             str2atom(l[34:39]), str2atom(l[39:44]),
             str2atom(l[73:78])]
    phase = l[44]
    lt = l[45:55]
    ht = l[55:65]
    mt = l[65:73]
    if lt.strip() == "": lt = def_lt
    else: lt = float(lt.strip())
    if ht.strip() == "": ht = def_ht
    else: ht = float(ht.strip())
    if mt.strip() == "": mt = def_mt
    else: mt = float(mt.strip())
    rotn = l[78]
    if rotn == " ": rotn = 0
    else: rotn = int(rotn)
    coeffs = []
    l = sl[1]
    coeffs.extend([l[0:15], l[15:30], l[30:45], l[45:60], l[60:75]])
    l = sl[2]
    coeffs.extend([l[0:15], l[15:30], l[30:45], l[45:60], l[60:75]])
    l = sl[3]
    coeffs.extend([l[0:15], l[15:30], l[30:45], l[45:60]])
    def parsecoeff(x):
        x = x.upper().replace("E ", "E")
        return float(x)
    coeffs = list(map(parsecoeff, coeffs))
    return Nasa7(name, comment, atoms, phase, lt, mt, ht, rotn, coeffs)
    

def atom2str(atom):
    if atom is None: return "     "
    else: return "%-2s%3d" % (atom[0], atom[1])

def str2atom(s):
    l = s.split()
    if len(l) == 0: return None
    elif len(l) != 2: return None

    if l[-1][-1] == ".": l[-1] = l[-1][:-1]
    if l[-1] == "0": return None
    elif l[0] == "0": return None
    else: return l[0], int(l[1])

def listen(name, h, s, cpl, comment, atoms, phase, rotn):
    s = " %-10s%8.2f %8.2f" % (name, h, s)
    s += " %8.2f%8.2f%8.2f%8.2f%8.2f%8.2f" \
         % (cpl[0], cpl[1], cpl[2], cpl[3], cpl[4], cpl[5])
    if len(cpl) == 6: s += "        "
    else: s += "%8.2f" % cpl[6]
    s += " %14s   %s %s %s %s" % (comment, atom2str(atoms[0]), atom2str(atoms[1]),
                                  atom2str(atoms[2]), atom2str(atoms[3]))
    s += " %s %s\n" % (phase, rotn)
    return s

def daten(name, comment, atoms, phase, lt, mt, ht, coeffs):
    if len(name) <= 10:
        namecomment = "%-10s%14s" % (name, comment)
    else:
        cfmt = "%%%ds" % (24 - len(name) - 1)
        namecomment = name + " " + (cfmt % comment)
    if len(namecomment) > 24:
        namecomment = namecomment[:24]
    s = "%s%s%s%s%s%s% 10.3f% 10.3f% 8.2f%s 1\n" \
        % (namecomment, atom2str(atoms[0]), atom2str(atoms[1]), atom2str(atoms[2]),
           atom2str(atoms[3]), phase, lt, ht, mt, atom2str(atoms[4]))
    s += "".join(map(formcoeff, coeffs[0:5])) + "    2\n"
    s += "".join(map(formcoeff, coeffs[5:10])) + "    3\n"
    s += "".join(map(formcoeff, coeffs[10:])) + "                   4\n"
    return s


class Lst (object):

    def __init__(self, name, h, s, cpl, comment, atoms, phase, rotn):
        self.name = name
        self.h = h
        self.s = s
        self.cpl = cpl
        self.comment = comment
        self.atoms = atoms
        self.phase = phase
        self.rotn = rotn

    def lst(self):
        return listen(self.name, self.h, self.s, self.cpl, self.comment, 
                      self.atoms, self.phase, self.rotn)

    # fit cp to extrap. (preliminary implementation: HOE and Wilhoit)
    
    def fitHOE(self, guess=[0., 1., 1000., 2000., 3000.]):
        natoms = 0
        for x in self.atoms:
            if x is not None:
                natoms += x[1]
        if natoms == 1: return self.fitatom()
        if natoms == 2: return self.fitdiatom()
        #nvib = 3 * natoms - 6
        nvib = 3 * natoms - 6 - self.rotn
        Tl = [300., 400., 500., 600., 800., 1000., 1500.]
        if len(self.cpl) == 6: Tl = Tl[:6]
        T = np.array(Tl)
        cp = np.array(self.cpl)
        
        def fvib(tdivT):
            return tdivT**2 * np.exp(-tdivT) * (np.exp(-tdivT) -1)**-2
        def cpHOE(p, x):
            b1, b2, t1, t2, t3 = p
            T = x
            Cvib = b1*fvib(t1/T) + b2*fvib(t2/T) + (nvib - (b1+b2)) * fvib(t3/T)
            #return (4. + Cvib) * _RU
            return (4. + self.rotn/2. + Cvib) * _RU
        
        res = lsq.NonLinear(cpHOE, T, cp, guess, maxfev=1000000)
        if not res.issolved(): return None

        def cpfunc(T): return cpHOE(res.params, T)
        return cpfunc
            
    def fitWilhoit(self, npar=3):
        natoms = 0
        for x in self.atoms:
            if x is not None:
                natoms += x[1]
        #if natoms == 1: return self.fitatom()
        #if natoms == 2: return self.fitdiatom()
        Tl = [300., 400., 500., 600., 800., 1000., 1500.]
        if len(self.cpl) == 6: Tl = Tl[:6]
        T = np.array(Tl)
        cp = np.array(self.cpl)
        #c0 = 4*_RU
        #cinf = (3*natoms-2)*_RU
        c0 = (4 + self.rotn/2.)*_RU
        cinf = (3*natoms-2-self.rotn/2.)*_RU
        def cpWil(T, S, al):
            y = T / (T+S)
            sm = 0.
            for i in range(len(al)):
                sm += al[i] * y**i
            return c0 + (cinf-c0)*y*y*(1. + (y-1.)*sm)
        
        minlsqerr = None
        for S in range(50, 5001, 10):
            S = float(S)
            y = T / (T+S)
            obj = ((cp-c0) / ((cinf-c0)*y*y) - 1.) / (y-1.)
            res = lsq.Poly(y, obj, npar)
            al = res.params
            fitcp = cpWil(T, S, al)
            lsqerr = sum(np.square(cp - fitcp))
            if minlsqerr is None or lsqerr < minlsqerr:
                minlsqerr = lsqerr
                Smin = S
                almin = al
        def cpfunc(T): return cpWil(T, Smin, almin)
        return cpfunc

    def fitatom(self):
        if abs(self.cpl[0] - self.cpl[-1]) <= 0.01:
            # ignore temperature dependence
            def cpfunc(T): return self.cpl[0]
            return cpfunc
        Tl = [300., 400., 500., 600., 800., 1000., 1500.]
        if len(self.cpl) == 6: Tl = Tl[:6]
        T = np.array(Tl)
        cp = np.array(self.cpl)
        # felec is preliminary 
        def felec(a, edivT1, edivT2):
            tmp = 1. + a*np.exp(-edivT1) + (1.-a)*np.exp(-edivT2)
            return ((a*edivT1**2*np.exp(-edivT1) + (1.-a)*edivT2**2*np.exp(-edivT2)) / tmp
                    - ((a*edivT1*np.exp(-edivT1) + (1.-a)*edivT2*np.exp(-edivT2)) / tmp)**2)
        def cpHOE(p, x):
            a, te1, te2 = p
            T = x
            return (2.5 + felec(a, te1/T, te2/T)) * _RU
        res = lsq.NonLinear(cpHOE, T, cp, guess=[0.8, 100., 200.], maxfev=1000000)
        if not res.issolved(): return None
        def cpfunc(T): return cpHOE(res.params, T)
        return cpfunc

    def fitdiatom(self):
        Tl = [300., 400., 500., 600., 800., 1000., 1500.]
        if len(self.cpl) == 6: Tl = Tl[:6]
        T = np.array(Tl)
        cp = np.array(self.cpl)

        # electronic contrib. is missing for diatomic molecules
        def fvib(tdivT):
            return tdivT**2 * np.exp(-tdivT) * (np.exp(-tdivT) -1)**-2
        def cpHOE(p, x):
            t1 = p[0]
            T = x
            return (3.5 + fvib(t1/T)) * _RU
        res = lsq.NonLinear(cpHOE, T, cp, guess=[1000.], maxfev=1000000)
        if not res.issolved(): return None
        def cpfunc(T): return cpHOE(res.params, T)
        return cpfunc
        
        

    def fiterr(self, cpfunc):
        Tl = [300., 400., 500., 600., 800., 1000., 1500.]
        if len(self.cpl) == 6: Tl = Tl[:6]
        T = np.array(Tl)
        cp = np.array(self.cpl)
        fitcp = cpfunc(T)
        rmsd = np.sqrt(sum(np.square(cp - fitcp)) / float(len(cp)))
        mad = sum(abs(cp - fitcp)) / float(len(cp))
        maxad = max(abs(cp - fitcp))
        return rmsd, mad, maxad
        

class Nasa7 (object):

    def __init__(self, name, comment, atoms, phase, lt, mt, ht, rotn, coeffs, err=None):
        self.name = name
        self.comment = comment
        self.atoms = atoms
        self.phase = phase
        self.lt = lt
        self.mt = mt
        self.ht = ht
        self.rotn = rotn
        self.coeffs = coeffs
        # err: (rmsd, mad, maxad)
        if err is None: self.err = (0., 0., 0.)
        else: self.err = err

    def __str__(self):
        return self.dat()
        
    def dat(self):
        return daten(self.name, self.comment, self.atoms, self.phase, 
                     self.lt, self.mt, self.ht, self.coeffs)

    def lst(self):
        cpl = [self.get_cp(300.), self.get_cp(400.), self.get_cp(500.),
               self.get_cp(600.), self.get_cp(800.), self.get_cp(1000.),
               self.get_cp(1500.)]
        return listen(self.name, self.get_h(298.15), self.get_s(298.15),
                      cpl, self.comment, self.atoms, self.phase, self.rotn)
        return s

    def get_coeff7(self, T):
        if T > self.mt: return self.coeffs[:7]
        else: return self.coeffs[7:]

    def get_cp(self, T):
        return nasa7_cp(T, self.get_coeff7(T))

    def get_h(self, T):
        return nasa7_h(T, self.get_coeff7(T))

    def get_s(self, T):
        return nasa7_s(T, self.get_coeff7(T))

    def get_g(self, T):
        return nasa7_g(T, self.get_coeff7(T))

    def validate(self, rtol=1.e-3):
        hcoeffs = self.coeffs[:7]
        lcoeffs = self.coeffs[7:]
        # s and cp must be positive
        chkpoints = 100
        for i in range(chkpoints+1):
            T = (self.ht - self.lt) * i / float(chkpoints) + self.lt
            s, cp = self.get_s(T), self.get_cp(T)
            if s < 0.: raise ValueError("%s: S = %s at T = %s K" % (self.name, s, T))
            if cp < 0.: raise ValueError("%s: Cp = %s at T = %s K" % (self.name, cp, T))
        # check continuity at mt
        lh, ls, lcp = nasa7_h(self.mt, lcoeffs), nasa7_s(self.mt, lcoeffs), nasa7_cp(self.mt, lcoeffs)
        hh, hs, hcp = nasa7_h(self.mt, hcoeffs), nasa7_s(self.mt, hcoeffs), nasa7_cp(self.mt, hcoeffs)
        if abs((hh-lh)/lh) > rtol:
            raise ValueError("%s: H(low,high) = (%s,%s) at T = %s K" % (self.name, hh, lh, self.mt))
        if abs((hs-ls)/ls) > rtol: 
            raise ValueError("%s: S(low,high) = (%s,%s) at T = %s K" % (self.name, hs, ls, self.mt))
        if abs((hcp-lcp)/lcp) > rtol:
            raise ValueError("%s: Cp(low,high) = (%s,%s) at T = %s K" % (self.name, hcp, lcp, self.mt))
        return True
    
    def checkerr(self, Tl, Sl, Cpl, Hl, atol=(0.5,0.5,0.5)):
        for i in range(len(Tl)):
            T, S, Cp, H = Tl[i], Sl[i], Cpl[i], Hl[i]
            nS, nCp, nH = self.get_s(T), self.get_cp(T), self.get_h(T)
            if abs(nH-H) > atol[0]:
                raise ValueError("%s: H(fitted,correct) = (%s,%s) at T = %s K" % (self.name, nH, H, T))
            if abs(nS-S) > atol[1]:
                raise ValueError("%s: S(fitted,correct) = (%s,%s) at T = %s K" % (self.name, nS, S, T))
            if abs(nCp-Cp) > atol[2]:
                raise ValueError("%s: Cp(fitted,correct) = (%s,%s) at T = %s K" % (self.name, nCp, Cp, T))
        return True
    
    def checkerrcp(self, Tl, Cpl, atol=0.5):
        for i in range(len(Tl)):
            T, Cp = Tl[i], Cpl[i]
            nCp = self.get_cp(T)
            if abs(nCp-Cp) > atol:
                raise ValueError("%s: Cp(fitted,correct) = (%s,%s) at T = %s K" % (self.name, nCp, Cp, T))
        return True


def readdat(fn, def_lt=300., def_mt=1000., def_ht=5000.):
    data = []
    fp = open(fn)
    for l in fp:
        if l.strip().upper().startswith("THERMO"): break
    for l in fp:
        if l.strip().upper().startswith("END"): break
        if l.strip() == "": continue
        if l.strip().startswith("!") or l.strip().startswith("#"): continue
        if l.find("!") != -1: l = l[:l.index("!")]
        ##if l.find("#") != -1: l = l[:l.index("#")]
        l = l.rstrip()
        if (len(l) < 80):
            try: def_lt, def_mt, def_ht = list(map(float, l.strip().split()))
            except: pass
        if (len(l) >= 80 and l[79] == "1"):
            data.append(readNasaChemfromstr([l,next(fp),next(fp),next(fp)],
                                        def_lt=def_lt, def_mt=def_mt, def_ht=def_ht))
    return data


def searchdat(fn, names, def_lt=300., def_mt=1000., def_ht=5000.):
    datadict= {}
    fp = open(fn)
    for l in fp:
        if l.strip().upper().startswith("THERMO"): break
    for l in fp:
        if l.strip().upper().startswith("END"): break
        if l.strip() == "": continue
        if l.strip().startswith("!") or l.strip().startswith("#"): continue
        if l.find("!") != -1: l = l[:l.index("!")]
        ##if l.find("#") != -1: l = l[:l.index("#")]
        l = l.rstrip()
        if (len(l) < 80):
            try: def_lt, def_mt, def_ht = list(map(float, l.strip().split()))
            except: pass
        if (len(l) >= 80 and l[79] == "1" and l[0:16].split()[0] in names):
            nasa7 = readNasaChemfromstr([l,next(fp),next(fp),next(fp)],
                                   def_lt=def_lt, def_mt=def_mt, def_ht=def_ht)
            datadict[nasa7.name] = nasa7
    return datadict

def readlst(fn):
    lsts = []
    fp = open(fn)
    l = next(fp)
    if l.strip().upper().startswith("UNITS:KCAL"): usecal = True
    if l.strip().upper().startswith("UNITS:KJ"): usecal = False
    title = next(fp).strip()
    l = next(fp)
    for l in fp:
        if l.strip() == "": continue
        if l.strip().startswith("!") or l.strip().startswith("#"): continue
        if l.find("!") != -1: l = l[:l.index("!")]
        if l.find("#") != -1: l = l[:l.index("#")]
        if l[:11].strip() == "": continue
        name = l[1:].split()[0]
        if len(name) >= 10: pnt = len(name) + 1
        else: pnt = 11
        values = list(map(float, l[pnt:pnt+74].split()))
        if usecal: values = list([x*_cal2J for x in values])
        h, s, cpl = values[0], values[1], values[2:]
        comment = l[pnt+75:pnt+89]
        atoms = [str2atom(l[pnt+92:pnt+97]), str2atom(l[pnt+98:pnt+103]),
                 str2atom(l[pnt+114:pnt+109]), str2atom(l[pnt+110:pnt+115]),
                 None]
        phase = l[pnt+116]
        rotn = int(l[pnt+118])
        lsts.append(Lst(name, h, s, cpl, comment, atoms, phase, rotn))
    return lsts
    
def writedat(data_list, fn=None):
    if fn is None: fp = sys.stdout
    else: fp = open(fn, "wb")
    fp.write("THERMO\n")
    fp.write("   300.000  1500.000  5000.000\n")
    for species in data_list:
        fp.write(species.dat())
    fp.write("END\n")
    if fn is not None: fp.close()
    

def writelst(data_list, fn=None, title=None):
    if fn is None: fp = sys.stdout
    else: fp = open(fn, "wb")
    fp.write(" UNITS:KJ\n ")
    if title is not None: fp.write(title)
    fp.write("\n")
    fp.write(_lstheader)
    fp.write("\n")
    for species in data_list:
        fp.write(species.lst())
    if fn is not None: fp.close()


# fit nasa7
# TODO: constraint to derivatives (DCp/Dt)

def fitNasa7(h, s, T, cp):
    T = np.array(T)
    cp = np.array(cp)
    lt = T[0]
    ht = T[-1]
    minlsqerr = None
    for i in range(len(T)-10):
        iT = i+5
        mt = T[iT]
        cpmid = cp[iT]
        clow = fitNasa7_1(T[:iT], cp[:iT], mt, cpmid)
        chigh = fitNasa7_1(T[iT:], cp[iT:], mt, cpmid)
        fitcplow = nasa7_cp(T[:iT], clow)
        fitcphigh = nasa7_cp(T[iT:], chigh)
        lsqerr = sum(np.square(cp[:iT] - fitcplow)) + sum(np.square(cp[iT:] - fitcphigh))
        if minlsqerr is None or lsqerr < minlsqerr:
            minlsqerr = lsqerr
            iTmin = iT
            clowmin = clow
            chighmin = chigh
    clow = clowmin
    chigh = chighmin
    iT = iTmin
    mt = T[iT]
    fitcplow = nasa7_cp(T[:iT], clow)
    fitcphigh = nasa7_cp(T[iT:], chigh)
    rmsd = np.sqrt((sum(np.square(cp[:iT] - fitcplow)) + sum(np.square(cp[iT:] - fitcphigh)))
                   / float(len(cp)))
    mad = (sum(abs(cp[:iT] - fitcplow)) + sum(abs(cp[iT:] - fitcphigh))) / float(len(cp))
    maxad = max(max(abs(cp[:iT] - fitcplow)), max(abs(cp[iT:] - fitcphigh)))
    tmp = 298.15
    a6low = h*1000./_RU - (clow[0]*tmp + clow[1]/2.*tmp**2 + clow[2]/3.*tmp**3 + clow[3]/4.*tmp**4 + clow[4]/5.*tmp**5)
    a7low = s/_RU - (clow[0]*np.log(tmp) + clow[1]*tmp + clow[2]/2.*tmp**2 + clow[3]/3.*tmp**3 + clow[4]/4.*tmp**4)
    lowcoeffs = [clow[0], clow[1], clow[2], clow[3], clow[4], a6low, a7low]
    tmp = mt
    mh = nasa7_h(mt, lowcoeffs)
    ms = nasa7_s(mt, lowcoeffs)
    a6high = mh*1000./_RU - (chigh[0]*tmp + chigh[1]/2.*tmp**2 + chigh[2]/3.*tmp**3 + chigh[3]/4.*tmp**4 + chigh[4]/5.*tmp**5)
    a7high = ms/_RU - (chigh[0]*np.log(tmp) + chigh[1]*tmp + chigh[2]/2.*tmp**2 + chigh[3]/3.*tmp**3 + chigh[4]/4.*tmp**4)
    coeffs = [chigh[0], chigh[1], chigh[2], chigh[3], chigh[4], a6high, a7high,
              clow[0], clow[1], clow[2], clow[3], clow[4], a6low, a7low]
    return lt, mt, ht, coeffs, (rmsd, mad, maxad)

def fitNasa7_1(T, cp, mt, cpmid):
    T = np.array(T)
    cp = np.array(cp)
    obj = (cp - cpmid) / _RU
    glist = [lambda x: x - mt,
             lambda x: x**2 - mt**2,
             lambda x: x**3 - mt**3,
             lambda x: x**4 - mt**4]
    res = lsq.Linear(T, obj, glist)
    a2, a3, a4, a5 = res.params
    a1 = cpmid / _RU - (a2*mt + a3*mt**2 + a4*mt**3 + a5*mt**4)
    return a1, a2, a3, a4, a5

def readGPOPcsv(fn):
    fp = open(fn)
    next(fp)
    T, S298, Cp298, H298m0 = list(map(float, next(fp).strip().split(",")))
    next(fp)
    Tl, Sl, Cpl, DHl = [], [], [], []
    for l in fp:
        T, S, Cp, DH = list(map(float, l.strip().split(",")))
        Tl.append(T)
        Sl.append(S)
        Cpl.append(Cp)
        DHl.append(DH)
    return S298, Cp298, H298m0, Tl, Sl, Cpl, DHl

def fitNasa7fromfile(fn, atol=(0.5,0.5,0.5)):
    data = []
    fp = open(fn)
    basedir = os.path.dirname(fn)
    for l in fp:
        if l.strip() == "": continue
        if l.strip().startswith("!") or l.strip().startswith("#"): continue
        if l.find("!") != -1: l = l[:l.index("!")]
        if l.find("#") != -1: l = l[:l.index("#")]
        s1, s2, s3 = l.split("$")
        l1 = s1.split()
        name, csvfn, H298, phase = l1[0], l1[1], float(l1[2]), l1[3]
        if len(l1) > 4: rotn = int(l1[4])
        else: rotn = 0
        atoms = []
        l2 = s2.split()
        n = 0
        while len(l2) > n:
            atoms.append((l2[n], int(l2[n+1])))
            n += 2
        if len(atoms) < 5:
            for i in range(5-len(atoms)):
                atoms.append(None)
        comment = s3.strip()
        if basedir == "": pass
        else: csvfn = basedir + "/" + csvfn
        S298, Cp298, H298m0, Tl, Sl, Cpl, DHl = readGPOPcsv(csvfn)
        Hl = [x + H298 for x in DHl]
        lt, mt, ht, coeffs, err = fitNasa7(H298, S298, Tl, Cpl)
        nasa7 = Nasa7(name, comment, atoms, phase, lt, mt, ht, rotn, coeffs, err)
        nasa7.validate()
        nasa7.checkerr(Tl, Sl, Cpl, Hl, atol)
        data.append(nasa7)
    return data

def fitNasa7fromlst(lsts, lt=300.0, ht=5000.0, atol=0.5):
    data = []
    for lst in lsts:
        func = lst.fitHOE()
        if func is None:
            print("!!!!")
            func = lst.fitWilhoit()
        rmsd, mad, maxad = lst.fiterr(func)
        Tl, Cpl = [], []
        T = lt
        while T <= ht:
            Tl.append(T)
            Cpl.append(func(T))
            T += 10.
        lt, mt, ht, coeffs, err = fitNasa7(lst.h, lst.s, Tl, Cpl)
        nasa7 = Nasa7(lst.name, lst.comment, lst.atoms, lst.phase, lt, mt, ht, 
                    lst.rotn, coeffs, err)
        nasa7.validate()
        #print lst.name, mt, rmsd, mad, maxad, err
        nasa7.checkerrcp(Tl, Cpl, atol)
        data.append(nasa7)
    return data

## gen

def readgencfg(cfgfn):
    groupvalues = {}
    cfgdir = os.path.dirname(cfgfn)
    fp = open(cfgfn)
    l = next(fp)
    n = int(l.strip().split()[1])
    for i in range(n):
        fn = next(fp).strip()
        fp1 = open(cfgdir + "/" + fn)
        next(fp1)
        tmp = next(fp1).strip().replace(",", " ")
        try:
            ngroup = int(tmp.split()[0])
        except:
            tmp = next(fp1).strip().replace(",", " ")
            ngroup = int(tmp.split()[0])
        for j in range(ngroup):
            tmp = next(fp1).strip().replace(",", " ").replace("- ", "-").split()
            if len(tmp) == 0: break
            name = tmp[0]
            h = float(tmp[1])*_cal2J
            s = float(tmp[2])*_cal2J
            try: cpl = list([float(x)*_cal2J for x in tmp[3:10]])
            except: cpl = list([float(x)*_cal2J for x in tmp[3:9]])
            groupvalues[name] = [h, s, cpl]
        fp1.close()
    fp.close()
    return groupvalues

def readgen(infn, groupvalues):
    lsts = []
    fp = open(infn)
    title = next(fp)
    comment = datetime.datetime.now().strftime("gen%Y%m")
    for l in fp:
        if l.strip() == "": continue
        if l.strip().startswith("!") or l.strip().startswith("#"): continue
        if l.find("!") != -1: l = l[:l.index("!")]
        if l.find("#") != -1: l = l[:l.index("#")]
        
        s1, s2, s3 = l.split("$")
        l1 = s1.split()
        name, typ = l1[0], l1[1]
        typ = typ.upper()
        if typ not in ["R", "M"]:
            raise ValueError("error in typ: %s" % typ)
        atoms = []
        n = 2
        while len(l1) > n:
            atoms.append((l1[n], int(l1[n+1])))
            n += 2
        if len(atoms) < 5:
            for i in range(5-len(atoms)):
                atoms.append(None)
        l2 = s2.split()
        groups = []
        n = 0
        while len(l2) > n:
            groups.append([l2[n], int(l2[n+1])])
            n += 2
        l3 = s3.split()
        if typ == "M": 
            rotn, rotsym = int(l3[0]), int(l3[1])
        else: 
            radgroup = l3[0]
            rotn, rotsym = int(l3[1]), int(l3[2])
            
        h, s = 0., 0.
        cpl = [0. for x in range(7)]
        no1500 = False
        for group, num in groups:
            values = groupvalues[group]
            h += values[0] * num
            s += values[1] * num
            for i in range(len(values[2])):
                cpl[i] += values[2][i] * num
            if len(values[2]) < 7: no1500 = True
        if typ == "R":
            values = groupvalues[radgroup]
            h += values[0] - _dHf_H
            s += values[1]
            for i in range(len(values[2])):
                cpl[i] += values[2][i]
            if len(values[2]) < 7: no1500 = True
        s -= _RU * np.log(float(rotsym))
        if no1500:
            cpl = cpl[:6]
        lsts.append(Lst(name, h, s, cpl, comment, atoms, "G", rotn))
    return lsts


if __name__ == "__main__":
    if 0:
        data = readdat("test-thermutils/therm.dat")
        print(data[0].get_h(298.15))
        print(data[0].get_s(298.15))
        print(data[0].get_cp(300.))
        print(data[0].lst())
        writelst(data)
        writedat(data)
    if 0:
        datadict = searchdat("test-thermutils/therm.dat", ["C3H3"])
        print(datadict["C3H3"])
    if 0:
        groupvalues = readgencfg("test-thermutils/THERM.CFG")
        lsts = readgen("test-thermutils/gen.inp", groupvalues)
        writelst(lsts)
    if 0:
        #lsts = readlst("test-thermutils/a.lst")
        lsts = readlst("test-thermutils/test.lst")
        data = readdat("test-thermutils/test.dat")
        dic = dict((x.name, x) for x in data)
        #writelst(data)
        for lst in lsts:
            func = lst.fitHOE()
            if func is None:
                print("!!!!")
                func = lst.fitWilhoit()
            #print lst.name, lst.fiterr(func1)#, lst.fiterr(func2)
            print("%-10s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" \
                  % (lst.name, func(300.), func(500.), func(1000.), func(1500.), 
                     func(2000.), func(3000.), func(4000.), func(5000.)))
            func_bk = func
            func = dic[lst.name].get_cp
            print("%-10s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" \
                  % (lst.name, func(300.), func(500.), func(1000.), func(1500.), 
                     func(2000.), func(3000.), func(4000.), func(5000.)))
            print("  ", lst.fiterr(func_bk))
            #print lst.name, func(300.), func(400.), func(500.), func(600.), func(800.), func(1000.), func(1500.),
            #print func(2000.), func(2500.), func(3500.), func(4000.), func(5000.)
    if 0:
        lsts = readlst("test-thermutils/test.lst")
        data = readdat("test-thermutils/test.dat")
        dic = dict((x.name, x) for x in data)
        datafit = fitNasa7fromlst(lsts)
        for x in datafit:
            func = x.get_cp
            print("%-10s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" \
                  % (x.name, func(300.), func(500.), func(1000.), func(1500.), 
                     func(2000.), func(3000.), func(4000.), func(5000.)))
            func = dic[x.name].get_cp
            print("%-10s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f" \
                  % (x.name, func(300.), func(500.), func(1000.), func(1500.), 
                     func(2000.), func(3000.), func(4000.), func(5000.)))
            
    if 1:
        data = fitNasa7fromfile("test-thermutils/fit.inp")
        nasa7 = data[0]
        
        datadict = searchdat("test-thermutils/therm.dat", ["C3H3"])
        ref = datadict["C3H3"]
        
        print(nasa7)
        print(ref)
        print(nasa7.lst())
        print(ref.lst())
        print(nasa7.err)
        for T in range(200, 5001, 10):
            print(T, nasa7.get_h(T)-nasa7.get_h(298.15), ref.get_h(T)-ref.get_h(298.15), nasa7.get_s(T), ref.get_s(T), nasa7.get_cp(T), ref.get_cp(T))
