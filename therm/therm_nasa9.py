#! /usr/bin/env python3

import os
import numpy as np
import thermutils

_RU = 8.3144621 # J / (K mole)

class Nasa9 (object):

    def __init__(self, name, atoms, mw, thdata):
        self.name = name
        self.atoms = atoms
        self.mw = mw
        self.thdata = thdata

    def get_coeffs(self, T):
        for x in self.thdata:
            Tl, Th, coeffs = x
            if Tl <= T and T <= Th:
                return coeffs
        raise ValueError("not in T range: T = %s" % T)
    
    def cp_divR (self, T):
        a = self.get_coeffs(T)
        return sum([a[i]*T**(i-2) for i in range(7)])

    def h_divR (self, T):
        a = self.get_coeffs(T)
        return -a[0]*T**-1 + a[1]*np.log(T) + a[2]*T + a[3]/2.*T**2 + a[4]/3.*T**3 + \
               a[5]/4.*T**4 + a[6]/5.*T**5 + a[7]

    def s_divR (self, T):
        a = self.get_coeffs(T)
        return -a[0]/2.*T**-2 - a[1]*T**-1+ a[2]*np.log(T) + a[3]*T + a[4]/2.*T**2 + \
               a[5]/3.*T**3 + a[6]/4.*T**4 + a[8]

    def get_cp(self, T):
        return self.cp_divR(T) * _RU

    def get_h(self, T):
        return self.h_divR(T) * _RU / 1000.

    def get_s(self, T):
        return self.s_divR(T) * _RU

    def get_Trange(self):
        Tmin, Tmax = 1e32, 0.0
        for x in self.thdata:
            Tl, Th, coeffs = x
            if Tl < Tmin: Tmin = Tl
            if Th > Tmax: Tmax = Th
        return Tmin, Tmax

    def convert_to_nasa7(self, fitdT=10.0, atol=0.5):
        Tmin, Tmax = self.get_Trange()
        Tlst = []
        Cpl = []
        T = Tmin
        while T <= Tmax + fitdT/1e16:
            Tlst.append(T)
            Cpl.append(self.get_cp(T))
            T += fitdT
        H298, S298 = self.get_h(298.15), self.get_s(298.15)
        lt, mt, ht, coeffs, err = thermutils.fitNasa7(H298, S298, Tlst, Cpl)
        nasa7 = thermutils.Nasa7(self.name, " conv9 ", self.atoms,
                                 "G", lt, mt, ht, "", coeffs, err)
        nasa7.validate()
        nasa7.checkerrcp(Tlst, Cpl, atol=atol)
        return nasa7

def readNASA9dat(fn):
    data = {}
    fp = open(fn)
    for l in fp:
        if l.strip().upper().startswith("THERMO"): break
    next(fp)
    for l in fp:
        if l.strip().upper().startswith("END"): break
        if l.strip() == "": continue
        if l.strip().startswith("!") or l.strip().startswith("#"): continue
        if l.find("!") != -1: l = l[:l.index("!")]
        if l.find("#") != -1: l = l[:l.index("#")]
        l = l.rstrip()
        name = l[0:24].strip().split()[0]
        l = next(fp)
        numT = int(l[0:2].strip().split()[0])
        mw = float(l[52:65].strip().split()[0])
        atoms = []
        for i in range(5):
            aname = l[10+8*i:10+8*i+2]
            anum = l[10+8*i+2:10+8*i+8]
            if aname.strip() == "": atoms.append(None)
            elif anum.strip() == "": atoms.append(None)
            else:
                tmp = anum.split(".")
                if len(tmp) == 1: anum = int(tmp[0].strip())
                elif len(tmp) == 2 and tmp[1].strip().strip("0") == "": anum = int(tmp[0].strip())
                else: 
                    print("Warning: non-int atom num ignored: %s" % anum)
                    anum = 0
                atoms.append((aname, anum))
                
        thdata = []
        for i in range(numT):
            l = next(fp)
            Tr = l[1:21].strip().split()
            Tl, Th = float(Tr[0]), float(Tr[1])
            coeffs = []
            l = next(fp)
            coeffs.extend([l[0:16], l[16:32], l[32:48], l[48:64], l[64:80]])
            l = next(fp)
            coeffs.extend([l[0:16], l[16:32], l[48:64], l[64:80]])
            def parsecoeff(x):
                x = x.strip()
                if x == "": return 0.0
                x = x.upper().replace("D", "E")
                return float(x)
            coeffs = list(map(parsecoeff, coeffs))
            thdata.append([Tl, Th, coeffs])
        data[name] = Nasa9(name, atoms, mw, thdata)
    return data


def atom2str_nasa9(atom):
    if atom is None: return "%-2s%6.2f" % ("", float(0))
    else: return "%-2s%6.2f" % (atom[0], float(atom[1]))

def formcoeff_nasa9(x):
    s = "% 17.10E" % x
    if s[-3] == "0": s = s[:-3] + s[-2:]
    else: s = s[:-5] + s[-4:]
    s = s.replace("E", "D")
    return s

def nasa7_to_nasa9(nasa7, mw, h298_0):
    s = "%-24s%-56s\n" % (nasa7.name, nasa7.comment)
    s += " 2   conv %s%s%s%s%s" % (atom2str_nasa9(nasa7.atoms[0]), atom2str_nasa9(nasa7.atoms[1]),
                                 atom2str_nasa9(nasa7.atoms[2]), atom2str_nasa9(nasa7.atoms[3]),
                                 atom2str_nasa9(nasa7.atoms[4]))
    s += " 0%13.5f%15.3f\n" % (mw, nasa7.get_h(298.15)*1000.)
    s += " %10.3f %10.3f7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0  %15.3f\n" % (nasa7.lt, nasa7.mt, h298_0*1000.)
    tmp = nasa7.coeffs[7:]
    s += "%s%s%s%s%s\n" % (formcoeff_nasa9(0.0), formcoeff_nasa9(0.0), formcoeff_nasa9(tmp[0]),
                           formcoeff_nasa9(tmp[1]), formcoeff_nasa9(tmp[2]))
    s += "%s%s                %s%s\n" % (formcoeff_nasa9(tmp[3]), formcoeff_nasa9(tmp[4]),
                                         formcoeff_nasa9(tmp[5]), formcoeff_nasa9(tmp[6]))
    
    s += " %10.3f %10.3f7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0  %15.3f\n" % (nasa7.mt, nasa7.ht, h298_0*1000.)
    tmp = nasa7.coeffs[:7]
    s += "%s%s%s%s%s\n" % (formcoeff_nasa9(0.0), formcoeff_nasa9(0.0), formcoeff_nasa9(tmp[0]),
                           formcoeff_nasa9(tmp[1]), formcoeff_nasa9(tmp[2]))
    s += "%s%s                %s%s\n" % (formcoeff_nasa9(tmp[3]), formcoeff_nasa9(tmp[4]),
                                         formcoeff_nasa9(tmp[5]), formcoeff_nasa9(tmp[6]))
    return s

if __name__ == "__main__":
    
    data = thermutils.readdat("test-thermutils/thermch3cf3.dat")
    for nasa7 in data:
        print(nasa7_to_nasa9(nasa7, 0.00, 0.00))
    
    nasa9data = readNASA9dat("test-thermutils/nasa9_thermo.inp")
    name = "FO"
    #name = "HOF"
    
    nasa9 = nasa9data[name]
    nasa7 = nasa9.convert_to_nasa7(atol=1.5)
    print(nasa7)
    #print nasa7_to_nasa9(nasa7, 1., 0.)
    Tmin, Tmax = nasa9.get_Trange()
    T = Tmin
    dT = 400.0
    while T <= Tmax + dT/1e16:
        print(T, nasa9.get_cp(T), nasa7.get_cp(T), nasa9.get_h(T), nasa7.get_h(T), nasa9.get_s(T), nasa7.get_s(T))
        T += dT
    
