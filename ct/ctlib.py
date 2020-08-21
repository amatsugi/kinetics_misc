#! /usr/bin/env python3

import sys
import os
import shutil
import time
import array

import cantera as ct
from cantera import FreeFlame
from cantera import BurnerFlame
from cantera import CounterflowDiffusionFlame
#import cantera.XML
from cantera import ck2cti
from cantera import ctml_writer



def filter_xvars_from_dat(datfn, atol=1e-12, rtol=0.01,
                          dxmin=1e-8, dxmax=1e1, xmax=None, alldata=False):
    xvars = []
    fp = open(datfn)
    next(fp)
    next(fp)
    flag = False
    vals = None
    vals_last = None
    for l in fp:
        vals = list(map(float, l.split()))
        if xmax is not None and xmax < vals[0]:
            flag = False
            break
        if vals_last is None:
            vals_last = [0. for x in vals]
            flag = True
        elif alldata: flag = True
        elif (vals[0] - vals_last[0]) < dxmin: flag = False
        elif (vals[0] - vals_last[0]) > dxmax: flag = True
        elif any((vals[i] > atol and abs((vals[i]-vals_last[i]) / vals[i]) > rtol)
                 for i in range(1,len(vals))):
            flag = True
        else: flag = False
        if flag:
            xvars.append(vals[0])
            vals_last = vals[:]
    if not flag and vals is not None:
        xvars.append(vals[0])
    return xvars



def getfrac(n2, phi, oratio, name_n2, name_fuel, name_oxid):
    fuel = (1. - n2) * phi / (phi + oratio)
    oxid = (1. - n2) * oratio / (phi + oratio)
    return {name_n2: n2, name_fuel: fuel, name_oxid: oxid}



class SimBase (object):

    def __init__(self, chemfn="chem.inp", thermfn="thermo.dat", tranfn=None,
                 idtag="", logfp=None):
        if logfp is None: self.logfp = sys.stdout
        else: self.logfp = logfp
        
        if idtag is None or idtag.strip() == "": idtag = os.path.splitext(chemfn)[0]
        idtag = idtag.strip().replace(" ", "_").replace("\t", "_").replace("\n", "_")
        self.idtag = idtag

        self.log(" Date: %s\n\n" % time.strftime('%Y-%m-%d %H:%M:%S'))
        self.log(" Mechanism: %s\n" % chemfn)
        if os.path.exists(thermfn):
            self.log(" Thermo: %s\n" % thermfn)
        else:
            self.log(" No thermo database\n")
            thermfn = None
        
        self.log("\n Converting Chemkin-format input file to CTI & CTML format.\n")
        
        self.name = os.path.splitext(chemfn)[0]
        ctifn = self.name + ".cti"
        ctmlfn = self.name + ".xml"

        options = ["--input=%s" % chemfn]
        if thermfn is not None: options.append("--thermo=%s" % thermfn)
        if tranfn is not None: options.append("--transport=%s" % tranfn)
        if self.idtag is not None: options.append("--id=%s" % self.idtag)
        options.append("--output=%s" % ctifn) 
        options.append("--permissive")
        ck2cti.main(options)
        
        ctml_writer.convert(ctifn)
        self.log(" ... Done.\n\n")
        self.log(" CTI: %s\n" % ctifn)
        self.log(" CTML: %s\n\n" % ctmlfn)
        
        self.gas = ct.Solution(ctmlfn)
        self.srcfn = ctmlfn

    def log(self, s):
        self.logfp.write(s)
        self.logfp.flush()


class ZeroD (SimBase):

    def __init__(self, chemfn="chem.inp", thermfn="thermo.dat", idtag="", logfp=None):
        SimBase.__init__(self, chemfn=chemfn, thermfn=thermfn, idtag=idtag, logfp=logfp)

    def set_gas(self, T, P, reac):
        if isinstance(reac, dict):
            X = ", ".join(["%s:%g" % (nm, reac[nm]) for nm in reac])
        else:
            X = reac
        self.log(" T = %g K\n" % T)
        self.log(" P = %g atm\n" % P)
        self.log(" Reactants = %s\n\n" % X)
        self.gas.TPX = T, P*ct.one_atm, X
        return

    def _initsim(self, outname, conP=False, energy="on", rtol=1e-9, atol=1e-15, nodat=False):
        
        # files
        if outname is None:
            count = 0
            while True:
                count += 1
                outname = self.name + ("_zd%03d" % count)
                if not os.path.exists(outname+".out"): break
        outfn = outname + ".out"
        datfn = outname + ".dat"
        origctmlfn = self.srcfn
        ctmlfn = outname + ".xml"
        shutil.copy(origctmlfn, ctmlfn)
        
        self.log(" Output: %s\n" % outfn)
        if nodat: self.log(" No full data output\n")
        else: self.log(" Full data output: %s\n" % datfn)
        self.log(" CTML file was copied from %s to %s\n\n" % (origctmlfn, ctmlfn))
        
        ofp = open(outfn, "w")
        if not nodat: dfp = open(datfn, "w")
        def outlog(s):
            self.log(s)
            ofp.write(s)
            ofp.flush()
        def dataout(l):
            if not nodat:
                for x in l: dfp.write("%s " % x)
                dfp.write("\n")
                dfp.flush()

        # reactor
        #r = ct.Reactor(self.gas)
        if conP:
            r = ct.IdealGasConstPressureReactor(self.gas, energy=energy)
        else:
            r = ct.IdealGasReactor(self.gas, energy=energy)
        sim = ct.ReactorNet([r])
        sim.rtol = rtol
        sim.atol = atol
        sim.set_initial_time(0.0)

        # msg
        outlog(" Date: %s\n\n" % time.strftime('%Y-%m-%d %H:%M:%S'))
        outlog(" Initial state:\n")
        outlog("  " + str(self.gas()) + "\n\n")
        if conP: outlog(" Pressure is held constant.\n")
        else: outlog(" Volume is held constant.\n")
        outlog(" Relative error tolerance = % 14.6E\n" % rtol)
        outlog(" Absolute error tolerance = % 14.6E\n" % atol)
        outlog("\n")
        
        return outname, outlog, dataout, r, sim
        

    def run(self, tmax, conP=False, energy="on", sens=False,
            rtol=1e-9, atol=1e-15, rtols=1e-6, atols=1e-9, ignDT=400., outname=None):
        outname, outlog, dataout, r, sim = \
                self._initsim(outname=outname, conP=conP, energy=energy, rtol=rtol, atol=atol)
        
        # write headers
        idtag = self.idtag
        if idtag.strip() == "": idtag = "(noid)"
        species = list(self.gas.species_names)
        
        if sens:
            if conP: dataout(["ZeroD", outname, idtag, "Sens_hp"])
            else: dataout(["ZeroD", outname, idtag, "Sens"])
        else:
            dataout(["ZeroD", outname, idtag])
        dataout(["Time(s)", "T(K)", "P(Pa)", "V(m3)", "rho(kg/m3)"] + species)
        
        if sens:
            sensfn = outname + ".sens"
            for i in range(self.gas.n_reactions):
                r.add_sensitivity_reaction(i)
            sim.atol_sensitivity = atols
            sim.rtol_sensitivity = rtols
            sensfp = open(sensfn, "wb")
            outlog(" Sensitivity analysis\n")
            outlog(" Relative sens tolerance = % 14.6E\n" % rtols)
            outlog(" Absolute sens tolerance = % 14.6E\n\n" % atols)
        outlog(" Time Integration:\n\n")
        outlog("  % 14s % 14s % 14s % 14s % 14s\n"
                 % ("Time(s)", "T(K)", "P(Pa)", "V(m3)", "rho(kg/m3)"))
        
        clock0 = time.clock()
        time0 = time.time()
        
        # t = 0
        if not sens:
            data1 = [sim.time, r.T, r.thermo.P, r.volume, r.density]
            data2 = data1 + list(self.gas.X)
            dataout(["% 18.10E" % x for x in data2])
            outlog("  % 14.6E % 14.6E % 14.6E % 14.6E % 14.6E\n" % tuple(data1))

        # integrate
        ignT = r.T + ignDT
        igntime = None
        ndata = 1
        delt = tmax / 100.
        tnow = 0.0
        tnext = 0.0
        told = 0.0
        tmpold = r.T
        while tnow < tmax:
            tnext += delt
            while tnow < tnext:
                tnow = sim.step()
                ndata += 1
                if igntime is None and r.T > ignT:
                    igntime = told + (tnow - told) * (ignT - tmpold) / (r.T - tmpold)
                data1 = [sim.time, r.T, r.thermo.P, r.volume, r.density]
                data2 = data1 + list(self.gas.X)
                dataout(["% 18.10E" % x for x in data2])
                if sens:
                    sim.sensitivities().T.tofile(sensfp)
                told = tnow
                tmpold = r.T
            outlog("  % 14.6E % 14.6E % 14.6E % 14.6E % 14.6E\n" % tuple(data1))
        
        outlog("\n\n  Integration completed.\n")
        if igntime is not None:
            outlog("  Ignition Time = % 14.6E s\n" % (igntime))
            outlog("  Temp criteria = % 14.6E K\n" % (ignT))
        outlog("\n  Output file has %d time datasets.\n\n" % (ndata))
        outlog("  CPU time     : % 12.3f s\n" % (time.clock() - clock0))
        outlog("  Time elapsed : % 12.3f s\n" % (time.time() - time0))

        if sens: sensfp.close()
        return outname


    ## TODO: other ignition criteria (dT/dt, CH max, etc...)
    def run_sens_ign_brute_force(self, tmax, conP=False, energy="on", rtol=1e-9, atol=1e-15,
                                 ignDT=400., diff=0.05, fwd_only=False, select=None, nodat=False, outname=None):
        outname, outlog, dataout, r, sim = \
                self._initsim(outname=outname, conP=conP, energy=energy, rtol=rtol, atol=atol, nodat=nodat)
        
        outlog(" Brute-force sensitivity analysis for ignition time will be performed.\n\n")
        
        # write headers
        idtag = self.idtag
        if idtag.strip() == "": idtag = "(noid)"
        species = list(self.gas.species_names)
        
        dataout(["ZeroD", outname, idtag, "SensIgn"])
        dataout(["Time(s)", "T(K)", "P(Pa)", "V(m3)", "rho(kg/m3)"] + species)
        
        outlog(" Time Integration:\n\n")
        outlog("  % 14s % 14s % 14s % 14s % 14s\n"
                 % ("Time(s)", "T(K)", "P(Pa)", "V(m3)", "rho(kg/m3)"))
        
        clock0 = time.clock()
        time0 = time.time()

        state = self.gas.TDY
        
        # t = 0
        data1 = [sim.time, r.T, r.thermo.P, r.volume, r.density]
        data2 = data1 + list(self.gas.X)
        dataout(["% 18.10E" % x for x in data2])
        outlog("  % 14.6E % 14.6E % 14.6E % 14.6E % 14.6E\n" % tuple(data1))
        data3 = data1 + list(self.gas.Y)
        data0 = [[x] for x in data3]
        
        # integrate
        ignT = r.T + ignDT
        igntime = None
        ndata = 1
        delt = tmax / 100.
        tnow = 0.0
        tnext = 0.0
        told = 0.0
        tmpold = r.T
        while tnow < tmax:
            tnext += delt
            while tnow < tnext:
                tnow = sim.step()
                ndata += 1
                if igntime is None and r.T > ignT:
                    igntime = told + (tnow - told) * (ignT - tmpold) / (r.T - tmpold)
                data1 = [sim.time, r.T, r.thermo.P, r.volume, r.density]
                data2 = data1 + list(self.gas.X)
                dataout(["% 18.10E" % x for x in data2])
                told = tnow
                tmpold = r.T
            outlog("  % 14.6E % 14.6E % 14.6E % 14.6E % 14.6E\n" % tuple(data1))
        
        outlog("\n\n  Integration completed.\n")
        if igntime is not None:
            outlog("  Ignition Time = % 14.6E s\n" % (igntime))
            outlog("  Temp criteria = % 14.6E K\n" % (ignT))
        outlog("\n  Output file has %d time datasets.\n\n" % (ndata))
        outlog("  CPU time     : % 12.3f s\n" % (time.clock() - clock0))
        outlog("  Time elapsed : % 12.3f s\n" % (time.time() - time0))
        time1 = time.time()

        if igntime is None: raise ValueError("No ignition.")
        igntime0 = igntime
        
        # sens
        outlog("\n Brute-force sensitivity analysis for ignition time\n\n")
        
        sensfn = outname + ".sens_ign"
        outlog("  Sensitivity output file: %s\n" % sensfn)
        
        outlog("\n")
        if fwd_only: outlog("  Forward differences: + %.2f %%\n" % (diff*100.))
        else: outlog("  Central differences: +/- %.2f %%\n" % (diff*100.))

        # do brute-force
        if select is None:
            nr = self.gas.n_reactions
            rids = list(range(nr))
        else:
            nr = len(select)
            rids = select[:]
        reactions = self.gas.reaction_equations(rids)
        outlog("  Number of reactions: %d\n\n" % (nr))

        def run1(rid, factor):
            self.gas.set_multiplier(1.0*factor, rid)
            self.gas.TDY = state
            if conP:
                r = ct.IdealGasConstPressureReactor(self.gas)
            else:
                r = ct.IdealGasReactor(self.gas)
            sim = ct.ReactorNet([r])
            sim.rtol = rtol
            sim.atol = atol
            sim.set_initial_time(0.0)
            
            igntime = None
            delt = tmax / 100.
            tnow = 0.0
            tnext = 0.0
            told = 0.0
            tmpold = r.T
            while tnow < tmax:
                tnext += delt
                while tnow < tnext:
                    tnow = sim.step()
                    if r.T > ignT:
                        igntime = told + (tnow - told) * (ignT - tmpold) / (r.T - tmpold)
                        break
                    told = tnow
                    tmpold = r.T
                if igntime is not None: break
            if igntime is None: raise ValueError("No ignition.")
            self.gas.set_multiplier(1.0)
            return igntime

        sens_ign = []
        for i in range(nr):
            time2 = time.time()
            rid = rids[i]
            outlog("  R%04d: %s\n" % (rid, reactions[i]))
            outlog("  Incrementing rate constant for R%04d by + %.2f %%\n" % (rid, diff*100.))
            igntime0f = run1(rid, 1. + diff)
            outlog("    Ignition Time = % 14.6E s\n" % (igntime0f))
            outlog("    Time elapsed : % 12.3f s\n" % (time.time() - time2))
            if not fwd_only:
                time2 = time.time()
                outlog("  Incrementing rate constant for R%04d by - %.2f %%\n" % (rid, diff*100.))
                igntime0b = run1(rid, 1. - diff)
                outlog("    Ignition Time = % 14.6E s\n" % (igntime0b))
                outlog("    Time elapsed : % 12.3f s\n" % (time.time() - time2))
            outlog("\n")

            if fwd_only: val = (igntime0f - igntime0) / igntime0 / diff
            else: val = (igntime0f - igntime0b) / igntime0 / (2.*diff)
            sens_ign.append((rid, val, reactions[i]))
            
            if i % 10 == 0 and i != 0:
                outlog("  %d / %d (%.1F %%) steps done.\n"
                       % (i, nr, 100.*i/float(nr)))
                outlog("  Total CPU time     : % 12.3f s\n" % (time.clock() - clock0))
                outlog("  Total time elapsed : % 12.3f s\n" % (time.time() - time0))
                outlog("  Estimated remaining time : % 12.3f s\n"
                       % ((nr - i)/float(i) * (time.time() - time1)))
                outlog("\n")

        outlog("  Writing sensitivity coefficients to %s\n" % sensfn)
        sens_ign = sorted(sens_ign, key=lambda x: abs(x[1]), reverse=True)
        sfp = open(sensfn, "w")
        sfp.write("Sensitivity coefficients (sorted) to the ignition time\n")
        sfp.write("  id                S    Reaction\n")
        for rid, val, rxn in sens_ign:
            sfp.write("% 4d   % 14.6E    %s\n" % (rid, val, rxn))
        outlog("  ... Done.\n\n")
        
        outlog("  Brute-force sensitivity analysis for ignition time completed.\n")
        outlog("  Total CPU time     : % 12.3f s\n" % (time.clock() - clock0))
        outlog("  Total time elapsed : % 12.3f s\n" % (time.time() - time0))
        
        return outname



class OneDFree (SimBase):

    def __init__(self, chemfn="chem.inp", thermfn="thermo.dat", tranfn="tran.dat",
                 idtag="", logfp=None):
        self.flame = None
        SimBase.__init__(self, chemfn=chemfn, thermfn=thermfn, tranfn=tranfn,
                         idtag=idtag, logfp=logfp)

    def init_flame(self, T, P, comp, initial_grid, vel):
        if isinstance(comp, dict):
            X = ", ".join(["%s:%g" % (nm, comp[nm]) for nm in comp])
        else:
            X = comp
        self.log(" T = %g K\n" % T)
        self.log(" P = %g atm\n" % P)
        self.log(" Inlet mole fractions = %s\n\n" % X)
        self.gas.TPX = T, P*ct.one_atm, X

        mdot = vel * self.gas.density      # m/s * kg/m3 = kg/m^2/s
        self.log(" u = %g m/s\n" % vel)
        self.log(" Grid = %s [m]\n" % initial_grid)
        self.flame = FreeFlame(gas=self.gas, grid=initial_grid)
        self.flame.inlet.mdot = mdot
        self.flame.inlet.X = X
        self.flame.inlet.T = T
        #self.flame.init()
        self.inletTPX = T, P, X
        return

    def update_flame(self, T, P, comp):
        if isinstance(comp, dict):
            X = ", ".join(["%s:%g" % (nm, comp[nm]) for nm in comp])
        else:
            X = comp
        self.log(" T = %g K\n" % T)
        self.log(" P = %g atm\n" % P)
        self.log(" Inlet mole fractions = %s\n\n" % X)
        self.gas.TPX = T, P*ct.one_atm, X
        self.flame.inlet.X = X
        self.flame.inlet.T = T
        self.inletTPX = T, P, X
        return
        

    def _initsim(self, outname, append=False):
        # files
        if outname is None:
            count = 0
            while True:
                count += 1
                outname = self.name + ("_free%03d" % count)
                if not os.path.exists(outname+".out"): break
        outfn = outname + ".out"
        origctmlfn = self.srcfn
        ctmlfn = outname + ".xml"
        shutil.copy(origctmlfn, ctmlfn)
        
        self.log("\n\n================================================================\n\n")
        self.log(" Output: %s\n" % outfn)
        self.log(" CTML file was copied from %s to %s\n\n" % (origctmlfn, ctmlfn))

        if append: ofp = open(outfn, "a")
        else: ofp = open(outfn, "w")
        def outlog(s):
            self.log(s)
            ofp.write(s)
            ofp.flush()
        
        outlog(" Date: %s\n\n" % time.strftime('%Y-%m-%d %H:%M:%S'))
        return outname, outlog


    def solve(self, outname=None, append=False, loglevel=3, refine_grid=1):
        outname, outlog = self._initsim(outname=outname, append=append)
        
        clock0 = time.clock()
        time0 = time.time()

        # TODO output initial condition to outlog
        
        self.flame.solve(loglevel, refine_grid)
    
        datfn = outname + ".dat"
        dfp = open(datfn, "w")
        def dataout(l):
            for x in l: dfp.write("%s " % x)
            dfp.write("\n")
            dfp.flush()
        
        # write headers
        idtag = self.idtag
        if idtag.strip() == "": idtag = "(noid)"
        species = list(self.gas.species_names)
        
        dataout(["OneDFree", outname, idtag])
        dataout(["z(m)", "T(K)", "P(Pa)", "u(m/s)", "V(1/s)", "rho(kg/m3)"] + species)
        
        z = self.flame.flame.grid
        T = self.flame.T
        u = self.flame.u
        V = self.flame.V
        for n in range(self.flame.flame.n_points):
            self.flame.set_gas_state(n)
            data = [z[n], T[n], self.gas.P, u[n], V[n], self.gas.density] \
                   + list(self.gas.X)
            dataout(["% 18.10E" % x for x in data])
        
        outlog("\n  Problem solved with  %d grid points.\n" % (len(self.flame.flame.grid)))
        outlog("  Flame speed = % 14.6E m/s\n" % (self.flame.u[0]))
        outlog("  Full data output: %s\n" % datfn)
        outlog("  CPU time     : % 12.3f s\n" % (time.clock() - clock0))
        outlog("  Time elapsed : % 12.3f s\n" % (time.time() - time0))
        
        return outname

    def get_flame_deltaT(self):
        z = self.flame.flame.grid
        T = self.flame.T
        dtmax = (0.0, 0.0)
        for i in range(1,len(z)-1):
            dt = (T[i+1] - T[i-1]) / (z[i+1] - z[i-1])
            if dt > dtmax[0]: 
                dtmax = (dt, z[i])
        deltaT = (T[-1] - T[0]) / dtmax[0]
        return deltaT

    def solve_sens_u0_brute_force(self, outname=None, loglevel=3, refine_grid=1,
                                  diff=0.05, fwd_only=False, select=None):
        outname, outlog = self._initsim(outname=outname)
        
        outlog(" Brute-force sensitivity analysis for flame speed will be performed.\n\n")
        
        clock0 = time.clock()
        time0 = time.time()
        
        self.flame.solve(loglevel, refine_grid)
    
        datfn = outname + ".dat"
        dfp = open(datfn, "w")
        def dataout(l):
            for x in l: dfp.write("%s " % x)
            dfp.write("\n")
            dfp.flush()
        
        # write headers
        idtag = self.idtag
        if idtag.strip() == "": idtag = "(noid)"
        species = list(self.gas.species_names)
        
        dataout(["OneDFree", outname, idtag, "Sensu0"])
        dataout(["z(m)", "T(K)", "P(Pa)", "u(m/s)", "V(1/s)", "rho(kg/m3)"] + species)
        
        z = self.flame.flame.grid
        T = self.flame.T
        u = self.flame.u
        V = self.flame.V
        for n in range(self.flame.flame.n_points):
            self.flame.set_gas_state(n)
            data = [z[n], T[n], self.gas.P, u[n], V[n], self.gas.density] \
                   + list(self.gas.X)
            dataout(["% 18.10E" % x for x in data])
            
        outlog("\n  Problem solved with  %d grid points.\n" % (len(self.flame.flame.grid)))
        outlog("  Flame speed = % 14.6E m/s\n" % (self.flame.u[0]))
        outlog("  Full data output: %s\n" % datfn)
        outlog("  CPU time     : % 12.3f s\n" % (time.clock() - clock0))
        outlog("  Time elapsed : % 12.3f s\n" % (time.time() - time0))
        time1 = time.time()
        dfp.close()
        
        # sens
        
        outlog("\n Brute-force sensitivity analysis for flame speed\n\n")
        
        sensfn = outname + ".sens_u0"
        outlog("  Sensitivity output file: %s\n" % sensfn)
        outlog("\n")
        outlog(" Re-computing flame speed.\n")
        
        # wrapper
        def run1(rid, factor):
            if rid is None:
                self.gas.set_multiplier(1.0)
            else:
                self.gas.set_multiplier(1.0*factor, rid)
            self.flame.solve(loglevel, 0)
            outlog("  Flame speed = % 14.6E m/s\n" % (self.flame.u[0]))
            self.gas.set_multiplier(1.0)
            return self.flame.u[0]
        
        # Re-compute flame speed using soln file
        u0 = run1(None, None)
        
        if fwd_only: outlog("  Forward differences: + %.2f %%\n" % (diff*100.))
        else: outlog("  Central differences: +/- %.2f %%\n" % (diff*100.))

        # do brute-force
        if select is None:
            nr = self.gas.n_reactions
            rids = list(range(nr))
        else:
            nr = len(select)
            rids = select[:]
        reactions = self.gas.reaction_equations(rids)
        outlog("  Number of reactions: %d\n\n" % (nr))
        
        
        sens_u0 = []
        for i in range(nr):
            rid = rids[i]
            outlog("  R%04d: %s\n" % (rid, reactions[i]))
            outlog("  Incrementing rate constant for R%04d by + %.2f %%\n" % (rid, diff*100.))
            u0f = run1(rid, 1. + diff)
            if not fwd_only:
                outlog("  Incrementing rate constant for R%04d by - %.2f %%\n" % (rid, diff*100.))
                u0b = run1(rid, 1. - diff)
            outlog("\n")

            if fwd_only: val = (u0f - u0) / u0 / diff
            else: val = (u0f - u0b) / u0 / (2.*diff)
            sens_u0.append((rid, val, reactions[i]))

            if i % 50 == 0 and i != 0:
                outlog("  %d / %d (%.1F %%) steps done.\n"
                       % (i, nr, 100.*i/float(nr)))
                outlog("  Total CPU time     : % 12.3f s\n" % (time.clock() - clock0))
                outlog("  Total time elapsed : % 12.3f s\n" % (time.time() - time0))
                outlog("  Estimated remaining time : % 12.3f s\n"
                       % ((nr - i)/float(i) * (time.time() - time1)))
                outlog("\n")
        
        outlog("  Writing sensitivity coefficients to %s\n" % sensfn)
        sens_u0 = sorted(sens_u0, key=lambda x: abs(x[1]), reverse=True)
        sfp = open(sensfn, "w")
        sfp.write("Sensitivity coefficients (sorted) to the flame speed\n")
        sfp.write("  id                S    Reaction\n")
        for rid, val, rxn in sens_u0:
            sfp.write("% 4d   % 14.6E    %s\n" % (rid, val, rxn))
        outlog("  ... Done.\n\n")
        
        outlog("  Brute-force sensitivity analysis for flame speed completed.\n")
        outlog("  Total CPU time     : % 12.3f s\n" % (time.clock() - clock0))
        outlog("  Total time elapsed : % 12.3f s\n" % (time.time() - time0))
        
        return outname


class OneDCDF (SimBase):

    def __init__(self, chemfn="chem.inp", thermfn="thermo.dat", tranfn="tran.dat",
                 idtag="", logfp=None):
        self.flame = None
        SimBase.__init__(self, chemfn=chemfn, thermfn=thermfn, tranfn=tranfn,
                         idtag=idtag, logfp=logfp)

    def init_flame(self, T_o, T_f, P, comp_o, comp_f, initial_grid, vel_o,
                   vel_f):
        if isinstance(comp_o, dict):
            X_o = ", ".join(["%s:%g" % (nm, comp_o[nm]) for nm in comp_o])
        else:
            X_o = comp_o
        if isinstance(comp_f, dict):
            X_f = ", ".join(["%s:%g" % (nm, comp_f[nm]) for nm in comp_f])
        else:
            X_f = comp_f
        self.log(" T(oxidizer inlet) = %g K\n" % T_o)
        self.log(" T(fuel inlet) = %g K\n" % T_f)
        self.log(" P = %g atm\n" % P)
        self.log(" Inlet oxidizer mole fractions = %s\n" % X_o)
        self.log(" Inlet fuel mole fractions = %s\n\n" % X_f)
        self.gas.TPX = T_o, P*ct.one_atm, X_o
        mdot_o = vel_o * self.gas.density
        self.gas.TPX = T_f, P*ct.one_atm, X_f
        mdot_f = vel_f * self.gas.density
        self.log(" u(oxidizer) = %g m/s\n" % vel_o)
        self.log(" u(fuel) = %g m/s\n" % vel_f)
        self.log(" Grid = %s [m]\n" % initial_grid)
        self.flame = CounterflowDiffusionFlame(gas=self.gas, grid=initial_grid)
        self.flame.fuel_inlet.mdot = mdot_f
        self.flame.fuel_inlet.X = X_f
        self.flame.fuel_inlet.T = T_f
        self.flame.oxidizer_inlet.mdot = mdot_o
        self.flame.oxidizer_inlet.X = X_o
        self.flame.oxidizer_inlet.T = T_o
        return

    def _initsim(self, outname):
        # files
        if outname is None:
            count = 0
            while True:
                count += 1
                outname = self.name + ("_cdf%03d" % count)
                if not os.path.exists(outname+".out"): break
        outfn = outname + ".out"
        origctmlfn = self.srcfn
        ctmlfn = outname + ".xml"
        shutil.copy(origctmlfn, ctmlfn)

        self.log("\n\n================================================================\n\n")
        self.log(" Output: %s\n" % outfn)
        self.log(" CTML file was copied from %s to %s\n\n" % (origctmlfn, ctmlfn))
        
        ofp = open(outfn, "w")
        def outlog(s):
            self.log(s)
            ofp.write(s)
            ofp.flush()
        
        outlog(" Date: %s\n\n" % time.strftime('%Y-%m-%d %H:%M:%S'))
        return outname, outlog

    def solve(self, outname=None, loglevel=3, refine_grid=1):
        outname, outlog = self._initsim(outname=outname)

        clock0 = time.clock()
        time0 = time.time()

        # TODO output initial condition to outlog

        self.flame.solve(loglevel=loglevel, refine_grid=refine_grid)

        datfn = outname + ".dat"
        dfp = open(datfn, "w")
        def dataout(l):
            for x in l: dfp.write("%s " % x)
            dfp.write("\n")
            dfp.flush()
        
        # write headers
        idtag = self.idtag
        if idtag.strip() == "": idtag = "(noid)"
        species = list(self.gas.species_names)
        
        dataout(["OneDCDF", outname, idtag])
        dataout(["z(m)", "T(K)", "P(Pa)", "u(m/s)", "V(1/s)", "rho(kg/m3)"] + species)
        
        z = self.flame.flame.grid
        T = self.flame.T
        u = self.flame.u
        V = self.flame.V
        for n in range(self.flame.flame.n_points):
            self.flame.set_gas_state(n)
            data = [z[n], T[n], self.gas.P, u[n], V[n], self.gas.density] \
                   + list(self.gas.X)
            dataout(["% 18.10E" % x for x in data])
        
        outlog("\n  Problem solved with  %d grid points.\n" % (len(self.flame.flame.grid)))
#        outlog("  Flame speed = % 14.6E m/s\n" % (self.flame.u[0]))
        outlog("  Full data output: %s\n" % datfn)
        outlog("  CPU time     : % 12.3f s\n" % (time.clock() - clock0))
        outlog("  Time elapsed : % 12.3f s\n" % (time.time() - time0))
        
        return outname

    def get_flame_deltaT(self):
        z = self.flame.flame.grid
        T = self.flame.T
        dtmax = (0.0, 0.0)
        for i in range(1,len(z)-1):
            dt = (T[i+1] - T[i-1]) / (z[i+1] - z[i-1])
            if dt > dtmax[0]: 
                dtmax = (dt, z[i])
        deltaT = (T[-1] - T[0]) / dtmax[0]
        return deltaT

    
def run0d(chemfn="chem.inp", thermfn="thermo.dat", idtag="", logfp=None,
           T=None, P=None, reac=None, 
           tmax=None, conP=False, energy="on", sens=False,
          rtol=1e-9, atol=1e-15, rtols=1e-6, atols=1e-9, ignDT=400., outname=None):
    zd = ZeroD(chemfn=chemfn, thermfn=thermfn, idtag=idtag, logfp=logfp)
    zd.set_gas(T=T, P=P, reac=reac)
    zd.run(tmax=tmax, conP=conP, energy=energy, sens=sens, rtol=rtol, atol=atol, 
           rtols=rtols, atols=atols, ignDT=ignDT, outname=outname)
    return zd

    
def run0d_sens_ign(chemfn="chem.inp", thermfn="thermo.dat", idtag="", logfp=None,
                   T=None, P=None, reac=None, 
                   tmax=None, conP=False, energy="on", rtol=1e-9, atol=1e-15,
                   ignDT=400., diff=0.05, fwd_only=False, select=None, nodat=False, outname=None):
    zd = ZeroD(chemfn=chemfn, thermfn=thermfn, idtag=idtag, logfp=logfp)
    zd.set_gas(T=T, P=P, reac=reac)
    zd.run_sens_ign_brute_force(tmax=tmax, conP=conP, energy=energy, rtol=rtol, atol=atol, ignDT=ignDT,
                                diff=diff, fwd_only=fwd_only, select=select, nodat=nodat, outname=outname)
    return zd
    



if __name__ == "__main__":
    pass
