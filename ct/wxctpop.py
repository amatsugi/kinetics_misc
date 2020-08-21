#! /usr/bin/env python3

import os
import sys
import time
import traceback

import wx
from wx.lib.filebrowsebutton import FileBrowseButton, DirBrowseButton
from wx.lib.mixins.listctrl import ColumnSorterMixin
from wx.lib.inspection import InspectionTool
from wx.lib.scrolledpanel import ScrolledPanel

import numpy as np

from wxmpl import PlotFrame, genColorAndStyles
from wxflux import FluxFrame
from ctpop import CanteraOutput

# some sub panels

class UnitChoices(wx.Panel):

    def __init__(self, parent, units,
                 orient=wx.HORIZONTAL, flags=wx.ALL, border=4,
                 staticbox=False, staticbox_label=''):
        wx.Panel.__init__(self, parent)
        # units: {"lavel": [unit1, unit2,...], ...}
        self.units = units
        self.choices = {}
        
        box = wx.StaticBoxSizer(wx.StaticBox(self, label=staticbox_label), orient) \
              if staticbox else wx.BoxSizer(orient)
        
        for label in self.units:
            choice = wx.Choice(self, choices=self.units[label])
            choice.SetSelection(0)
            box.Add(choice, 0, flags, border)
            self.choices[label] = choice
        
        self.SetSizer(box)

    def getString(self, label):
        return self.choices[label].GetStringSelection()
        
    def getCurrentSelection(self, label):
        return self.choices[label].GetCurrentSelection()


class SpeciesChoice(wx.Panel):

    def __init__(self, parent, species,
                 orient=wx.HORIZONTAL, flags=wx.ALL, border=4,
                 staticbox=False, staticbox_label='', add_num=False, temp=False, init=0):
        wx.Panel.__init__(self, parent)
        self.species = list(species)
        if add_num:
            if temp: labels = ['%d: %s' % (i, x) for i,x in enumerate(species)]
            else: labels = ['%d: %s' % (i+1, x) for i,x in enumerate(species)]
            
        else: labels = self.species[:]
        self.choices = []
        
        box = wx.StaticBoxSizer(wx.StaticBox(self, label=staticbox_label), orient) \
              if staticbox else wx.BoxSizer(orient)
        
        self.choice = wx.Choice(self, choices=labels)
        self.choice.SetSelection(0)
        box.Add(self.choice, 0, flags, border)
        self.SetSizer(box)

    def getCurrentString(self):
        return self.species[self.choice.GetCurrentSelection()]
        
    def getCurrentIndex(self):
        return self.choice.GetCurrentSelection()
        

class MultiChoices(wx.Panel):

    def __init__(self, parent, choices, orient=wx.VERTICAL, flags=wx.ALL, border=4,
                 staticbox=False, staticbox_label='', add_num=False, top=None,
                 init=None, small=None):
        """init: True for check all, list of indices for check them, None for none"""
        wx.Panel.__init__(self, parent)
        self.choices = choices
        if add_num: labels = ['%d: %s' % (i+1, x) for i,x in enumerate(self.choices)]
        else: labels = self.choices[:]
        self.top = top
        self.small = small
        
        self.checklist = wx.CheckListBox(self, choices=labels)
        if init is True: self.checkAll()
        elif init is not None: self.checkSome(init)
        
        self.Bind(wx.EVT_CHECKLISTBOX, self.onCheck, self.checklist)
        self.checkallbutton = wx.Button(self, label='Check all')
        self.Bind(wx.EVT_BUTTON, self.onCheckAll, self.checkallbutton)
        self.uncheckallbutton = wx.Button(self, label='Uncheck all')
        self.Bind(wx.EVT_BUTTON, self.onUncheckAll, self.uncheckallbutton)
        if top is not None:
            self.checktopbutton = wx.Button(self, label='Check top')
            self.Bind(wx.EVT_BUTTON, self.onCheckTop, self.checktopbutton)
            self.top_spin = wx.SpinCtrl(self, size=(60, -1),
                                        min=min(len(self.choices), 1),
                                        max=len(self.choices),
                                        initial=min(len(self.choices), 10))
        if small is not None:
            self.checksmassbutton = wx.Button(self, label='Check C0-C1')
            self.Bind(wx.EVT_BUTTON, self.onCheckSmall, self.checksmassbutton)
            
        buttonbox = wx.BoxSizer(wx.VERTICAL)
        buttonbox.Add(self.checkallbutton, 0, wx.TOP | wx.BOTTOM, border)
        buttonbox.Add(self.uncheckallbutton, 0, wx.TOP | wx.BOTTOM, border)
        if top is not None:
            topbox = wx.BoxSizer(wx.HORIZONTAL)
            topbox.Add(self.checktopbutton, 0, border=0)
            topbox.Add(self.top_spin, 0, wx.LEFT, border)
            buttonbox.Add(topbox, 0, wx.TOP | wx.BOTTOM, border)
        if small is not None:
            buttonbox.Add(self.checksmassbutton, 0, wx.TOP | wx.BOTTOM, border)
        box = wx.StaticBoxSizer(wx.StaticBox(self, label=staticbox_label), orient) \
              if staticbox else wx.BoxSizer(orient)
        box.Add(self.checklist, 0, flags, border)
        box.Add(buttonbox, 0, flags, border)
        
        self.SetSizer(box)
        
    def onCheck(self, event):
        self.checklist.SetSelection(event.GetSelection())
        
    def onCheckAll(self, event):
        self.checkAll()
            
    def onUncheckAll(self, event):
        self.checkAll(False)
        
    def onCheckTop(self, event):
        self.checkSome(self.top(self.top_spin.GetValue()), False)
        
    def onCheckSmall(self, event):
        self.checkSome(self.small(), False)

    def checkAll(self, check=True):
        for i in range(self.checklist.Count):
            self.checklist.Check(i, check=check)

    def checkSome(self, checks=None, uncheck=True):
        """checks: list of indices"""
        if checks is None: checks = []
        for i in range(self.checklist.Count):
            if i in checks: self.checklist.Check(i, True)
            elif uncheck: self.checklist.Check(i, False)

    def getCheckedIndices(self):
        return [i for i in range(self.checklist.Count) \
                if self.checklist.IsChecked(i)]

    def getCheckedStrings(self):
        return [self.choices[i] for i in range(self.checklist.Count) \
                if self.checklist.IsChecked(i)]


class XvarSlider(wx.Panel):

    def __init__(self, parent, ctout, init=None):
        wx.Panel.__init__(self, parent)
        self.ctout = ctout
        if init is None:
            init = int((self.ctout.ndata-1)/2)
        self.index_text = wx.StaticText(self, label='Select %s' % ctout.xvarname)
        self.index_slider = wx.Slider(self, minValue=0, maxValue=self.ctout.ndata-1,
                                      value=init, style=wx.SL_HORIZONTAL)
        self.Bind(wx.EVT_COMMAND_SCROLL, self.onScrollIndex, self.index_slider)
        slider_box = wx.BoxSizer(wx.VERTICAL)
        slider_box.Add(self.index_text, 0, border=0)
        slider_box.Add(self.index_slider, 0, border=0)
        self.SetSizer(slider_box)
        
    def onScrollIndex(self, event):
        index = self.index_slider.GetValue()
        xvar = self.ctout.xvarAt(index)
        self.index_text.SetLabel('%d: %.4e' % (index, xvar))

    def getIndex(self):
        return self.index_slider.GetValue()

    def getXvar(self):
        return self.ctout.xvarAt(self.index_slider.GetValue())


# panels for notebook pages

class InfoPanel(ScrolledPanel):
    
    def __init__(self, parent, mainframe):
        ScrolledPanel.__init__(self, parent)
        self.ctout = mainframe.ctout
        self.message = mainframe.message
        self.create()
        self.SetupScrolling()

    def create(self):
        
        info_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Info'), wx.VERTICAL)
        text = '\n'.join(self.ctout.info)
        info_vbox.Add(wx.StaticText(self, label=text), 0, wx.ALL, 4)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(info_vbox, 0, wx.ALL, 4)
        
        self.SetSizer(vbox)

        
class ConcPanel(wx.Panel, ColumnSorterMixin):
    
    def __init__(self, parent, mainframe):
        wx.Panel.__init__(self, parent)
        self.ctout = mainframe.ctout
        self.message = mainframe.message
        
        self.create()
        self.updateList()
        ColumnSorterMixin.__init__(self, 4)
        
    def create(self):
        self.list = wx.ListCtrl(self, style=wx.LC_REPORT)
        self.list.InsertColumn(0, "ID", width=40)
        self.list.InsertColumn(1, "Speceis", width=100)
        self.list.InsertColumn(2, "Weight", wx.LIST_FORMAT_RIGHT, width=60)
        self.list.InsertColumn(3, "Concencration", wx.LIST_FORMAT_RIGHT, width=100)

        self.xvar_slider = XvarSlider(self, self.ctout, init=0)
        self.unit_choices = UnitChoices(self, {'conc': ["molefrac", "massfrac", "molecules/cm3", "mole/cm3"]}, wx.VERTICAL)
        self.update_button = wx.Button(self, label='Update')
        self.Bind(wx.EVT_BUTTON, self.onUpdate, self.update_button)
        self.plot_button = wx.Button(self, label='Plot')
        self.Bind(wx.EVT_BUTTON, self.onPlot, self.plot_button)
        self.copy_button = wx.Button(self, label='Copy')
        self.Bind(wx.EVT_BUTTON, self.onCopy, self.copy_button)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.xvar_slider, 0, wx.ALL, 4)
        vbox.Add(self.unit_choices, 0, border=0)
        vbox.Add(self.update_button, 0, wx.ALL, 4)
        vbox.Add(self.plot_button, 0, wx.ALL, 4)
        vbox.Add(self.copy_button, 0, wx.ALL, 4)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(vbox, 0, wx.ALL, 4)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        
    # Used by the ColumnSorterMixin
    def GetListCtrl(self):
        return self.list

    def onUpdate(self, event):
        self.updateList()

    def onPlot(self, event):
        self.updateList()
        self.plot()

    def onCopy(self, event):
        self.copy()

    def updateList(self):
        index = self.xvar_slider.getIndex()
        conc_unit = self.unit_choices.getString('conc')
        xvar, press, temp, concs, addv = \
              self.ctout.getProfileAt(index, conc_unit=conc_unit)
        
        # ColumnSorterMixin looks 'itemDataMap' variable to sort items 
        self.itemDataMap = {}
        for i in range(self.ctout.ns):
            self.itemDataMap[i+1] = (i+1, self.ctout.species[i], 
                                     self.ctout.weights[i], concs[i])
            
        self.list.DeleteAllItems()
        for key, data in self.itemDataMap.items():
            index = self.list.InsertItem(self.ctout.ns, '%d' % (data[0],))
            self.list.SetItem(index, 1, data[1])
            self.list.SetItem(index, 2, '%.2f' % (data[2],))
            self.list.SetItem(index, 3, '%.3e' % (data[3],))
            self.list.SetItemData(index, key)

    def plot(self):
        index = self.xvar_slider.getIndex()
        conc_unit = self.unit_choices.getString('conc')
        xvar, press, temp, concs, addv = \
              self.ctout.getProfileAt(index, conc_unit=conc_unit)
        
        self.message('Plotting mass spectrum at %s...' % (xvar,))
        frame = PlotFrame(self, title='Mass Spectrum at %s' % (xvar,))
        ax = frame.addSubplot(111, xlabel="m",
                              ylabel='Concentration / %s' % (conc_unit,))
        d = {}
        for i in range(self.ctout.ns):
            intmass = int(round(self.ctout.weights[i]))
            if intmass in d: d[intmass] += concs[i]
            else: d[intmass] = concs[i]
        x, y = [], []
        for intmass in d:
            x.append(intmass)
            y.append(d[intmass])
        if len(y) > 20: topy = sorted(y)[-20]
        else: topy = 0.0
        
        ax.axhline(color='k', linewidth=1)
        ax.vlines(x, 0., y, color='k', linestyles='solid', linewidth=2)
        for i in range(len(x)):
            if y[i] > topy:
                ax.text(x[i], 1.05*y[i], "%d" % x[i], ha='center', va='bottom')
        
        frame.draw()
        frame.Show()
        self.message('Plotting mass spectrum at %s... done.' % (xvar,))

    def copy(self):
        if not wx.TheClipboard.Open():
            wx.MessageBox("failed to open clipboard.")
            return
        text = ""
        for i in range(self.ctout.ns):
            d = self.itemDataMap[self.list.GetItemData(i)]
            first = True
            for x in d:
                if first:
                    text += "%s" % x
                    first = False
                else:
                    text += "\t%s" % x
            text += "\n"
        wx.TheClipboard.SetData(wx.TextDataObject(text))
        wx.TheClipboard.Close()
        return

        
class RatesPanel(wx.Panel, ColumnSorterMixin):
    
    def __init__(self, parent, mainframe):
        wx.Panel.__init__(self, parent)
        self.ctout = mainframe.ctout
        self.message = mainframe.message

        self.create()
        self.updateList()
        ColumnSorterMixin.__init__(self, 3)
        
    def create(self):
        self.list = wx.ListCtrl(self, style=wx.LC_REPORT)
        self.list.InsertColumn(0, "ID", width=60)
        self.list.InsertColumn(1, "Reactions", width=280)
        self.list.InsertColumn(2, "Rates", wx.LIST_FORMAT_RIGHT, width=100)

        self.xvar_slider = XvarSlider(self, self.ctout, init=0)
        self.unit_choices = UnitChoices(self, {'rate': ["molecules/cm3-s", "mole/cm3-s"]}, wx.VERTICAL)
        self.update_button = wx.Button(self, label='Update')
        self.Bind(wx.EVT_BUTTON, self.onUpdate, self.update_button)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.xvar_slider, 0, wx.ALL, 4)
        vbox.Add(self.unit_choices, 0, border=0)
        vbox.Add(self.update_button, 0, wx.ALL, 4)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(vbox, 0, wx.ALL, 4)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        
    def GetListCtrl(self):
        return self.list

    def onUpdate(self, event):
        self.updateList()

    def updateList(self):
        index = self.xvar_slider.getIndex()
        rate_unit = self.unit_choices.getString('rate')
        rates = [float(x) for x in self.ctout.getReactionRatesAt(index, unit=rate_unit)]
        
        self.itemDataMap = {}
        for i in range(self.ctout.nr):
            self.itemDataMap[i+1] = (i+1,  self.ctout.reactions[i], rates[i])

        self.list.DeleteAllItems()
        for key, data in self.itemDataMap.items():
            index = self.list.InsertItem(self.ctout.nr, '%d' % (data[0],))
            self.list.SetItem(index, 1, data[1])
            self.list.SetItem(index, 2, '%.3e' % (data[2],))
            self.list.SetItemData(index, key)


class ProfilePanel(ScrolledPanel):
    
    header_list = ('None', 'Add', 'With "#"')
    header_values = (None, '', '#')
    header_default_index = 2
    wildcard = 'Space separated files (*.dat)|*.dat|' \
               'Comma separated files (*.csv)|*.csv|' \
               'Tab separated files (*.tsv)|*.tsv|' \
               'All files (space separated) (*.*)|*.*'
    sep_values = (' ', ', ', '\t', ' ')
    
    def __init__(self, parent, mainframe):
        ScrolledPanel.__init__(self, parent)
        self.ctout = mainframe.ctout
        self.message = mainframe.message
        self.create()
        self.SetupScrolling()

    def create(self):
        # configs
        self.unit_choices = UnitChoices(self, {'conc': ["molefrac", "massfrac", "molecules/cm3", "mole/cm3"]},
                                        wx.VERTICAL, staticbox=True,
                                        staticbox_label='Units')
        
        # for exporting
        self.choice_header = wx.Choice(self, choices=self.header_list)
        self.choice_header.SetSelection(self.header_default_index)
        self.spin_precision = wx.SpinCtrl(self, size=(60, -1),
                                          min=2, max=16, initial=6)
        self.export_button = wx.Button(self, label='Export')
        self.Bind(wx.EVT_BUTTON, self.onExport, self.export_button)
        
        # species
        self.species_choices = MultiChoices(self, self.ctout.species, wx.VERTICAL,
                                            wx.HORIZONTAL, staticbox=True,
                                            staticbox_label='Species',
                                            add_num=True, top=self.getTopFrac)
        
        # plot buttons
        self.conc_plot_button = wx.Button(self, label='Conc')
        self.Bind(wx.EVT_BUTTON, self.onPlotConc, self.conc_plot_button)
        self.press_plot_button = wx.Button(self, label='Press')
        self.Bind(wx.EVT_BUTTON, self.onPlotPress, self.press_plot_button)
        self.temp_plot_button = wx.Button(self, label='Temp')
        self.Bind(wx.EVT_BUTTON, self.onPlotTemp, self.temp_plot_button)
        if self.ctout.naddv > 0:
            self.choice_addv = wx.Choice(self, choices=self.ctout.addvnames)
            self.choice_addv.SetSelection(0)
            self.addv_plot_button = wx.Button(self, label='Addv')
            self.Bind(wx.EVT_BUTTON, self.onPlotAddv, self.addv_plot_button)
        
        # Layout
        btn_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Plot'), wx.VERTICAL)
        btn_vbox.Add(self.conc_plot_button, 0, wx.ALL, 4)
        btn_vbox.Add(self.press_plot_button, 0, wx.ALL, 4)
        btn_vbox.Add(self.temp_plot_button, 0, wx.ALL, 4)
        if self.ctout.naddv > 0:
            btn_vbox.AddSpacer(2)
            btn_vbox.Add(self.choice_addv, 0, wx.ALL, 4)
            btn_vbox.Add(self.addv_plot_button, 0, wx.ALL, 4)
        
        export_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Export'), wx.VERTICAL)
        export_vbox.Add(wx.StaticText(self, label='Header:'), 0, wx.ALL, 4)
        export_vbox.Add(self.choice_header, 0, wx.ALL, 4)
        export_vbox.AddSpacer(2)
        export_vbox.Add(wx.StaticText(self, label='Precision:'), 0, wx.ALL, 4)
        export_vbox.Add(self.spin_precision, 0, wx.ALL, 4)
        export_vbox.AddSpacer(2)
        export_vbox.Add(self.export_button, 0, wx.ALL, 4)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.unit_choices, 0, border=0)
        vbox.AddSpacer(10)
        vbox.Add(export_vbox, 0, border=0)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(vbox, 0, wx.ALL, 4)
        hbox.Add(self.species_choices, 0, wx.ALL, 4)
        hbox.Add(btn_vbox, 0, wx.ALL, 4)
        
        self.SetSizer(hbox)
        
    def onPlotConc(self, event): self.plotConc()
    def onPlotPress(self, event): self.plotPress()
    def onPlotTemp(self, event): self.plotTemp()
    def onPlotAddv(self, event): self.plotAddv()
    def onExport(self, event): self.export()
        
    def plotConc(self):
        self.message('Plotting concentration profiles...')
        conc_unit = self.unit_choices.getString('conc')
        species = self.species_choices.getCheckedStrings()
        if len(species) == 0:
            wx.MessageBox('Select spieces', 'alert', 
                          wx.OK | wx.CENTRE | wx.ICON_EXCLAMATION)
            return
        self.message('Plotting %d species...' % (len(species),))
        xvars = np.asarray(self.ctout.getXvarList())
        concs = list(map(np.asarray, self.ctout.getConcLists(species, conc_unit)))

        frame = PlotFrame(self, title='Concentration profiles')
        frame.setLegend(1)
        ax = frame.addSubplot(111, xlabel=self.ctout.xvarname,
                              ylabel='Concentration / %s' % (conc_unit,))
        
        clrsty = genColorAndStyles()
        for s,c in zip(species, concs):
            clr, sty = next(clrsty)
            ax.plot(xvars, c, color=clr, ls=sty, lw=1.5, label=s)
        frame.draw()
        frame.Show()
        self.message('Plotting %d species... done.' % (len(species),))
        
    def plotPress(self):
        self.message('Plotting pressure profile...')
        xvars = np.asarray(self.ctout.getXvarList())
        presses = np.asarray(self.ctout.getPressList())

        frame = PlotFrame(self, title='Pressure profile')
        frame.setLegend(1)
        ax = frame.addSubplot(111, xlabel=self.ctout.xvarname, ylabel='P(Pa)')
        ax.plot(xvars, presses, lw=1.5, label='Pressure')
        frame.draw()
        frame.Show()
        self.message('Plotting pressure profile... done.')
        
    def plotTemp(self):
        self.message('Plotting temperature profile...')
        xvars = np.asarray(self.ctout.getXvarList())
        temps = np.asarray(self.ctout.getTempList())

        frame = PlotFrame(self, title='Temperature profile')
        frame.setLegend(1)
        ax = frame.addSubplot(111, xlabel=self.ctout.xvarname, ylabel='T(K)')
        ax.plot(xvars, temps, lw=1.5, label='Temperature')
        frame.draw()
        frame.Show()
        self.message('Plotting temperature profile... done.')
        
    def plotAddv(self):
        self.message('Plotting profile...')
        xvars = np.asarray(self.ctout.getXvarList())
        addv = self.ctout.addvnames[self.choice_addv.GetCurrentSelection()]
        addvs = np.asarray(self.ctout.getAddvList(addv))

        frame = PlotFrame(self, title='profile')
        frame.setLegend(1)
        ax = frame.addSubplot(111, xlabel=self.ctout.xvarname, ylabel=addv)
        ax.plot(xvars, addvs, lw=1.5, label=addv)
        frame.draw()
        frame.Show()
        self.message('Plotting profile... done.')

    def export(self):
        inidir = os.path.dirname(self.ctout.datfile)
        dlg = wx.FileDialog(self, message='Export profiles',
                            defaultDir=inidir, defaultFile='profile',
                            wildcard=self.wildcard,
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            #dlg.Destroy()
        else:
            #dlg.Destroy()
            return
            
        sep = self.sep_values[dlg.GetFilterIndex()]
        header = self.header_values[self.choice_header.GetCurrentSelection()]
        precision = self.spin_precision.GetValue(),
        species = self.species_choices.getCheckedStrings()
        
        self.message('Exporting profiles...')
        wx.BeginBusyCursor()
        fp = open(path, 'w')
        self.ctout.exportProfiles(fp, species=False if len(species) == 0 else species,
                                  conc_unit=self.unit_choices.getString('conc'),
                                  precision=precision, sep=sep, header=header)
        fp.close()
        wx.EndBusyCursor()
        self.message('Exporting profiles... done.')

    def getTopFrac(self, num):
        conc_unit = self.unit_choices.getString('conc')
        fracs = self.ctout.getConcLists(self.ctout.species, conc_unit)
        tmp = []
        for i, x in enumerate(fracs): tmp.append((max(x), i))
        tmp = sorted(tmp, reverse=True)
        num = min(num, len(tmp))
        return [x[1] for x in tmp[:num]]

class ProductionPanel(ScrolledPanel):
    
    header_list = ('None', 'Add', 'With "#"')
    header_values = (None, '', '#')
    header_default_index = 2
    wildcard = 'Space separated files (*.dat)|*.dat|' \
               'Comma separated files (*.csv)|*.csv|' \
               'Tab separated files (*.tsv)|*.tsv|' \
               'All files (space separated) (*.*)|*.*'
    sep_values = (' ', ', ', '\t', ' ')
    
    def __init__(self, parent, mainframe):
        ScrolledPanel.__init__(self, parent)
        self.ctout = mainframe.ctout
        self.message = mainframe.message
        self.create()
        self.SetupScrolling()
        
    def create(self):
        self.species_choice = SpeciesChoice(self, species=self.ctout.species,
                                            add_num=True, init=0)
        self.rate_unit_choices = UnitChoices(self, {'rate': ["molecules/cm3-s", "mole/cm3-s"]}, wx.VERTICAL)

        # rate at selected xvar
        self.xvar_slider = XvarSlider(self, self.ctout)
        self.maxnum_spin = wx.SpinCtrl(self, size=(60, -1),
                                       min=min(self.ctout.nr, 5),
                                       max=self.ctout.nr,
                                       initial=min(self.ctout.nr, 10))
        self.shorten_labels_ctrl = wx.CheckBox(self, label="Shorten labels")
        self.shorten_labels_ctrl.SetValue(False)
        self.rates_draw_button = wx.Button(self, label='Draw')
        self.Bind(wx.EVT_BUTTON, self.onDrawRates, self.rates_draw_button)

        # export
        self.sort_ctrl = wx.CheckBox(self, label="Sort by rates")
        self.sort_ctrl.SetValue(True)
        self.spin_precision_rates = wx.SpinCtrl(self, size=(60, -1),
                                               min=2, max=16, initial=6)
        self.export_rates_button = wx.Button(self, label='Export')
        self.Bind(wx.EVT_BUTTON, self.onExportRates, self.export_rates_button)

        # rate profiles
        self.reaction_choices = MultiChoices(self, self.ctout.reactions,
                                             wx.VERTICAL, add_num=True,
                                             top=self.getTopRate)
        self.shorten_legends_ctrl = wx.CheckBox(self, label="Shorten labels")
        self.shorten_legends_ctrl.SetValue(False)
        self.profile_plot_button = wx.Button(self, label='Plot')
        self.Bind(wx.EVT_BUTTON, self.onPlotProfile, self.profile_plot_button)
        
        # export profiles
        self.choice_header = wx.Choice(self, choices=self.header_list)
        self.choice_header.SetSelection(self.header_default_index)
        self.spin_precision_profile = wx.SpinCtrl(self, size=(60, -1),
                                                  min=2, max=16, initial=6)
        self.export_profile_button = wx.Button(self, label='Export')
        self.Bind(wx.EVT_BUTTON, self.onExportProfile, self.export_profile_button)
        
        # Layout
        sp_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Species && rate unit'),
                                    wx.VERTICAL)
        sp_vbox.Add(self.species_choice, 0, border=0)
        sp_vbox.Add(self.rate_unit_choices, 0, border=0)
        
        rates_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Production rates'),
                                       wx.VERTICAL)
        rates_vbox.Add(self.xvar_slider, 0, wx.ALL, 4)
        rates_vbox.Add(wx.StaticText(self, label='Max. num.:'), 0, wx.ALL, 4)
        rates_vbox.Add(self.maxnum_spin, 0, wx.ALL, 4)
        rates_vbox.Add(self.shorten_labels_ctrl, 0, wx.ALL, 4)
        rates_vbox.Add(self.rates_draw_button, 0, wx.ALL, 4)
        rates_vbox.AddSpacer(8)
        rates_vbox.Add(wx.StaticText(self, label='Precision:'), 0, wx.ALL, 4)
        rates_vbox.Add(self.spin_precision_rates, 0, wx.ALL, 4)
        rates_vbox.Add(self.sort_ctrl, 0, wx.ALL, 4)
        rates_vbox.Add(self.export_rates_button, 0, wx.ALL, 4)
        
        sp_rates_vbox = wx.BoxSizer(wx.VERTICAL)
        sp_rates_vbox.Add(sp_vbox, 0, border=0)
        sp_rates_vbox.AddSpacer(10)
        sp_rates_vbox.Add(rates_vbox, 0, border=0)
        
        prof_vbox = wx.BoxSizer(wx.VERTICAL)
        prof_vbox.Add(self.shorten_legends_ctrl, 0, wx.ALL, 4)
        prof_vbox.Add(self.profile_plot_button, 0, wx.ALL, 4)
        prof_vbox.AddSpacer(12)
        prof_vbox.Add(wx.StaticText(self, label='Header:'), 0, wx.ALL, 4)
        prof_vbox.Add(self.choice_header, 0, wx.ALL, 4)
        prof_vbox.AddSpacer(2)
        prof_vbox.Add(wx.StaticText(self, label='Precision:'), 0, wx.ALL, 4)
        prof_vbox.Add(self.spin_precision_profile, 0, wx.ALL, 4)
        prof_vbox.AddSpacer(2)
        prof_vbox.Add(self.export_profile_button, 0, wx.ALL, 4)
        
        prof_hbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Profile'), wx.HORIZONTAL)
        prof_hbox.Add(self.reaction_choices, 0, border=0)
        prof_hbox.Add(prof_vbox, 0, border=0)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(sp_rates_vbox, 0, wx.ALL, 4)
        hbox.Add(prof_hbox, 0, wx.ALL, 4)
        
        self.SetSizer(hbox)

    def onDrawRates(self, event): self.drawRates()
    def onPlotProfile(self, event): self.plotProfile()
    def onExportRates(self, event): self.exportRates()
    def onExportProfile(self, event): self.exportProfile()
        
    def drawRates(self):
        name = self.species_choice.getCurrentString()
        maxnum = self.maxnum_spin.GetValue()
        index = self.xvar_slider.getIndex()
        xvar = self.xvar_slider.getXvar()
        rate_unit = self.rate_unit_choices.getString('rate')
        shorten = self.shorten_labels_ctrl.IsChecked()
        
        self.message('Reading production rates for %s...' % (name,))
        totrate, rate_list = self.ctout.getProductionRatesOfAt(name, index, rate_unit)
        d = {}
        for i,x in enumerate(rate_list): d[i] = x
        items = sorted(d.items(), key=lambda x: abs(x[1]), reverse=True)[:maxnum]
        rinds, rates = list(zip(*items))

        self.message('Drawing production rates for %s...' % (name,))
        
        frame = PlotFrame(self, title='Top %d production rates for %s at %s = %s' \
                          % (maxnum, name, self.ctout.xvarname, xvar))
        ax = frame.addSubplot(111,
                              xlabel='Production rates for %s / %s' % (name, rate_unit))
        
        pos = list(reversed(range(maxnum)))
        ax.barh(pos[0]+1, [totrate], 0.8, color='#bb8888', label='Total')
        ax.barh(pos, rates, 0.8, color='#bbbbbb', label='Production rates')
        if shorten:
            ax.text(0.0, pos[0]+1.4, 'Total', ha='center', va='center')
            for p,i in zip(pos, rinds):
                ax.text(0.0, p+0.4, '%s' % (self.ctout.reactions[i]),
                        ha='center', va='center')
        else:
            ax.text(0.0, pos[0]+1.4, 'Total', ha='center', va='center')
            for p,i in zip(pos, rinds):
                ax.text(0.0, p+0.4, '%d: %s' % (i+1, self.ctout.reactions[i]),
                        ha='center', va='center')
        ax.axvline(x=0, color='#444444')
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylim(-0.2, pos[0]+2.0)
        xmin, xmax = ax.get_xlim()
        xabsmax = max(abs(xmin), abs(xmax))
        ax.set_xlim(-xabsmax, xabsmax)
        frame.draw()
        frame.Show()
        self.message('Drawing production rates for %s... done.' % (name,))
        
    def plotProfile(self):
        name = self.species_choice.getCurrentString()
        rate_unit = self.rate_unit_choices.getString('rate')
        indices = self.reaction_choices.getCheckedIndices()
        reactions = self.reaction_choices.getCheckedStrings()
        shorten = self.shorten_legends_ctrl.IsChecked()
        if len(indices) == 0:
            wx.MessageBox('Select reactions', 'alert', 
                          wx.OK | wx.CENTRE | wx.ICON_EXCLAMATION)
            return
        if shorten:
            labels = ['R%d' % (i+1,) for i in indices]
        else:
            labels = ['%d: %s' % (i+1, x) for i,x in zip(indices, reactions)]
        
        self.message('Reading production rates profiles for %s...' % (name,))
        xvars = np.asarray(self.ctout.getXvarList())
        totrates, cntrb_lists = \
                  self.ctout.getProductionRatesProfileOf(name, indices, rate_unit)
        totrates = np.asarray(totrates)
        cntrb_lists = list(map(np.asarray, cntrb_lists))
        
        self.message('Ploting production rates for %s...' % (name,))
        frame = PlotFrame(self, title='Production rates for %s' % (name,))
        frame.setLegend(1)
        ax = frame.addSubplot(111, xlabel=self.ctout.xvarname,
                              ylabel='Production rates for %s / %s' % 
                              (name, rate_unit))
        
        clrsty = genColorAndStyles()
        clr, sty = next(clrsty)
        ax.plot(xvars, totrates, color=clr, ls=sty, lw=1.5, label='Total')
        for s,c in zip(labels, cntrb_lists):
            clr, sty = next(clrsty)
            ax.plot(xvars, c, color=clr, ls=sty, lw=1.5, label=s)
        frame.draw()
        frame.Show()
        self.message('Ploting production rates profiles for %s... done' % (name,))

    def exportRates(self):
        name = self.species_choice.getCurrentString()
        index = self.xvar_slider.getIndex()
        inidir = os.path.dirname(self.ctout.datfile)
        dlg = wx.FileDialog(self, message='Export production rates for %s' % (name,),
                            defaultDir=inidir, 
                            defaultFile='prod_%s_%d.txt' % (name, index),
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            #dlg.Destroy()
        else:
            #dlg.Destroy()
            return
        
        rate_unit = self.rate_unit_choices.getString('rate')
        precision = self.spin_precision_rates.GetValue()
        self.message('Exporting production rates for %s...' % (name,))
        wx.BeginBusyCursor()
        fp = open(path, 'w')
        self.ctout.exportProductionRatesOfAt(fp, name, index, rate_unit=rate_unit,
                                            sort=self.sort_ctrl.IsChecked(),
                                            precision=precision)
        fp.close()
        wx.EndBusyCursor()
        self.message('Exporting production rates for %s... done' % (name,))

    def exportProfile(self):
        name = self.species_choice.getCurrentString()
        inidir = os.path.dirname(self.ctout.datfile)
        dlg = wx.FileDialog(self, 
                            message='Export production rates profiles for %s' % (name,),
                            defaultDir=inidir, 
                            defaultFile='prod_profile_%s' % (name,),
                            wildcard=self.wildcard,
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            #dlg.Destroy()
        else:
            #dlg.Destroy()
            return
        rate_unit = self.rate_unit_choices.getString('rate')
        indices = self.reaction_choices.getCheckedIndices()
        sep = self.sep_values[dlg.GetFilterIndex()]
        header = self.header_values[self.choice_header.GetCurrentSelection()]
        precision = self.spin_precision_profile.GetValue()
        
        self.message('Exporting production rates profiles for %s...' % (name,))
        wx.BeginBusyCursor()
        fp = open(path, 'w')
        self.ctout.exportProductionRatesProfileOf(fp, name, indices,
                                                 rate_unit=rate_unit,
                                                 sep=sep, header=header,
                                                 precision=precision)
        fp.close()
        wx.EndBusyCursor()
        self.message('Exporting production rates profiles for %s... done.' % (name,))

    def getTopRate(self, num):
        name = self.species_choice.getCurrentString()
        index = self.xvar_slider.getIndex()
        rate_unit = self.rate_unit_choices.getString('rate')
        totrate, rate_list = self.ctout.getProductionRatesOfAt(name, index, rate_unit)
        tmp = []
        for i, x in enumerate(rate_list): tmp.append((abs(x), i))
        tmp = sorted(tmp, reverse=True)
        num = min(num, len(tmp))
        return [x[1] for x in tmp[:num]]

class HeatPanel(ScrolledPanel):

    header_list = ('None', 'Add', 'With "#"')
    header_values = (None, '', '#')
    header_default_index = 2
    wildcard = 'Space separated files (*.dat)|*.dat|' \
               'Comma separated files (*.csv)|*.csv|' \
               'Tab separated files (*.tsv)|*.tsv|' \
               'All files (space separated) (*.*)|*.*'
    sep_values = (' ', ', ', '\t', ' ')
    
    def __init__(self, parent, mainframe):
        ScrolledPanel.__init__(self, parent)
        self.ctout = mainframe.ctout
        self.message = mainframe.message
        self.create()
        self.SetupScrolling()
        
    def create(self):
        self.heat_choice = SpeciesChoice(self, species=["u", "h"],
                                         add_num=True, init=0)

        # rate at selected xvar
        self.xvar_slider = XvarSlider(self, self.ctout)
        self.maxnum_spin = wx.SpinCtrl(self, size=(60, -1),
                                       min=min(self.ctout.nr, 5),
                                       max=self.ctout.nr,
                                       initial=min(self.ctout.nr, 10))
        self.shorten_labels_ctrl = wx.CheckBox(self, label="Shorten labels")
        self.shorten_labels_ctrl.SetValue(False)
        self.rates_draw_button = wx.Button(self, label='Draw')
        self.Bind(wx.EVT_BUTTON, self.onDrawRates, self.rates_draw_button)

        # export
        self.sort_ctrl = wx.CheckBox(self, label="Sort by rates")
        self.sort_ctrl.SetValue(True)
        self.spin_precision_rates = wx.SpinCtrl(self, size=(60, -1),
                                               min=2, max=16, initial=6)
        self.export_rates_button = wx.Button(self, label='Export')
        self.Bind(wx.EVT_BUTTON, self.onExportRates, self.export_rates_button)

        # rate profiles
        self.reaction_choices = MultiChoices(self, self.ctout.reactions,
                                             wx.VERTICAL, add_num=True,
                                             top=self.getTopRate)
        self.shorten_legends_ctrl = wx.CheckBox(self, label="Shorten labels")
        self.shorten_legends_ctrl.SetValue(False)
        self.profile_plot_button = wx.Button(self, label='Plot')
        self.Bind(wx.EVT_BUTTON, self.onPlotProfile, self.profile_plot_button)
        
        # export profiles
        self.choice_header = wx.Choice(self, choices=self.header_list)
        self.choice_header.SetSelection(self.header_default_index)
        self.spin_precision_profile = wx.SpinCtrl(self, size=(60, -1),
                                                  min=2, max=16, initial=6)
        self.export_profile_button = wx.Button(self, label='Export')
        self.Bind(wx.EVT_BUTTON, self.onExportProfile, self.export_profile_button)
        
        # Layout
        sp_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='heat'),
                                    wx.VERTICAL)
        sp_vbox.Add(self.heat_choice, 0, border=0)
        
        rates_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Production rates'),
                                       wx.VERTICAL)
        rates_vbox.Add(self.xvar_slider, 0, wx.ALL, 4)
        rates_vbox.Add(wx.StaticText(self, label='Max. num.:'), 0, wx.ALL, 4)
        rates_vbox.Add(self.maxnum_spin, 0, wx.ALL, 4)
        rates_vbox.Add(self.shorten_labels_ctrl, 0, wx.ALL, 4)
        rates_vbox.Add(self.rates_draw_button, 0, wx.ALL, 4)
        rates_vbox.AddSpacer(8)
        rates_vbox.Add(wx.StaticText(self, label='Precision:'), 0, wx.ALL, 4)
        rates_vbox.Add(self.spin_precision_rates, 0, wx.ALL, 4)
        rates_vbox.Add(self.sort_ctrl, 0, wx.ALL, 4)
        rates_vbox.Add(self.export_rates_button, 0, wx.ALL, 4)
        
        sp_rates_vbox = wx.BoxSizer(wx.VERTICAL)
        sp_rates_vbox.Add(sp_vbox, 0, border=0)
        sp_rates_vbox.AddSpacer(10)
        sp_rates_vbox.Add(rates_vbox, 0, border=0)
        
        prof_vbox = wx.BoxSizer(wx.VERTICAL)
        prof_vbox.Add(self.shorten_legends_ctrl, 0, wx.ALL, 4)
        prof_vbox.Add(self.profile_plot_button, 0, wx.ALL, 4)
        prof_vbox.AddSpacer(12)
        prof_vbox.Add(wx.StaticText(self, label='Header:'), 0, wx.ALL, 4)
        prof_vbox.Add(self.choice_header, 0, wx.ALL, 4)
        prof_vbox.AddSpacer(2)
        prof_vbox.Add(wx.StaticText(self, label='Precision:'), 0, wx.ALL, 4)
        prof_vbox.Add(self.spin_precision_profile, 0, wx.ALL, 4)
        prof_vbox.AddSpacer(2)
        prof_vbox.Add(self.export_profile_button, 0, wx.ALL, 4)
        
        prof_hbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Profile'), wx.HORIZONTAL)
        prof_hbox.Add(self.reaction_choices, 0, border=0)
        prof_hbox.Add(prof_vbox, 0, border=0)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(sp_rates_vbox, 0, wx.ALL, 4)
        hbox.Add(prof_hbox, 0, wx.ALL, 4)
        
        self.SetSizer(hbox)
    def onDrawRates(self, event): self.drawRates()
    def onPlotProfile(self, event): self.plotProfile()
    def onExportRates(self, event): self.exportRates()
    def onExportProfile(self, event): self.exportProfile()
        
    def drawRates(self):
        name = self.heat_choice.getCurrentString()
        maxnum = self.maxnum_spin.GetValue()
        index = self.xvar_slider.getIndex()
        xvar = self.xvar_slider.getXvar()
        shorten = self.shorten_labels_ctrl.IsChecked()
        
        self.message('Reading heat (%s) production rates...' % (name,))
        totrate, rate_list = self.ctout.getHeatProductionRatesOfAt(name, index)
        d = {}
        for i,x in enumerate(rate_list): d[i] = x
        items = sorted(d.items(), key=lambda x: abs(x[1]), reverse=True)[:maxnum]
        rinds, rates = list(zip(*items))

        self.message('Drawing heat (%s) production rates...' % (name,))
        
        frame = PlotFrame(self, title='Top %d production rates for %s at %s = %s' \
                          % (maxnum, name, self.ctout.xvarname, xvar))
        ax = frame.addSubplot(111,
                              xlabel='Production rates for %s / J/cm3-s' % name)
        
        pos = list(reversed(range(maxnum)))
        ax.barh(pos[0]+1, [totrate], 0.8, color='#bb8888', label='Total')
        ax.barh(pos, rates, 0.8, color='#bbbbbb', label='Production rates')
        if shorten:
            ax.text(0.0, pos[0]+1.4, 'Total', ha='center', va='center')
            for p,i in zip(pos, rinds):
                ax.text(0.0, p+0.4, '%s' % (self.ctout.reactions[i]),
                        ha='center', va='center')
        else:
            ax.text(0.0, pos[0]+1.4, 'Total', ha='center', va='center')
            for p,i in zip(pos, rinds):
                ax.text(0.0, p+0.4, '%d: %s' % (i+1, self.ctout.reactions[i]),
                        ha='center', va='center')
        ax.axvline(x=0, color='#444444')
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylim(-0.2, pos[0]+2.0)
        xmin, xmax = ax.get_xlim()
        xabsmax = max(abs(xmin), abs(xmax))
        ax.set_xlim(-xabsmax, xabsmax)
        frame.draw()
        frame.Show()
        self.message('Drawing heat (%s) production rates done.' % (name,))
        
    def plotProfile(self):
        name = self.heat_choice.getCurrentString()
        indices = self.reaction_choices.getCheckedIndices()
        reactions = self.reaction_choices.getCheckedStrings()
        shorten = self.shorten_legends_ctrl.IsChecked()
        if len(indices) == 0:
            wx.MessageBox('Select reactions', 'alert', 
                          wx.OK | wx.CENTRE | wx.ICON_EXCLAMATION)
            return
        if shorten:
            labels = ['R%d' % (i+1,) for i in indices]
        else:
            labels = ['%d: %s' % (i+1, x) for i,x in zip(indices, reactions)]
        
        self.message('Reading heat (%s) production rates profiles...' % (name,))
        xvars = np.asarray(self.ctout.getXvarList())
        totrates, cntrb_lists = \
                  self.ctout.getHeatProductionRatesProfileOf(name, indices)
        totrates = np.asarray(totrates)
        cntrb_lists = list(map(np.asarray, cntrb_lists))
        
        self.message('Ploting heat (%s) production rates...' % (name,))
        frame = PlotFrame(self, title='Production rates for %s' % (name,))
        frame.setLegend(1)
        ax = frame.addSubplot(111, xlabel=self.ctout.xvarname,
                              ylabel='Heat (%s) production rates / J/cm3-s' % 
                              name)
        
        clrsty = genColorAndStyles()
        clr, sty = next(clrsty)
        ax.plot(xvars, totrates, color=clr, ls=sty, lw=1.5, label='Total')
        for s,c in zip(labels, cntrb_lists):
            clr, sty = next(clrsty)
            ax.plot(xvars, c, color=clr, ls=sty, lw=1.5, label=s)
        frame.draw()
        frame.Show()
        self.message('Ploting heat (%s) production rates profiles done' % (name,))

    def exportRates(self):
        name = self.heat_choice.getCurrentString()
        index = self.xvar_slider.getIndex()
        inidir = os.path.dirname(self.ctout.datfile)
        dlg = wx.FileDialog(self, message='Export heat (%s) production rates' % (name,),
                            defaultDir=inidir, 
                            defaultFile='prod_%s_%d.txt' % (name, index),
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            #dlg.Destroy()
        else:
            #dlg.Destroy()
            return
        
        precision = self.spin_precision_rates.GetValue()
        self.message('Exporting heat (%s) production rates for...' % (name,))
        wx.BeginBusyCursor()
        fp = open(path, 'w')
        self.ctout.exportHeatProductionRatesOfAt(fp, name, index,
                                            sort=self.sort_ctrl.IsChecked(),
                                            precision=precision)
        fp.close()
        wx.EndBusyCursor()
        self.message('Exporting heat (%s) production rates done' % (name,))

    def exportProfile(self):
        name = self.heat_choice.getCurrentString()
        inidir = os.path.dirname(self.ctout.datfile)
        dlg = wx.FileDialog(self, 
                            message='Export heat (%s) production rates profiles' % (name,),
                            defaultDir=inidir, 
                            defaultFile='prod_profile_%s' % (name,),
                            wildcard=self.wildcard,
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            #dlg.Destroy()
        else:
            #dlg.Destroy()
            return

        indices = self.reaction_choices.getCheckedIndices()
        sep = self.sep_values[dlg.GetFilterIndex()]
        header = self.header_values[self.choice_header.GetCurrentSelection()]
        precision = self.spin_precision_profile.GetValue()
        
        self.message('Exporting heat (%s) production rates profiles' % (name,))
        wx.BeginBusyCursor()
        fp = open(path, 'w')
        self.ctout.exportHeatProductionRatesProfileOf(fp, name, indices,
                                                 sep=sep, header=header,
                                                 precision=precision)
        fp.close()
        wx.EndBusyCursor()
        self.message('Exporting heat (%s) production rates profiles for done.' % (name,))

    def getTopRate(self, num):
        name = self.heat_choice.getCurrentString()
        index = self.xvar_slider.getIndex()
        totrate, rate_list = self.ctout.getHeatProductionRatesOfAt(name, index)
        tmp = []
        for i, x in enumerate(rate_list): tmp.append((abs(x), i))
        tmp = sorted(tmp, reverse=True)
        num = min(num, len(tmp))
        return [x[1] for x in tmp[:num]]




class SensPanel(ScrolledPanel):
    
    header_list = ('None', 'Add', 'With "#"')
    header_values = (None, '', '#')
    header_default_index = 2
    wildcard = 'Space separated files (*.dat)|*.dat|' \
               'Comma separated files (*.csv)|*.csv|' \
               'Tab separated files (*.tsv)|*.tsv|' \
               'All files (space separated) (*.*)|*.*'
    sep_values = (' ', ', ', '\t', ' ')
    
    def __init__(self, parent, mainframe):
        ScrolledPanel.__init__(self, parent)
        self.ctout = mainframe.ctout
        self.message = mainframe.message
        self.create()
        self.SetupScrolling()
        
    def create(self):
        # var
        self.variables_choice = SpeciesChoice(self, species=self.ctout.sens_vars,
                                              add_num=True, init=0)

        # sens at selected xvar
        self.xvar_slider = XvarSlider(self, self.ctout)
        self.maxnum_spin = wx.SpinCtrl(self, size=(60, -1),
                                       min=min(self.ctout.nr, 5),
                                       max=self.ctout.nr,
                                       initial=min(self.ctout.nr, 10))
        self.shorten_labels_ctrl = wx.CheckBox(self, label="Shorten labels")
        self.shorten_labels_ctrl.SetValue(False)
        self.sens_draw_button = wx.Button(self, label='Draw')
        self.Bind(wx.EVT_BUTTON, self.onDrawSens, self.sens_draw_button)

        # export
        self.sort_ctrl = wx.CheckBox(self, label="Sort by coeffs")
        self.sort_ctrl.SetValue(True)
        self.spin_precision_sens = wx.SpinCtrl(self, size=(60, -1),
                                               min=2, max=16, initial=6)
        self.export_sens_button = wx.Button(self, label='Export')
        self.Bind(wx.EVT_BUTTON, self.onExportSens, self.export_sens_button)

        # sens profiles
        self.reaction_choices = MultiChoices(self, self.ctout.reactions,
                                             wx.VERTICAL, add_num=True,
                                             top=self.getTopSens)
        self.shorten_legends_ctrl = wx.CheckBox(self, label="Shorten labels")
        self.shorten_legends_ctrl.SetValue(False)
        self.profile_plot_button = wx.Button(self, label='Plot')
        self.Bind(wx.EVT_BUTTON, self.onPlotProfile, self.profile_plot_button)
        
        # export profiles
        self.choice_header = wx.Choice(self, choices=self.header_list)
        self.choice_header.SetSelection(self.header_default_index)
        self.spin_precision_profile = wx.SpinCtrl(self, size=(60, -1),
                                                  min=2, max=16, initial=6)
        self.export_profile_button = wx.Button(self, label='Export')
        self.Bind(wx.EVT_BUTTON, self.onExportProfile, self.export_profile_button)
        
        # Layout
        var_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Variable'),
                                     wx.VERTICAL)
        var_vbox.Add(self.variables_choice, 0, border=0)
        
        sens_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Sensitivity'), wx.VERTICAL)
        sens_vbox.Add(self.xvar_slider, 0, wx.ALL, 4)
        sens_vbox.Add(wx.StaticText(self, label='Max. num.:'), 0, wx.ALL, 4)
        sens_vbox.Add(self.maxnum_spin, 0, wx.ALL, 4)
        sens_vbox.Add(self.shorten_labels_ctrl, 0, wx.ALL, 4)
        sens_vbox.Add(self.sens_draw_button, 0, wx.ALL, 4)
        sens_vbox.AddSpacer(8)
        sens_vbox.Add(wx.StaticText(self, label='Precision:'), 0, wx.ALL, 4)
        sens_vbox.Add(self.spin_precision_sens, 0, wx.ALL, 4)
        sens_vbox.Add(self.sort_ctrl, 0, wx.ALL, 4)
        sens_vbox.Add(self.export_sens_button, 0, wx.ALL, 4)
        
        var_sens_vbox = wx.BoxSizer(wx.VERTICAL)
        var_sens_vbox.Add(var_vbox, 0, border=0)
        var_sens_vbox.AddSpacer(10)
        var_sens_vbox.Add(sens_vbox, 0, border=0)
        
        prof_vbox = wx.BoxSizer(wx.VERTICAL)
        prof_vbox.Add(self.shorten_legends_ctrl, 0, wx.ALL, 4)
        prof_vbox.Add(self.profile_plot_button, 0, wx.ALL, 4)
        prof_vbox.AddSpacer(12)
        prof_vbox.Add(wx.StaticText(self, label='Header:'), 0, wx.ALL, 4)
        prof_vbox.Add(self.choice_header, 0, wx.ALL, 4)
        prof_vbox.AddSpacer(2)
        prof_vbox.Add(wx.StaticText(self, label='Precision:'), 0, wx.ALL, 4)
        prof_vbox.Add(self.spin_precision_profile, 0, wx.ALL, 4)
        prof_vbox.AddSpacer(2)
        prof_vbox.Add(self.export_profile_button, 0, wx.ALL, 4)
        
        prof_hbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Profile'), wx.HORIZONTAL)
        prof_hbox.Add(self.reaction_choices, 0, border=0)
        prof_hbox.Add(prof_vbox, 0, border=0)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(var_sens_vbox, 0, wx.ALL, 4)
        hbox.Add(prof_hbox, 0, wx.ALL, 4)
        
        self.SetSizer(hbox)

    def onDrawSens(self, event): self.drawSens()
    def onPlotProfile(self, event): self.plotProfile()
    def onExportSens(self, event): self.exportSens()
    def onExportProfile(self, event): self.exportProfile()
        
    def drawSens(self):
        name = self.variables_choice.getCurrentString()
        ivar = self.variables_choice.getCurrentIndex()
        maxnum = self.maxnum_spin.GetValue()
        index = self.xvar_slider.getIndex()
        xvar = self.xvar_slider.getXvar()
        shorten = self.shorten_labels_ctrl.IsChecked()
        
        self.message('Reading sensitivity for %s...' % (name,))
        sens_list = self.ctout.getSensCoeffsOfAt(ivar, index)
        d = {}
        for i,x in enumerate(sens_list): d[i] = x
        items = sorted(d.items(), key=lambda x: abs(x[1]), reverse=True)[:maxnum]
        rinds, sens = list(zip(*items))

        self.message('Drawing sensitivity for %s...' % (name,))
        frame = PlotFrame(self, title='Top %d sensitivity for %s at %s = %s' 
                          % (maxnum, name, self.ctout.xvarname, xvar))
        ax = frame.addSubplot(111, xlabel='Sensitivity for %s' % (name,))
        
        pos = list(reversed(range(maxnum)))
        ax.barh(pos, sens, 0.8, color='#bbbbbb', label='Sensitivity')
        if shorten:
            for p,i in zip(pos, rinds):
                ax.text(0.0, p+0.4, '%s' % (self.ctout.reactions[i]),
                        ha='center', va='center')
        else:
            for p,i in zip(pos, rinds):
                ax.text(0.0, p+0.4, '%d: %s' % (i+1, self.ctout.reactions[i]),
                        ha='center', va='center')
        ax.axvline(x=0, color='#444444')
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylim(-0.2, pos[0]+1)
        xmin, xmax = ax.get_xlim()
        xabsmax = max(abs(xmin), abs(xmax))
        ax.set_xlim(-xabsmax, xabsmax)
        frame.draw()
        frame.Show()
        self.message('Drawing sensitivity for %s... done.' % (name,))
        
    def plotProfile(self):
        name = self.variables_choice.getCurrentString()
        ivar = self.variables_choice.getCurrentIndex()
        indices = self.reaction_choices.getCheckedIndices()
        reactions = self.reaction_choices.getCheckedStrings()
        shorten = self.shorten_legends_ctrl.IsChecked()
        if len(indices) == 0:
            wx.MessageBox('Select reactions', 'alert', 
                          wx.OK | wx.CENTRE | wx.ICON_EXCLAMATION)
            return
        if shorten:
            labels = ['R%d' % (i+1,) for i in indices]
        else:
            labels = ['%d: %s' % (i+1, x) for i,x in zip(indices, reactions)]
        
        self.message('Reading sensitivity profiles for %s...' % (name,))
        xvars = np.asarray(self.ctout.getXvarList())
        senses = list(map(np.asarray,
                     self.ctout.getSensCoeffsProfileOf(ivar, indices)))
        
        self.message('Ploting sensitivity profiles for %s...' % (name,))
        frame = PlotFrame(self, title='Sensitivity profiles for %s' % (name,))
        frame.setLegend(1)
        ax = frame.addSubplot(111, xlabel=self.ctout.xvarname,
                              ylabel='Sensitivity for %s' % (name,))
        
        clrsty = genColorAndStyles()
        for s,c in zip(labels, senses):
            clr, sty = next(clrsty)
            ax.plot(xvars, c, color=clr, ls=sty, lw=1.5, label=s)
        frame.draw()
        frame.Show()
        self.message('Ploting sensitivity profiles for %s... done' % (name,))

    def exportSens(self):
        name = self.variables_choice.getCurrentString()
        index = self.xvar_slider.getIndex()
        inidir = os.path.dirname(self.ctout.datfile)
        dlg = wx.FileDialog(self, message='Export sensitivity for %s' % (name,),
                            defaultDir=inidir, 
                            defaultFile='sens_%s_%d.txt' % (name, index),
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            #dlg.Destroy()
        else:
            #dlg.Destroy()
            return
        ivar = self.variables_choice.getCurrentIndex()
        
        self.message('Exporting sensitivity for %s...' % (name,))
        wx.BeginBusyCursor()
        fp = open(path, 'w')
        self.ctout.exportSensCoeffsOfAt(fp, ivar, index,
                                       sort=self.sort_ctrl.IsChecked(),
                                       precision=self.spin_precision_sens.GetValue())
        fp.close()
        wx.EndBusyCursor()
        self.message('Exporting sensitivity for %s... done' % (name,))

    def exportProfile(self):
        name = self.variables_choice.getCurrentString()
        inidir = os.path.dirname(self.ctout.datfile)
        dlg = wx.FileDialog(self, message='Export sensitivity profiles for %s' % (name,),
                            defaultDir=inidir, 
                            defaultFile='sens_profile_%s' % (name,),
                            wildcard=self.wildcard,
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            #dlg.Destroy()
        else:
            #dlg.Destroy()
            return
        ivar = self.variables_choice.getCurrentIndex()
        indices = self.reaction_choices.getCheckedIndices()
        sep = self.sep_values[dlg.GetFilterIndex()]
        header = self.header_values[self.choice_header.GetCurrentSelection()]
        precision = self.spin_precision_profile.GetValue()
        
        self.message('Exporting sensitivity profiles for %s...' % (name,))
        wx.BeginBusyCursor()
        fp = open(path, 'w')
        self.ctout.exportSensCoeffsProfileOf(fp, ivar, indices,
                                            sep=sep, header=header,
                                            precision=precision)
        fp.close()
        wx.EndBusyCursor()
        self.message('Exporting sensitivity profiles for %s... done.' % (name,))

    def getTopSens(self, num):
        ivar = self.variables_choice.getCurrentIndex()
        index = self.xvar_slider.getIndex()
        sens_list = self.ctout.getSensCoeffsOfAt(ivar, index)
        tmp = []
        for i, x in enumerate(sens_list): tmp.append((abs(x), i))
        tmp = sorted(tmp, reverse=True)
        num = min(num, len(tmp))
        return [x[1] for x in tmp[:num]]


class PathFluxPanel(ScrolledPanel):
    
    default_rtol = 0.01
    warn_species_num = 200
    warn_connex_num = 800
    sdirections = ["Both", "Consumption", "Formation", "None"]
    layouts = ["Grid", "Random", "dot"]
    layouts_tooltips = ["Grid layout",
                        "Random layout",
                        "Hierarchical layout (Graphviz-dot)"]
    wildcard = 'Text files (*.txt)|*.txt|' \
               'DOT formatted files (*.dot)|*.dot|' \
               'All files (space separated) (*.*)|*.*'
    formats = ('txt', 'dot', 'txt')
    
    def __init__(self, parent, mainframe):
        ScrolledPanel.__init__(self, parent)
        self.ctout = mainframe.ctout
        self.message = mainframe.message
        self.create()
        self.SetupScrolling()

    def create(self):
        # species
        self.include_choices = MultiChoices(self, self.ctout.species, wx.VERTICAL,
                                            wx.HORIZONTAL, staticbox=True,
                                            staticbox_label='Include',
                                            add_num=True, top=self.getTopFrac)
        self.exclude_choices = MultiChoices(self, self.ctout.species, wx.VERTICAL,
                                            wx.HORIZONTAL, staticbox=True,
                                            staticbox_label='Exclude',
                                            add_num=True, small=self.getSmall)
        # configs
        self.fluxtypes = ['Mass'] + self.ctout.elements
        self.fluxtype_choice = wx.Choice(self, choices=self.fluxtypes)
        self.fluxtype_choice.SetSelection(0)
        self.xvar_slider = XvarSlider(self, self.ctout)
        self.prior_included_ctrl = wx.CheckBox(self, label="Included-prior")
        self.prior_included_ctrl.SetValue(False)
        self.remove_unconnected_ctrl = wx.CheckBox(self, label="Remove uncon.")
        self.remove_unconnected_ctrl.SetValue(True)
        self.rtol_ctrl = wx.TextCtrl(self, value=str(self.default_rtol))
        self.Bind(wx.EVT_TEXT, self.onRtolChange, self.rtol_ctrl)
        self.search_direction = wx.RadioBox(self, label="Search", choices=self.sdirections,
                                            style=wx.RA_VERTICAL)
        
        # draw
        self.draw_layout = wx.RadioBox(self, label="layout", choices=self.layouts,
                                       style=wx.RA_VERTICAL)
        for i, s in enumerate(self.layouts_tooltips):
            self.draw_layout.SetItemToolTip(i, s)
        self.draw_button = wx.Button(self, label='Draw')
        self.Bind(wx.EVT_BUTTON, self.onDraw, self.draw_button)
        self.export_button = wx.Button(self, label='Export')
        self.Bind(wx.EVT_BUTTON, self.onExport, self.export_button)
        
        # Layout
        config_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Config'), wx.VERTICAL)
        config_vbox.Add(wx.StaticText(self, label='Flux type:'), 0, wx.ALL, 4)
        config_vbox.Add(self.fluxtype_choice, 0, wx.ALL, 4)
        config_vbox.Add(self.xvar_slider, 0, wx.ALL, 4)
        config_vbox.Add(self.prior_included_ctrl, 0, wx.ALL, 4)
        config_vbox.Add(self.remove_unconnected_ctrl, 0, wx.ALL, 4)
        config_vbox.Add(wx.StaticText(self, label='rtol:'), 0, wx.ALL, 4)
        config_vbox.Add(self.rtol_ctrl, 0, wx.ALL, 4)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(config_vbox, 0, border=0)
        vbox.Add(self.search_direction, 0, wx.ALL, 4)
        
        draw_vbox = wx.StaticBoxSizer(wx.StaticBox(self, label='Draw'), wx.VERTICAL)
        draw_vbox.Add(self.draw_button, 0, wx.ALL, 4)
        draw_vbox.Add(self.export_button, 0, wx.ALL, 4)

        vbox2 = wx.BoxSizer(wx.VERTICAL)
        vbox2.Add(self.draw_layout, 0, wx.ALL, 4)
        vbox2.Add(draw_vbox, 0, border=0)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.include_choices, 0, wx.ALL, 4)
        hbox.Add(self.exclude_choices, 0, wx.ALL, 4)
        hbox.Add(vbox, 0, wx.ALL, 4)
        hbox.Add(vbox2, 0, wx.ALL, 4)
        
        self.SetSizer(hbox)
        
    def onDraw(self, event): self.draw()
    def onExport(self, event): self.export()
    def onRtolChange(self, event): self.validate()
        
    def validate(self):
        try:
            rtol = float(self.rtol_ctrl.GetValue())
            if rtol < 0:
                self.draw_button.Enable(False)
                self.export_button.Enable(False)
            else:
                self.draw_button.Enable(True)
                self.export_button.Enable(True)
        except:
            self.draw_button.Enable(False)
            self.export_button.Enable(False)
            
    def draw(self):
        fracs, mat, spdict, connex, concontrib = self.getConnex()
        name = self.fluxtypes[self.fluxtype_choice.GetCurrentSelection()]
        layout = self.draw_layout.GetStringSelection()
        index = self.xvar_slider.getIndex()
        xvar = self.xvar_slider.getXvar()
        if name == "Mass": 
            str1 = "mass"
            str2 = "Mass Flux"
            unit = "g/cm3-s"
        else:
            str1 = "element"
            str2 = "Element Flux for %s" % (name,)
            unit = "atoms/cm3-s"
        if len(spdict) > self.warn_species_num or len(connex) > self.warn_connex_num:
            msg = "Visualizing a large network " + \
                  " (%d species, %d connections)\n" % (len(spdict), len(connex))+ \
                  "may take a long time"
            dlg = wx.MessageDialog(self, msg, "Warning",
                                   wx.OK | wx.CANCEL | wx.ICON_EXCLAMATION)
            result = dlg.ShowModal()
            #dlg.Destroy()
            if result != wx.ID_OK: return
        
        self.message('Drawing net %s flux (%d species, %d connections)...'
                     % (str1, len(spdict), len(connex)))
        frame = FluxFrame(self, title='Net %s at %s = %s' % (str2, self.ctout.xvarname, xvar))
        frame.setMisc(os.path.dirname(self.ctout.datfile), 'flux_%s_%d' % (name, index))
        frame.setFlux(spdict, self.ctout.species, fracs, connex, concontrib,
                      self.ctout.reactions, unit=unit)
        frame.draw(layout)
        frame.Show(True)
        self.message('Drawing net %s flux (%d species, %d connections)... done.'
                     % (str1, len(spdict), len(connex)))
        
    def export(self):
        fracs, mat, spdict, connex, concontrib = self.getConnex()
        name = self.fluxtypes[self.fluxtype_choice.GetCurrentSelection()]
        index = self.xvar_slider.getIndex()
        inidir = os.path.dirname(self.ctout.datfile)
        if name == "Mass": 
            str1 = "mass"
            str2 = "Mass Flux"
            unit = "g/cm3-s"
        else:
            str1 = "element"
            str2 = "Element Flux for %s" % (name,)
            unit = "atoms/cm3-s"
        dlg = wx.FileDialog(self, message='Export %s' % (str2,),
                            defaultDir=inidir,
                            defaultFile='flux_%s_%d' % (name, index),
                            wildcard=self.wildcard,
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            #dlg.Destroy()
        else:
            #dlg.Destroy()
            return
        self.message('Exporting net %s flux...' % (str1,))
        fp = open(path, 'w')
        format = self.formats[dlg.GetFilterIndex()]
        if format == "txt":
            xvar = self.xvar_slider.getXvar()
            if name == "Mass": header = "Net mass flux at %s = %s" % (self.ctout.xvarname, xvar)
            else: header = "Net element flux for %s at %s = %s" % (name, self.ctout.xvarname, xvar)
            self.ctout.dumpFluxToTxt(fp, fracs, spdict, connex, concontrib,
                                    header=header, unit=unit)
        elif format == "dot":
            self.ctout.dumpFluxToDot(fp, fracs, spdict, connex)
        fp.close()
        self.message('Exporting net %s flux... done.' % (str1,))

    def getTopFrac(self, num):
        name = self.fluxtypes[self.fluxtype_choice.GetCurrentSelection()]
        index = self.xvar_slider.getIndex()
        if name == "Mass": fracs = self.ctout.getMassFractionAt(index)
        else: fracs = self.ctout.getElementFractionAt(index, name)
        tmp = []
        for i, x in enumerate(fracs): tmp.append((x, i))
        tmp = sorted(tmp, reverse=True)
        num = min(num, len(tmp))
        return [x[1] for x in tmp[:num]]

    def getSmall(self):
        if "c" in self.ctout.elements:
            ielem = self.ctout.elements.index("c")
        elif "C" in self.ctout.elements:
            ielem = self.ctout.elements.index("C")
        else:
            return []
        smalls = []
        compositions = self.ctout.compositions
        for i in range(self.ctout.ns):
            if compositions[i][ielem] <= 1:
                smalls.append(i)
        return smalls
        
    def getConnex(self):
        name = self.fluxtypes[self.fluxtype_choice.GetCurrentSelection()]
        index = self.xvar_slider.getIndex()
        prior_included = self.prior_included_ctrl.IsChecked()
        remove_unconnected = self.remove_unconnected_ctrl.IsChecked()
        rtol = float(self.rtol_ctrl.GetValue())
        search_direction_str = self.search_direction.GetStringSelection()
        if name == "Mass": 
            str1 = "mass"
            str2 = "Mass Flux"
            unit = "g/cm3-s"
        else:
            str1 = "element"
            str2 = "Element Flux for %s" % (name,)
            unit = "atoms/cm3-s"
        
        self.message('Constructing net %s flux matrix...' % (str1,))
        include = self.include_choices.getCheckedStrings()
        exclude = self.exclude_choices.getCheckedStrings()
        included_only = False
        if search_direction_str == "Both": search_direction = 0
        elif search_direction_str == "Consumption": search_direction = 1
        elif search_direction_str == "Formation": search_direction = -1
        elif search_direction_str == "None": search_direction = 0; included_only = True
        
        if name == "Mass":
            fracs = self.ctout.getMassFractionAt(index)
            mat = self.ctout.getNetMassFluxMatrixAt(index)
        else:
            fracs = self.ctout.getElementFractionAt(index, name)
            mat = self.ctout.getNetElementFluxMatrixAt(index, name)
        
        spdict, connex = \
                self.ctout.filterFluxFromMatrix(mat, include=include,
                                               exclude=exclude,
                                               prior_included=prior_included,
                                               included_only=included_only,
                                               search_direction=search_direction,
                                               remove_unconnected=remove_unconnected,
                                               rtol=rtol)
        if name == "Mass":
            concontrib = self.ctout.getNetMassFluxContrib(index, mat, connex)
        else:
            concontrib = self.ctout.getNetElementFluxContrib(index, name, mat, connex)
        
        self.message('Constructing net %s flux matrix... done.: ' % (str1,) + 
                     '%d species, %d connections' % (len(spdict), len(connex)))
        return fracs, mat, spdict, connex, concontrib
        

        
# main controls

class FileDropTarget(wx.FileDropTarget):
    
    def __init__(self, mainframe):
        wx.FileDropTarget.__init__(self)
        self.mainframe = mainframe
        
    def OnDropFiles(self, x, y, filenames):
        if len(filenames) != 1:
            wx.MessageBox('Multiple file drop is not yet implemented', 'Error',
                          wx.ICON_ERROR)
            return False
        else:
            self.mainframe.load(filenames[0])
        return True


menuID_LOAD = wx.NewIdRef(count=1)
menuID_RELOAD = wx.NewIdRef(count=1)
menuID_INSPECTOR = wx.NewIdRef(count=1)

class MainFrame(wx.Frame):

    title = 'Cantera Post Processer'
    first_call = True
    
    def __init__(self):
        wx.Frame.__init__(self, None, title=self.title)
        self.SetMinSize((400, 300))
        self.SetSize((720, 520))

        self.ctout = None

        self.createMenuBar()
        self.statusbar = wx.StatusBar(self)
        self.SetStatusBar(self.statusbar)

        self.panel = wx.Panel(self)
        self.notebook = wx.Notebook(self.panel)
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.onPageChanged, self.notebook)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.notebook, 1, wx.EXPAND)
        self.panel.SetSizer(sizer)
        self.panel.Layout()
        
        self.SetDropTarget(FileDropTarget(self))
        self.updataUIs()
        if self.first_call and len(sys.argv) > 1:
            fn = sys.argv[1]
            self.first_call = False
            self.load(fn)

    def createMenuBar(self):
        self.menu_file = wx.Menu()
        
        self.menu_file.Append(menuID_LOAD, '&Load...\tCtrl-L', 'load')
        self.Bind(wx.EVT_MENU, self.onLoad, id = menuID_LOAD)
        self.menu_file.Append(menuID_RELOAD, '&Reload...\tCtrl-R', 'reload')
        self.Bind(wx.EVT_MENU, self.onReload, id = menuID_RELOAD)
        self.menu_file.Append(menuID_INSPECTOR, '&Widget Inspector...\tF6', 
                              'Show the wxPython Widget Inspection Tool')
        self.Bind(wx.EVT_MENU, self.onWidgetInspector, id = menuID_INSPECTOR)
        self.menu_file.AppendSeparator()
        
        self.menu_file.Append(wx.ID_EXIT, 'E&xit', 'exit')
        self.Bind(wx.EVT_MENU, self.onExit, id = wx.ID_EXIT)
        
        self.menubar = wx.MenuBar()
        self.menubar.Append(self.menu_file, '&File')
        self.SetMenuBar(self.menubar)
        
    def updateMenu(self):
        self.menu_file.Enable(menuID_RELOAD, self.ctout is not None)

    def updataUIs(self):
        self.notebook.DeleteAllPages()
        if self.ctout is None:
            pass
        else:
            self.notebook.AddPage(InfoPanel(self.notebook, self), 'Info')
            self.notebook.AddPage(ConcPanel(self.notebook, self), 'Concentration')
            self.notebook.AddPage(RatesPanel(self.notebook, self), 'Rates')
            self.notebook.AddPage(ProfilePanel(self.notebook, self), 'Profile')
            self.notebook.AddPage(ProductionPanel(self.notebook, self), 'Production')
            self.notebook.AddPage(HeatPanel(self.notebook, self), 'Heat')
            if self.ctout.sens:
                self.notebook.AddPage(SensPanel(self.notebook, self), 'Sensitivity')
            self.notebook.AddPage(PathFluxPanel(self.notebook, self), 'PathFlux')
            self.SetTitle('%s - %s' % (self.ctout.datfile, self.title))
        self.updateMenu()
        

    def message(self, text): self.statusbar.SetStatusText(text)
    def clearMessage(self): self.message('')
            
    def onPageChanged(self, event):
        self.clearMessage()

    def onLoad(self, event):
        dlg = wx.FileDialog(self,
                            wildcard='Data file (*.dat)|*.dat|All files (*.*)|*.*',
                            style=wx.FD_FILE_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK: self.load(dlg.GetPath())
        #dlg.Destroy()

    def onReload(self, event):
        self.load(self.ctout.datfile)
        
    def onWidgetInspector(self, event):
        InspectionTool().Show()
        
    def onExit(self, event):
        self.Destroy()

    def load(self, datfile):
        wx.BeginBusyCursor()
        bi = wx.BusyInfo('Loading, please wait...', self)
        try:
            self.ctout = CanteraOutput(datfile)
            self.updataUIs()
            self.message('Loaded.')
        except:
            wx.MessageBox(traceback.format_exc(), 'Traceback', wx.ICON_ERROR)
        #bi.Destroy()
        wx.EndBusyCursor()


class App(wx.App):
    
    def __init__(self):
        wx.App.__init__(self)
        
    def OnInit(self):
        frame = MainFrame()
        frame.Show(True)
        return True

if __name__ == '__main__':
    app = App()
    app.MainLoop()
