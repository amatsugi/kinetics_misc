#! /usr/bin/env python3

import wx

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties 
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2WxAgg


def genColorAndSymbols():
    count = 0
    symbols = "os^Dhpv<>Hd"
    numsymbols = len(symbols)
    colors = "bgrcmy"
    numcolors = len(colors)
    while True:
        yield symbols[count % numsymbols] + colors[count % numcolors]
        count += 1

def genColorAndStyles():
    count = 0
    colors = ['black', 'blue', 'green', 'red', 'cyan', 'magenta', 'yellow']
    numcolors = len(colors)
    styles = ['-', '--', '-.', ':']
    numstyles = len(styles)
    while True:
        div, mod = divmod(count, numcolors)
        yield colors[mod], styles[div % numstyles]
        count += 1


class PlotFrame(wx.Frame):
    
    legend_locs = ['No legend', 
                   'upper right', 'upper left', 'lower left', 'lower right',
                   'right', 'center left', 'center right', 
                   'lower center', 'upper center', 'center']
    scales = ['lin-lin', 'xlog', 'ylog', 'log-log']
    
    def __init__(self, parent, id=wx.ID_ANY, title='Plot', figsize=(6,4)):
        wx.Frame.__init__(self, parent=parent, id=-id, title=title)
        self.panel = wx.Panel(self)
        
        self.figure = Figure(figsize=figsize, tight_layout=True)
        self.canvas = FigureCanvasWxAgg(self.panel, wx.ID_ANY, self.figure)
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)

        # additional tools
        tsize = (16,16)
        _NTB_copy = wx.NewIdRef(count=1)
        _NTB_grid = wx.NewIdRef(count=1)
        _NTB_legend = wx.NewIdRef(count=1)
        _NTB_scale = wx.NewIdRef(count=1)
        _NTB_msg = wx.NewIdRef(count=1)
        copy_bmp = wx.ArtProvider.GetBitmap(wx.ART_COPY, wx.ART_TOOLBAR, tsize)
        grid_bmp = wx.ArtProvider.GetBitmap(wx.ART_CROSS_MARK, wx.ART_TOOLBAR, tsize)
        
        self.toolbar.AddSeparator()
        self.toolbar.AddTool(_NTB_copy, "copy", copy_bmp, 'Copy to clipboard')
        self.toolbar.Bind(wx.EVT_TOOL, self.onCopyToClipboard, id=_NTB_copy)
        self.toolbar.AddCheckTool(_NTB_grid, "grid", grid_bmp, shortHelp='Grid', longHelp='Toggle grid')
        self.toolbar.Bind(wx.EVT_TOOL, self.onToggleGrid, id=_NTB_grid)

        self.toolbar.AddSeparator()
        self.legend_choice = wx.Choice(self.toolbar, _NTB_legend,
                                       choices = self.legend_locs)
        self.legend_choice.SetSelection(0)
        self.loc = self.legend_locs[0]
        self.toolbar.AddControl(self.legend_choice)
        self.toolbar.Bind(wx.EVT_CHOICE, self.onChoiceLegend, id=_NTB_legend)

        self.toolbar.AddSeparator()
        self.scale_choice = wx.Choice(self.toolbar, _NTB_scale,
                                      choices = self.scales)
        self.scale_choice.SetSelection(0)
        self.scale = self.scales[0]
        self.toolbar.AddControl(self.scale_choice)
        self.toolbar.Bind(wx.EVT_CHOICE, self.onChoiceScale, id=_NTB_scale)

        self.statusBar = wx.StatusBar(self, -1)
        self.SetStatusBar(self.statusBar)
        self.canvas.mpl_connect('motion_notify_event', self.UpdateStatusBarXY)
        
        self.toolbar.Realize()
        self.toolbar.Show()
        
        box = wx.BoxSizer(wx.VERTICAL)
        box.Add(self.canvas, 1, wx.EXPAND)
        box.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.panel.SetSizer(box)
        self.panel.Layout()
        self.panel.Fit()
        w, h = self.panel.GetEffectiveMinSize()
        self.SetSize((w+12,h+60))

    def UpdateStatusBarXY(self, event):
        if event.inaxes:
            self.statusBar.SetStatusText("x=%g  y=%g" %(event.xdata, event.ydata))
            
    def onCopyToClipboard(self, event):
        self.canvas.Copy_to_Clipboard(event)
        self.statusBar.SetStatusText('Copied.')

    def onToggleGrid(self, event):
        on = event.IsChecked()
        for a in self.figure.get_axes():
            a.grid(on)
        self.statusBar.SetStatusText('Grid %s.' % ('on' if on else 'off',))
        self.draw()
    
    def onChoiceLegend(self, event):
        self.loc = event.GetString()
        self.draw()
    
    def onChoiceScale(self, event):
        self.scale = event.GetString()
        self.draw()
        
    def getFigure(self):
        return self.figure

    def addSubplot(self, *args, **kwargs):
        return self.figure.add_subplot(*args, **kwargs)

    def setLegend(self, num):
        self.legend_choice.SetSelection(num)
        self.loc = self.legend_locs[num]

    def setScale(self, num):
        self.scale_choice.SetSelection(num)
        self.scale = self.scales[num]
    
    def draw(self):
        for a in self.figure.get_axes():
            if self.loc == 'No legend':
                a.legend_ = None
            else:
                a.legend(loc=self.loc, prop=FontProperties(size='smaller'))
            if self.scale == 'lin-lin':
                a.set_xscale('linear')
                a.set_yscale('linear')
            elif self.scale == 'xlog':
                a.set_xscale('log', nonposx='mask')
                a.set_yscale('linear')
            elif self.scale == 'ylog':
                a.set_xscale('linear')
                a.set_yscale('log', nonposy='mask')
            elif self.scale == 'log-log':
                a.set_xscale('log', nonposx='mask')
                a.set_yscale('log', nonposy='mask')
        self.canvas.draw()


class TestFrame(wx.Frame):
    def __init__(self, parent=None, id=wx.ID_ANY, title='Plot'):
        wx.Frame.__init__(self, parent=parent, id=id, title=title)
        
        self.panel = wx.Panel(self)
        
        btn = wx.Button(self.panel, label='Plot')
        self.Bind(wx.EVT_BUTTON, self.onPlot, btn)
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(btn, 0, wx.ALL, 4)

        self.panel.SetSizer(vbox)
        vbox.Fit(self)

    def onPlot(self, event):
        import numpy as np
        x = np.arange(0.0, 5.0, 0.10)
        
        frame = PlotFrame(self, title='sin cos')
        ax = frame.addSubplot(111, xlabel='x', ylabel='y')
        
        clrsty = genColorAndStyles()
        clr, sty = next(clrsty)
        ax.plot(x, np.sin(x), color=clr, ls=sty, lw=1.5, label='sin(x)')
        clr, sty = next(clrsty)
        ax.plot(x, np.cos(x), color=clr, ls=sty, lw=1.5, label='cos(x)')

        frame.draw()
        frame.Show()


if __name__ == "__main__":
    app = wx.App()
    frame = TestFrame(title='test plot flame')
    frame.Show(True)
    app.MainLoop()


