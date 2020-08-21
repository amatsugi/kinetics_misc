#! /usr/bin/env python3

import random
import subprocess
try: from cStringIO import StringIO
except: from io import StringIO
import wx
import wx.lib.ogl as ogl

class SpeciesEvtHandler(ogl.ShapeEvtHandler):
    def __init__(self, fluxflame, msg_st="", msg_pop=""):
        ogl.ShapeEvtHandler.__init__(self)
        self.fluxflame = fluxflame
        self.msg_st = msg_st
        self.msg_pop = msg_pop

    def update(self, clearmsg=False, adjustpos=False):
        if clearmsg: self.fluxflame.clearMessage()
        else: self.fluxflame.message("%s" % self.msg_st)
        if adjustpos: self.fluxflame.adjustCurrentPos()
        self.fluxflame.updateScrollbars()

    def select(self):
        shape = self.GetShape()
        canvas = shape.GetCanvas()
        dc = wx.ClientDC(canvas)
        canvas.PrepareDC(dc)
        if shape.Selected():
            shape.Select(False, dc)
            self.update(True)
        else:
            shapeList = canvas.GetDiagram().GetShapeList()
            toUnselect = []
            for s in shapeList:
                if s.Selected(): toUnselect.append(s)
            shape.Select(True, dc)
            if toUnselect:
                for s in toUnselect: s.Select(False, dc)
            self.update()

    def OnRightClick(self, x, y, keys=0, attachment=0):
        shape = self.GetShape()
        win = wx.PopupTransientWindow(self.fluxflame, wx.SIMPLE_BORDER)
        win.SetBackgroundColour("#cceeee")
        st = wx.StaticText(win, label=self.msg_pop, pos=(5,5))
        sz = st.GetBestSize()
        win.SetSize((sz.width+10, sz.height+10))
        win.Position(wx.GetMousePosition(), (0, 0))
        win.Popup()
        if shape.Selected(): self.update()
        else: self.select()

    def OnLeftClick(self, x, y, keys=0, attachment=0):
        self.select()

    def OnEndDragLeft(self, x, y, keys=0, attachment=0):
        shape = self.GetShape()
        ogl.ShapeEvtHandler.OnEndDragLeft(self, x, y, keys, attachment)
        if not shape.Selected():
            self.OnLeftClick(x, y, keys, attachment)
        if ((shape.GetX() - shape.GetWidth()/2. <= 0) or 
            (shape.GetY() - shape.GetHeight()/2. <= 0)):
            self.update(adjustpos=True)
        else:
            self.update()

    def OnSizingEndDragLeft(self, pt, x, y, keys, attch):
        shape = self.GetShape()
        ogl.ShapeEvtHandler.OnSizingEndDragLeft(self, pt, x, y, keys, attch)
        canvas = shape.GetCanvas()
        dc = wx.ClientDC(canvas)
        canvas.PrepareDC(dc)
        shape.Recentre(dc)
        self.update()


class ConnexEvtHandler(ogl.ShapeEvtHandler):
    
    def __init__(self, fluxflame, msg_st="", msg_pop=""):
        ogl.ShapeEvtHandler.__init__(self)
        self.fluxflame = fluxflame
        self.msg_st = msg_st
        self.msg_pop = msg_pop

    def update(self, clearmsg=False):
        if clearmsg: self.fluxflame.clearMessage()
        else: self.fluxflame.message("%s" % self.msg_st)
        self.fluxflame.updateScrollbars()

    def select(self):
        shape = self.GetShape()
        canvas = shape.GetCanvas()
        dc = wx.ClientDC(canvas)
        canvas.PrepareDC(dc)
        if shape.Selected():
            shape.Select(False, dc)
            self.update(True)
        else:
            shapeList = canvas.GetDiagram().GetShapeList()
            toUnselect = []
            for s in shapeList:
                if s.Selected(): toUnselect.append(s)
            shape.Select(True, dc)
            if toUnselect:
                for s in toUnselect: s.Select(False, dc)
            self.update()

    def OnRightClick(self, x, y, keys=0, attachment=0):
        shape = self.GetShape()
        win = wx.PopupTransientWindow(self.fluxflame, wx.SIMPLE_BORDER)
        win.SetBackgroundColour("#eeeecc")
        st = wx.StaticText(win, label=self.msg_pop, pos=(5,5))
        sz = st.GetBestSize()
        win.SetSize((sz.width+10, sz.height+10))
        win.Position(wx.GetMousePosition(), (0, 0))
        win.Popup()
        if shape.Selected(): self.update()
        else: self.select()

    def OnLeftClick(self, x, y, keys=0, attachment=0):
        self.select()

    def OnDraw(self, dc):
        line = self.GetShape()
        size = line.GetPen().GetWidth()
        # adjust arrow head
        x1, y1, x2, y2 = line.FindLineEndPoints()
        # shift end point by size * 2/3 to eliminate arrow overlap
        shift = size * 1.5
        # also, offset both ends by 2 point
        offset = 2
        lx, ly = x2-x1, y2-y1
        ld = (lx*lx + ly*ly)**0.5
        if ld < 0.1:
            pass
        else:
            x1 = x1 + offset * (lx / ld)
            y1 = y1 + offset * (ly / ld)
            x2 = x2 - (shift + offset) * (lx / ld)
            y2 = y2 - (shift + offset) * (ly / ld)
        line.SetEnds(x1, y1, x2, y2)
        ogl.ShapeEvtHandler.OnDraw(self, dc)
        

class FluxFrame(wx.Frame):

    colors = ["Color1", "Color2", "Grayscale", "Black"]
    colors_tooltips = ["Color species by fractions",
                       "Color species by net fluxes",
                       "Grayscale species by fractions",
                       "Black"]
    layouts = ["Grid", "Random", "dot"]
    layouts_tooltips = ["Grid layout",
                        "Random layout",
                        "Hierarchical layout (Graphviz-dot)"]
    vmargin = 10
    hmargin = 10
    default_height = 24
    default_width = 80
    
    def __init__(self, parent, id=wx.ID_ANY, title='Flux'):
        ogl.OGLInitialize()
        wx.Frame.__init__(self, parent=parent, id=id, title=title)
        self.SetMinSize((400, 300))
        self.SetSize((640, 520))

        self.statusbar = wx.StatusBar(self)
        self.SetStatusBar(self.statusbar)
        
        self.panel = wx.Panel(self)
        self.canvas = ogl.ShapeCanvas(self.panel)
        self.canvas.SetBackgroundColour(wx.WHITE)
        self.diagram = ogl.Diagram()
        self.canvas.SetDiagram(self.diagram)
        self.diagram.SetCanvas(self.canvas)

        self.directory = "."
        self.imgfilename = "flux_image"
        self.spdict = {}
        self.species_names = []
        self.fracs = []
        self.connex = {}
        self.concontrib = {}
        self.reactions_names = []
        self.unit = ""
        self.maxfrac = 0.0
        self.maxnetflux = 0.0
        self.maxflux = 0.0

        self.font = wx.Font(10, wx.SWISS, wx.NORMAL, wx.FONTWEIGHT_BOLD)
        self.species_shapes = {}
        self.connex_lines = {}
        self.connex_arrows = {}
        
        # control
        
        self.color_ctrl = wx.RadioBox(self.panel, label="color", choices=self.colors,
                                      style=wx.RA_VERTICAL)
        for i, s in enumerate(self.colors_tooltips):
            self.color_ctrl.SetItemToolTip(i, s)
        self.Bind(wx.EVT_RADIOBOX, self.onChangeColorSize, self.color_ctrl)
        self.penwidth_text = wx.StaticText(self.panel, label='Penwidth: %d' % (10,))
        self.penwidth_slider = wx.Slider(self.panel, minValue=1, maxValue=30, value=10,
                                         size=(80, -1), style = wx.SL_HORIZONTAL)
        self.Bind(wx.EVT_COMMAND_SCROLL, self.onChangeColorSize, self.penwidth_slider)
        
        self.draw_layout = wx.RadioBox(self.panel, label="layout", choices=self.layouts,
                                       style=wx.RA_VERTICAL)
        for i, s in enumerate(self.layouts_tooltips):
            self.draw_layout.SetItemToolTip(i, s)
        self.draw_button = wx.Button(self.panel, label='Redraw')
        self.Bind(wx.EVT_BUTTON, self.onReDraw, self.draw_button)
        
        self.scale_xm10_button = wx.Button(self.panel, label='-x', size=(36, -1))
        self.Bind(wx.EVT_BUTTON, self.onScaleXM10, self.scale_xm10_button)
        self.scale_xp10_button = wx.Button(self.panel, label='+x', size=(36, -1))
        self.Bind(wx.EVT_BUTTON, self.onScaleXP10, self.scale_xp10_button)
        self.scale_ym10_button = wx.Button(self.panel, label='-y', size=(36, -1))
        self.Bind(wx.EVT_BUTTON, self.onScaleYM10, self.scale_ym10_button)
        self.scale_yp10_button = wx.Button(self.panel, label='+y', size=(36, -1))
        self.Bind(wx.EVT_BUTTON, self.onScaleYP10, self.scale_yp10_button)
        self.sshot_button = wx.Button(self.panel, label="SShot")
        self.Bind(wx.EVT_BUTTON, self.onTakeScreenShot, self.sshot_button)

        # layout
        
        penwidth_box = wx.BoxSizer(wx.VERTICAL)
        penwidth_box.Add(self.penwidth_text, 0, border=0)
        penwidth_box.Add(self.penwidth_slider, 0, border=0)
        xscale_box = wx.BoxSizer(wx.HORIZONTAL)
        xscale_box.Add(self.scale_xm10_button, 0, wx.RIGHT, 2)
        xscale_box.Add(self.scale_xp10_button, 0, wx.LEFT, 2)
        yscale_box = wx.BoxSizer(wx.HORIZONTAL)
        yscale_box.Add(self.scale_ym10_button, 0, wx.RIGHT, 2)
        yscale_box.Add(self.scale_yp10_button, 0, wx.LEFT, 2)
        
        config_box = wx.BoxSizer(wx.VERTICAL)
        config_box.Add(self.color_ctrl, 0, wx.ALL, 2)
        config_box.Add(penwidth_box, 0, wx.ALL, 2)
        config_box.Add(self.draw_layout, 0, wx.ALL, 2)
        config_box.Add(self.draw_button, 0, wx.ALL, 2)
        config_box.AddSpacer(4)
        config_box.Add(xscale_box, 0, wx.ALL, 2)
        config_box.Add(yscale_box, 0, wx.ALL, 2)
        config_box.AddSpacer(4)
        config_box.Add(self.sshot_button, 0, wx.ALL, 2)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(config_box, 0, wx.ALL, 4)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        
        self.panel.SetSizer(sizer)
        self.panel.Layout()

    def onChangeColorSize(self, event):
        self.setColorSize()
        self.panel.Layout()

    def onReDraw(self, event): 
        layout = self.draw_layout.GetStringSelection()
        self.redraw(layout)
        self.panel.Layout()

    def onTakeScreenShot(self, event):
        self.takeScreenShot()
        self.panel.Layout()

    def onScaleXP10(self, event): self.scaleCurrentPos(1.1, 1.); self.panel.Layout()
    def onScaleXM10(self, event): self.scaleCurrentPos(0.9, 1.); self.panel.Layout()
    def onScaleYP10(self, event): self.scaleCurrentPos(1., 1.1); self.panel.Layout()
    def onScaleYM10(self, event): self.scaleCurrentPos(1., 0.9); self.panel.Layout()
    
    def setMisc(self, directory, imgfilename):
        self.directory = directory
        self.imgfilename = imgfilename
        
    def setFlux(self, spdict, species_names, fracs, connex, concontrib, 
                reactions_names, unit="atoms/cm3-s"):
        #  spdict[id]: (sumfor, sumback, [(forind1, frac), ...], [(backind1, frac, ...)])
        self.spdict = spdict
        self.species_names = species_names
        self.fracs = fracs
        self.connex = connex   # connex[(from,to)] = (flux, consump frac, form frac)
        self.concontrib = concontrib   # concontrib[(from,to)]: [(rind, frac)]
        self.reactions_names = reactions_names
        self.unit = unit
        self.maxfrac = 0.0
        self.maxnetflux = 0.0
        self.maxflux = 0.0
        for i in self.spdict:
            if self.fracs[i] > self.maxfrac: self.maxfrac = self.fracs[i]
            netflux = self.spdict[i][1] - self.spdict[i][0]
            if abs(netflux) > self.maxnetflux: self.maxnetflux = abs(netflux)
        for x in self.connex:
            if self.connex[x][0] > self.maxflux: self.maxflux = self.connex[x][0]
    
    def message(self, text): self.statusbar.SetStatusText(text)
    def clearMessage(self): self.message('')

    def updateScrollbars(self):
        posx, posy = self.canvas.GetViewStart()
        maxx = 0
        maxy = 0
        for i in self.species_shapes:
            shape = self.species_shapes[i]
            x, y = shape.GetX(), shape.GetY()
            w, h = shape.GetWidth(), shape.GetHeight()
            xx = x + w/2. + self.hmargin
            yy = y + h/2. + self.vmargin
            if maxx < xx: maxx = xx
            if maxy < yy: maxy = yy
        self.canvas.SetScrollbars(20, 20, maxx/20, maxy/20, posx, posy)
        self.canvas.Refresh(True)
        self.panel.Layout()
        
    def draw(self, layout):
        self.diagram.DeleteAllShapes()
        self.message("Drawing flux with %s layout method..." % (layout,))
        if layout == "Grid":
            pos = self.genPosGrid()
        elif layout == "Random":
            pos = self.genPosRandom()
        elif layout == "dot":
            pos = self.genPosGraphvizDot()
        else:
            raise ValueError("Unknown layout layout: %s" % layout)
        self.draw_layout.SetStringSelection(layout)
        pos = self.adjustPos(pos)
        
        for i in self.spdict:
            x, y = pos[i]
            name = self.species_names[i]
            frac = self.fracs[i]
            
            msg_st = "%d: %s  (%.3g %%)" % (i+1, name, frac*100)
            msg_pop = "%s [id=%d]\n" % (name, i+1)
            msg_pop += "Fraction = %.3g %%\n" % (frac*100,)
            msg_pop += "Consumption flux = %.4e %s\n" % (self.spdict[i][0], self.unit)
            for l in self.spdict[i][2]:
                msg_pop += "    -> %s [id=%d]  (%.3g %%)\n" % \
                           (self.species_names[l[0]], l[0]+1, l[1]*100)
            msg_pop += "Formation flux = %.4e %s\n" % (self.spdict[i][1], self.unit)
            for l in self.spdict[i][3]:
                msg_pop += "    <- %s [id=%d]  (%.3g %%)\n" % \
                           (self.species_names[l[0]], l[0]+1, l[1]*100)
            msg_pop += "Total flux = %.4e %s" % \
                       (self.spdict[i][1] - self.spdict[i][0], self.unit)
            
            shape = self.addSpeciesShape(x, y, text=name, msg_st=msg_st, msg_pop=msg_pop)
            self.species_shapes[i] = shape
            
        for x in self.connex:
            k, l = x
            flux, cfr, ffr = self.connex[x]
            contrib = self.concontrib[x]
            froms, tos = self.species_names[k], self.species_names[l]
            msg_st = "%s (%.3g %%) -> %s (%.3g %%)  Flux = %.2e %s" % \
                     (froms, cfr*100, tos, ffr*100, flux, self.unit)
            msg_pop = "%s [id=%d] -> %s [id=%d]\n" % (froms, k+1, tos, l+1)
            msg_pop += "Flux = %.4e %s\n" % (flux, self.unit)
            msg_pop += "Account for:\n    %.3g %% of %s consumption\n" % (cfr*100, froms)
            msg_pop += "    %.3g %% of %s formation\n" % (ffr*100, tos)
            if len(contrib) > 0:
                msg_pop += "Contributions:"
                for rind, rfrac in contrib:
                    # ignore < 0.1% contribution
                    if abs(rfrac) < 0.001: continue
                    msg_pop += "\n    %d: %s  (%.3g %%)" % \
                               (rind+1, self.reactions_names[rind], rfrac*100)
            line, arrow = self.addConnexLine(self.species_shapes[k], self.species_shapes[l],
                                             msg_st=msg_st, msg_pop=msg_pop)
            self.connex_lines[x] = line
            self.connex_arrows[x] = arrow
        self.setColorSize()
        self.diagram.ShowAll(True)
        self.message("Drawing flux with %s layout method... done" % (layout,))
        self.updateScrollbars()

    def addSpeciesShape(self, x, y, text="", msg_st="", msg_pop=""):
        shape = ogl.EllipseShape(self.default_width, self.default_height)
        shape.SetX(x)
        shape.SetY(y)
        shape.AddText(text)
        shape.SetFont(self.font)
        self.diagram.AddShape(shape)
        
        evthandler = SpeciesEvtHandler(self, msg_st, msg_pop)
        evthandler.SetShape(shape)
        evthandler.SetPreviousHandler(shape.GetEventHandler())
        shape.SetEventHandler(evthandler)

        return shape

    def addConnexLine(self, start, end, msg_st="", msg_pop=""):
        line = ogl.LineShape()
        line.MakeLineControlPoints(2)
        arrow = line.AddArrow(ogl.ARROW_ARROW)
        start.AddLine(line, end)
        self.diagram.AddShape(line)

        evthandler = ConnexEvtHandler(self, msg_st, msg_pop)
        evthandler.SetShape(line)
        evthandler.SetPreviousHandler(line.GetEventHandler())
        line.SetEventHandler(evthandler)
        
        return line, arrow
        
    def setColorSize(self):
        color = self.color_ctrl.GetStringSelection()
        penwidth_width = self.penwidth_slider.GetValue()
        self.penwidth_text.SetLabel('Penwidth: %d' % (penwidth_width,))

        def rgbinterp(rgbmin, rgbmax, frac):
            rgb = []
            for i in range(3):
                tmp = int(frac * (rgbmax[i] - rgbmin[i]) + rgbmin[i])
                if tmp < 0: tmp = 0
                elif tmp > 255: tmp = 255
                rgb.append(tmp)
            return tuple(rgb)
        
        for i in self.species_shapes:
            shape = self.species_shapes[i]
            if color == "Color1":
                frac = self.fracs[i]
                rgb = rgbinterp((255, 255, 255), (255, 0, 0), (frac/self.maxfrac)**0.5)
            elif color == "Color2":
                netflux = self.spdict[i][1] - self.spdict[i][0]
                if netflux > 0.: maxrgb = (255, 50, 50)
                else: maxrgb = (100, 100, 255)
                rgb = rgbinterp((255, 255, 255), maxrgb, abs(netflux/self.maxnetflux)**0.5)
            elif color == "Grayscale":
                frac = self.fracs[i]
                rgb = rgbinterp((255, 255, 255), (180, 180, 180), (frac/self.maxfrac)**0.5)
            else:
                rgb = (255, 255, 255)
            pen = wx.Pen(wx.BLACK, 1)
            pen.SetCap(wx.CAP_BUTT)
            pen.SetJoin(wx.JOIN_MITER)
            shape.SetPen(pen)
            shape.SetBrush(wx.Brush(rgb, wx.SOLID))
                
        for x in self.connex_lines:
            flux, cfr, ffr = self.connex[x]
            size = max(penwidth_width/20., penwidth_width*abs(flux/self.maxflux)**0.5)
            arrsize = 10. + size/5.
            line = self.connex_lines[x]
            arrow = self.connex_arrows[x]
            arrow.SetSize(arrsize)
            if color == "Color1" or color == "Color2":
                rgb = rgbinterp((120, 120, 180), (240, 20, 20), (flux/self.maxflux)**0.5)
            elif color == "Grayscale":
                rgb = rgbinterp((150, 150, 150), (0, 0, 0), (flux/self.maxflux)**0.5)
            else:
                rgb = (0, 0, 0)
            pen = wx.Pen(rgb, size)
            pen.SetCap(wx.CAP_BUTT)
            pen.SetJoin(wx.JOIN_MITER)
            line.SetPen(pen)
            line.SetBrush(wx.Brush(rgb, wx.SOLID))
        self.updateScrollbars()
        return

    def redraw(self, layout):
        self.message("Re-Drawing flux with %s layout method..." % (layout,))
        if layout == "Grid":
            pos = self.genPosGrid()
        elif layout == "Random":
            pos = self.genPosRandom()
        elif layout == "dot":
            pos = self.genPosGraphvizDot()
        else:
            raise ValueError("Unknown layout layout: %s" % layout)
        self.draw_layout.SetStringSelection(layout)
        pos = self.adjustPos(pos)
        self.changePos(pos)
        self.message("Re-Drawing flux with %s layout method... done" % (layout,))
        
    def genPosGrid(self):
        x, y = 0, 0
        pos = {}
        for i in self.spdict:
            pos[i] = [x, y]
            if x > 400: x = 0; y += 50
            else: x += 120
        return pos

    def genPosRandom(self):
        pos = {}
        if len(self.spdict) < 15: maxx, maxy = 500, 400
        else: maxx, maxy = 800, 600
        for i in self.spdict:
            x = random.random() * maxx
            y = random.random() * maxy
            pos[i] = [x, y]
        return pos
    
    def genPosGraphvizDot(self):
        si = StringIO()
        
        si.write("digraph g {\n")
        for i in self.spdict: si.write('  "%d";\n' % (i,))
        for x in self.connex:
            flux, cfr, ffr = self.connex[x]
            # weight is integer
            # otherwise graphviz-dot sometimes ends with seg. fault
            weight = max(1, int(100.*abs(flux/self.maxflux)))
            si.write('  "%d" -> "%d" [weight=%d];\n' % (x[0], x[1], weight))
        si.write("}\n")
        si.write("\n")
        instr = si.getvalue()
        si.close()
        
        p = subprocess.Popen(["dot", "-Tplain"], shell=True, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res = p.communicate(instr.encode())[0]
        so = StringIO()
        so.write(res.decode())
        so.seek(0)
        pos = {}
        try:
            ll = next(so).split()
            maxx, maxy = float(ll[2]), float(ll[3])
            if ll[2] == "0": maxx = 1
            if ll[3] == "0": maxy = 1
            scale = 96. # inch -> point
            scalex, scaley = scale, scale*0.8
            for l in so:
                ll = l.split()
                if ll[0] == "node":
                    ind = int(ll[1])
                    x, y = float(ll[2]), float(ll[3])
                    pos[ind] = [x*scalex, (maxy-y)*scaley]
        except Exception as err:
            raise ValueError("Error while parsing output: %s\n[IN]\n%s\n[OUT]\n%s" % 
                             (err, instr, res))
        so.close()
        return pos

    def getCurrentPos(self):
        pos = {}
        for i in self.species_shapes:
            shape = self.species_shapes[i]
            pos[i] = [shape.GetX(), shape.GetY()]
        return pos

    def changePos(self, pos):
        dc = wx.ClientDC(self.canvas)
        self.canvas.PrepareDC(dc)
        for i in pos:
            x, y = pos[i]
            self.species_shapes[i].Move(dc, x, y)
        self.updateScrollbars()
        
    def adjustPos(self, pos):
        minx, miny = None, None
        for i in pos:
            x, y = pos[i]
            if minx is None: minx = x
            if miny is None: miny = y
            xx = x - self.default_width/2. - self.hmargin
            yy = y - self.default_height/2. - self.vmargin
            if minx > xx: minx = xx
            if miny > yy: miny = yy
        for i in pos:
            x, y = pos[i]
            pos[i] = [x-minx, y-miny]
        return pos

    def adjustCurrentPos(self):
        pos = self.getCurrentPos()
        self.changePos(self.adjustPos(pos))

    def scaleCurrentPos(self, scalex, scaley):
        pos = self.getCurrentPos()
        for i in pos:
            x, y = pos[i]
            pos[i] = x*scalex, y*scaley
        self.changePos(self.adjustPos(pos))

    def takeScreenShot(self):
        bmp_types = {".bmp": wx.BITMAP_TYPE_BMP,
                     ".gif": wx.BITMAP_TYPE_GIF,
                     ".jpg": wx.BITMAP_TYPE_JPEG,
                     ".pcx": wx.BITMAP_TYPE_PCX,
                     ".png": wx.BITMAP_TYPE_PNG,
                     ".pnm": wx.BITMAP_TYPE_PNM,
                     ".tif": wx.BITMAP_TYPE_TIF,
                     ".xbm": wx.BITMAP_TYPE_XBM,
                     ".xpm": wx.BITMAP_TYPE_XPM}
        file_types = sorted(bmp_types.keys())
        wildcard = "|".join(["%s files (*%s)|*%s" % 
                             (x[1:].upper(),x,x) for x in file_types])
        dlg = wx.FileDialog(self, message="Save to file",
                            defaultDir=self.directory, 
                            defaultFile=self.imgfilename,
                            wildcard=wildcard,
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        dlg.SetFilterIndex(file_types.index(".png"))
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            format = bmp_types[file_types[dlg.GetFilterIndex()]]
            dlg.Destroy()
        else:
            dlg.Destroy()
            return

        #w, h = self.canvas.GetVirtualSizeTuple()
        w, h = self.canvas.GetVirtualSize()
        dc = wx.MemoryDC()
        #bitmap = wx.EmptyBitmap(w+2*self.hmargin, h+2*self.vmargin)
        bitmap = wx.Bitmap(w+2*self.hmargin, h+2*self.vmargin)
        dc.SelectObject(bitmap)
        dc.SetBackground(wx.WHITE_BRUSH)
        dc.Clear()
        self.canvas.Redraw(dc)
        bitmap.SaveFile(filename, format)

def dummydata(N=10):
    spinds = list(range(N))
    names = ["s%d" % i for i in spinds]
    reactions = ["r%d" % i for i in range(10)]
    fracs = [random.random()**5. for i in spinds]
    
    connex = {}
    concontrib = {}
    for i in range(N):
        for j in range(i+1, N):
            if random.random() > 0.7:
                if random.random() > 0.2:
                    connex[(i, j)] = (random.random()**4, random.random(), random.random())
                    concontrib[(i, j)] = [(1, 0.2), (3,0.8)]
                else:
                    connex[(j, i)] = (random.random()**4, random.random(), random.random())
                    concontrib[(j, i)] = [(2, 0.2), (3,0.8)]
    spdicts = {}
    for k in spinds:
        spdicts[k] = (random.random(), random.random(), 
                      [(8, 0.2), (4, 0.8)], [(6, 0.1), (1, 0.2)])
    return spdicts, names, fracs, connex, concontrib, reactions


if __name__ == "__main__":
    app = wx.App()
    frame = FluxFrame(None, title='test flame')

    frame.setFlux(*dummydata(12))
    frame.draw("dot")
    
    frame.Show(True)
    
    app.MainLoop()



