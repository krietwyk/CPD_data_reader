"""
CPD Visualisation and Interpolation Module - 210214

A simple script for visualising and interpolating CPD measurements
performed on a KP Technology scanning Kelvin probe. The script interpolates
2D grids of uniform spacing using an inverse distance approach.

Data can be opened directly from a file or as a tab-delimited array using
the input array text box and button.
"""

from tkinter import (Frame, N, E, W, NW, StringVar, DISABLED, Text,
                     IntVar, DoubleVar, END, Canvas, PhotoImage, ttk,
                     filedialog, Tk, messagebox, Spinbox, Checkbutton)

from os import path

import tempfile

import csv

from numpy import (transpose, asarray, sqrt, append, savetxt, ones,
                   amin, amax)

import matplotlib.pyplot as plt

def power_(list_, exponent_):
    """Performs power operation of each element in a list."""
    return [i ** exponent_ for i in list_]

# Prevent ipython from creating a figure window
plt.ioff()

# Use temporary files to create images in the GUI canvas
tmpdir = tempfile.mkdtemp()
im_name1 = 'example01.png'
im_name2 = 'example02.png'

class Application(Frame):
    def __init__(self, master=None):
        """ Main GUI."""
        Frame.__init__(self, master)
        Frame.grid(self, padx = (40, 40), pady = (10, 10))

        self.path1 = StringVar()
        path1_lbl = ttk.Label(self, text='Path:')
        path1_lbl.grid(column = 1, columnspan = 4, row = 1, sticky=W)
        path1_entry = ttk.Entry(self, width = 30, textvariable = self.path1)
        path1_entry.grid(column = 1, columnspan = 4, row = 2, sticky = (W, E))
        path1_btn = ttk.Button(self, text = "Browse", command = self.open1)
        path1_btn.grid(column = 1, row = 3, sticky = W, padx = (0, 0))

        self.graph1_btn = ttk.Button(self, text = "Graph", state = DISABLED,
                                     command = lambda: 
                                         self.graph1(self.C, 1))
        self.graph1_btn.grid(column = 2, row = 3, sticky=W)

        self.info1_lbl = ttk.Label(self, text = 'Scan Info/Comments:')
        self.info1_lbl.grid(column = 1, columnspan = 4, row = 5,
                            sticky = W, pady = (10, 0))
        self.info1_Entry = Text(self, width = 30, height = 12)
        self.info1_Entry.grid(column = 1, columnspan = 4, row = 6,
                              sticky = (N, W, E))

        in_C = StringVar()
        self.input_entry = ttk.Entry(self, width = 10, textvariable=in_C)
        self.input_entry.grid(column = 4, row = 3, sticky = (W, E))
        self.int_btn = ttk.Button(self, text = "Input Array",
                     command=lambda: self.Input_C())
        self.int_btn.grid(column = 4, row = 5, sticky = W)     

        self.int_lbl = ttk.Label(self, text = "Inverse Distance "
                                 "Weighting Interpolation", font = 'bold')
        self.int_lbl.grid(column = 1, columnspan = 4, row = 7, sticky = W,
                          pady = (10, 2))        
        
        self.int_btn = ttk.Button(self, text = "Interpolate", state = DISABLED, 
                     command=lambda: self.interpolation(self.C))
        self.int_btn.grid(column = 1, row = 8, sticky = W)
        
        self.NP1 = IntVar()    
        self.NP1.set(5)
        n_lbl = ttk.Label(self, text = 'No. pts') 
        n_lbl.grid(column = 1, row = 9, sticky = (N, W))
        self.n_box = Spinbox(self, format='%.0f', from_ = 0, 
                             to = self.NP1.get(),
                             width = 5, increment = 1) #1E10
        self.n_box.grid(column = 1, row = 10, sticky = (N, W))
        self.n_box.insert(1,5)
        
        p_lbl = ttk.Label(self, text = 'Power')
        p_lbl.grid(column = 2, row = 9, sticky = (N, W))
        self.p_box = Spinbox(self, format = '%.0f', from_ = 0, to = 1E10,
                             width = 5, increment = 1)
        self.p_box.grid(column = 2, row = 10, sticky = (N, W))
        self.p_box.insert(1,2)

        self.transpose_chk = IntVar()
        self.trans_btn = Checkbutton(self,  text="Tranpose", 
                                     variable = self.transpose_chk, 
                                     state = DISABLED, 
                                     command = lambda: self.check_tr(self.C)) 
        self.trans_btn.grid(column = 3, row = 10, sticky = (N, W))

        self.Refl_x = IntVar()
        self.Refl_x.set(1)
        self.Refl_x_btn = Checkbutton(self,  text = "Reflect_x", 
                                      variable = self.Refl_x,
                                      state = DISABLED, 
                                      command = lambda: self.check_tr(self.C))
        self.Refl_x_btn.grid(column = 3, row = 9, sticky = (N, W))    
        self.info1_Entry.insert(END,"Reflect x = %s" %(self.Refl_x.get()))   
        self.N_out_btn = ttk.Button(self, text = "Save-Grid", state = DISABLED, 
                                    command = lambda: self.N_out(self.N))
        self.N_out_btn.grid(column = 1, row = 12, sticky = W)
        self.M_out_btn = ttk.Button(self, text = "Save-xyz", state = DISABLED, 
                                    command = lambda: self.N_out(self.M))
        self.M_out_btn.grid(column = 2, row = 12, sticky = W)
        self.plot_btn = ttk.Button(self,  text = "Plot", state = DISABLED,
                                   command = lambda: self.plot(self.C, self.N))
        self.plot_btn.grid(column = 3, row = 12, sticky = (N, W))

        # Default values for the interpolated grid
        self.x01_ = StringVar()
        self.y01_ = StringVar()
        self.x01_.set(0.600001) #default
        self.y01_.set(0.600001) #default

        self.dx1_ = DoubleVar() # from file
        self.dy1_ = DoubleVar() # from file

        self.y02_ = StringVar()
        self.x02_ = StringVar()
        self.x02_.set(0.60) #default
        self.y02_.set(0.60) #default

        self.dx2_ = StringVar()
        self.dy2_ = StringVar()
        self.dx2_.set(0.5) #default
        self.dy2_.set(0.5) #default

        self.np_x2_ = StringVar()
        self.np_y2_ = StringVar()
        self.np_x2_.set(13) #default
        self.np_y2_.set(13) #default
  
        self.grid2_lbl = ttk.Label(self, text = "Interpolated Grid")
        self.grid2_lbl.grid(column = 1, columnspan = 4, row = 22, 
                            sticky = (W,E), pady = (5,2))        
  
        x01_entry = ttk.Entry(self, width = 6, textvariable = self.x01_)
        x01_entry.grid(column = 2, row = 17, sticky = W)      
        x01_lbl = ttk.Label(self, text = 'x0')
        x01_lbl.grid(column = 1, row = 17, sticky = (N, W), padx = (40, 0))

        y01_entry = ttk.Entry(self, width = 6, textvariable = self.y01_)
        y01_entry.grid(column = 4, row = 17, sticky = W)
        y01_lbl = ttk.Label(self, text = 'y0')
        y01_lbl.grid(column = 3, row = 17, sticky = (N, W), padx = (40, 0))
        
        dx1_entry = ttk.Entry(self, width = 6, textvariable = self.dx1_)
        dx1_entry.grid(column = 2, row = 18, sticky = W)
        dx1_lbl = ttk.Label(self, text = 'dx1')
        dx1_lbl.grid(column = 1, row = 18, sticky = (N, W), padx = (40, 0))

        dy1_entry = ttk.Entry(self, width = 6, textvariable = self.dy1_)
        dy1_entry.grid(column = 4, row = 18, sticky = W)
        dy1_lbl = ttk.Label(self, text = 'dy1')
        dy1_lbl.grid(column = 3, row = 18, sticky = (N, W), padx = (40, 0))
        
        x02_entry = ttk.Entry(self, width = 6, textvariable = self.x02_)
        x02_entry.grid(column = 2, row = 23, sticky = W)
        x02_lbl = ttk.Label(self, text = 'x0 (Int)')
        x02_lbl.grid(column = 1, row = 23, sticky = (N, W), padx = (40, 0))

        y02_entry = ttk.Entry(self, width = 6, textvariable = self.y02_)
        y02_entry.grid(column = 4, row = 23, sticky = W)
        y02_lbl = ttk.Label(self, text = 'y0 (Int)')
        y02_lbl.grid(column = 3, row = 23, sticky = (N, W), padx =(40, 0))

        dx2_entry = ttk.Entry(self, width=6, textvariable=self.dx2_)
        dx2_entry.grid(column = 2, row = 24, sticky = W)
        dx2_lbl = ttk.Label(self, text = 'dx2')
        dx2_lbl.grid(column = 1, row = 24, sticky = (N, W), padx = (40, 0))

        dy2_entry = ttk.Entry(self, width = 6, textvariable = self.dy2_)
        dy2_entry.grid(column = 4, row = 24, sticky = W)
        dy2_lbl = ttk.Label(self, text = 'dy2')
        dy2_lbl.grid(column = 3, row = 24, sticky = (N, W), padx = (40, 0))

        np_x2_entry = ttk.Entry(self, width = 6, textvariable = self.np_x2_)
        np_x2_entry.grid(column = 2, row = 25, sticky = W)
        np_x2_lbl = ttk.Label(self, text = 'np_x2')
        np_x2_lbl.grid(column = 1, row = 25, sticky = (N, W), padx = (40, 0))

        np_y2_entry = ttk.Entry(self, width = 6, textvariable = self.np_y2_)
        np_y2_entry.grid(column = 4, row = 25, sticky = W)
        np_y2_lbl = ttk.Label(self, text = 'np_y2')
        np_y2_lbl.grid(column = 3, row = 25, sticky = (N, W), padx = (40, 0))
        
        self.tipwf = DoubleVar()
        self.tipwf.set(4220) 
        self.tip = IntVar()
        self.tip_btn = Checkbutton(self,  text = "Tip correction (mV)", 
                                   variable = self.tip, state = DISABLED, 
                                   command = lambda: 
                                       self.workfunction(self.C))
        self.tip_btn.grid(column = 3, columnspan = 2,row = 26, 
                          sticky = (N, W), pady = (15, 0))
        
        self.tip_entry = ttk.Entry(self, width = 6, textvariable = self.tipwf)
        self.tip_entry.grid(column = 2, row = 26, sticky = E, pady = (15, 0))
        
        minmax_lbl = ttk.Label(self, text = 'Min and Max x/y') 
        minmax_lbl.grid(column = 2, columnspan = 3, row = 19, sticky = (N, W))
        
        self.xi = IntVar()    
        self.xi.set(0)
        self.xf= IntVar()    
        self.xf.set(0)                
        xdata_lbl = ttk.Label(self, text = 'x-data_subset') 
        xdata_lbl.grid(column = 1, row = 20, sticky = (N, W))
        self.xidata_box = Spinbox(self, format = '%.0f', from_ = 0, 
                                  to = self.xf.get(), width = 5, increment = 1) 
        self.xidata_box.grid(column = 2, row = 20, sticky = (N, W), 
                             padx = (15, 0))
        self.xidata_box.insert(1, 0)        
        self.xfdata_box = Spinbox(self, format = '%.0f', from_ = 0, 
                                  to = self.xf.get(), width = 5, increment = 1) 
        self.xfdata_box.grid(column = 3, row = 20, sticky = (N, W))
        self.xfdata_box.insert(1, 0)
        
        self.yi = IntVar()    
        self.yi.set(0)
        self.yf = IntVar()    
        self.yf.set(0)      
        ydata_lbl = ttk.Label(self, text = 'y-data_subset') 
        ydata_lbl.grid(column = 1, row = 21, sticky = (N, W))
        self.yidata_box = Spinbox(self, format = '%.0f', from_ = 0, 
                                  to = self.yf.get(), width = 5, increment = 1) 
        self.yidata_box.grid(column = 2, row = 21, sticky = (N, W), 
                             padx = (15,0))
        self.yidata_box.insert(1,0)
        self.yfdata_box = Spinbox(self, format = '%.0f', from_ = 0, 
                                  to = self.yf.get(), width = 5, increment = 1)
        self.yfdata_box.grid(column = 3, row = 21, sticky = (N, W))
        self.yfdata_box.insert(1,0)        
    
    def enable_btn(self):
        self.graph1_btn.config(state = 'ENABLED')
        self.int_btn.config(state = 'ENABLED')
        self.trans_btn.config(state = 'normal')
        self.Refl_x_btn.config(state = 'normal')
        self.N_out_btn.config(state = 'ENABLED')
        self.M_out_btn.config(state = 'ENABLED')
        self.tip_btn.config(state = 'normal')
        
    def Input_C(self):
        """"Loads tab-delimited data from the corresponding text-box"""
        A = self.input_entry.get()
        self.info1_Entry.delete(1.0, END)
        
        B = []
        reader = csv.reader(A.split('\n'), delimiter = '\t')
        for row in reader:
            row = row[:]
            B.append(row)
        B = B[0:len(B)-1][0:len(B[0])]   
        for i in range(0, len(B)):
            for j in range(0, len(B[i])):
                B[i][j] = float(B[i][j])
        B = asarray(B)        
        print(B)
        self.dx1_.set(0.18)
        self.dy1_.set(0.18)
        self.NP1.set(len(B[0]) * len(B) - 1)
        self.n_box.config(to = self.NP1.get())
        
        self.xf.set(len(B[0]))
        self.xfdata_box.delete(0,"end")      
        self.xfdata_box.insert(1, self.xf.get())
        self.xidata_box.config(to = self.xf.get())
        self.xfdata_box.config(to = self.xf.get())
                
        self.yf.set(len(B))
        self.yfdata_box.delete(0,"end") 
        self.yfdata_box.insert(1, self.yf.get())
        self.yidata_box.config(to = self.yf.get())
        self.yfdata_box.config(to = self.yf.get())
        
        self.C = B
        self.enable_btn()
        self.graph1(B, 1)
        return self.C
    
    def open1(self):
        """"Loads data from file and displays footer information"""
        fileName = filedialog.askopenfilename(title = "Open CPD file")
        try:
            f = open(fileName, 'r')
        except IOError:
            pass

        i = f.readline().strip()
        if i == "Work Function data":
            pass
        else:
            i = messagebox.askyesno(title = "Open file?", 
                                    message = "CPD file header not present!\n"
                                    "Proceed to open file?", default = 'no')
            if i == 1:
                pass
            else:
                pass
        self.path1.set(fileName)
        g = f.read()
        f.close()
        
        # Only CPD data is loaded but there are other measurements 
        Tracking = g.find('Tracking')
        Grad = g.find('Grad')
        # Time = g.find('Time')
        # User = g.find('User')
        Footer = g.find('/') - 2 
        # Footer = Footer - 2

        CPD_data = g[0:Tracking]
        # Tracking_data = g[Tracking:Grad]
        # Grad_data = g[Grad:Time]
        # Time_data = g[Time:User]
        # User_data = g[User:Footer]
        Footer_data = g[Footer:]

        a = Footer_data.split("\n")
        Date = a[0]
        Grad = [int(i) for i in a[-10].split() if i.isdigit()]
        DAS_Averages = [int(i) for i in a[-18].split() if i.isdigit()]
        WF_Averages = [int(i) for i in a[-14].split() if i.isdigit()]
        Xpoints = [int(i) for i in a[-5].split() if i.isdigit()]
        Ypoints = [int(i) for i in a[-4].split() if i.isdigit()]
        Xsteps = [int(i) for i in a[-3].split() if i.isdigit()]
        Ysteps = [int(i) for i in a[-2].split() if i.isdigit()]
        Xstepscm = Xsteps[0] * 0.0015 #xyscale
        Ystepscm = Ysteps[0] * 0.0015
        
        # C is the array used for the raw data
        C = []
        reader = csv.reader(CPD_data.split('\n'), delimiter=',')
        for row in reader:
            row = row[:-1]
            C.append(row)
        C = C[0:len(C)-1][0:len(C[0])]   
        for i in range(0, len(C)):
            for j in range(0, len(C[i])):
                C[i][j] = float(C[i][j])
        C = asarray(C)
        
        # Clear the textbox before putting in data for new file
        self.info1_Entry.delete(1.0, END) 
        info_0 = "%s\n%s\n" % (fileName, Date)
        info_1 = ("Xpoints = %s\nYpoints = %s\nXsteps = %s (%s cm)\n"
                  "Ysteps= %s (%s cm)" \
            % (Xpoints[0], Ypoints[0], Xsteps[0], Xstepscm, Ysteps[0], 
               Ystepscm))
        info_2 = "\nDAS_Aves = %s\nWF_Averages = %s\nGrad = %s"\
            % (DAS_Averages[0],WF_Averages[0],Grad[0])
        self.info1_Entry.insert(END,info_0)
        self.info1_Entry.insert(END,info_1)
        self.info1_Entry.insert(END,info_2)
        if self.Refl_x.get() == 1:
            info_3 = "\nReflect x = %s" %(self.Refl_x.get())
            self.info1_Entry.insert(END,info_3)

        if self.transpose_chk.get() == 1:
            info_4 = "\nTranspose = %s" %(self.transpose_chk.get())
            self.info1_Entry.insert(END,info_4)
        
        self.enable_btn()
        
        # To prevent divide by 0 issue 1E-8 is added 
        if Xpoints[0] == len(C[0]):
            self.dx1_.set(Xstepscm + 1E-8) 
            self.dy1_.set(Ystepscm + 1E-8)           
        else:
            self.dy1_.set(Xstepscm + 1E-8)
            self.dx1_.set(Ystepscm + 1E-8)
        
        self.NP1.set(len(C[0])*len(C) - 1)
        self.n_box.config(to = self.NP1.get())
        
        self.xf.set(len(C[0]))
        self.xfdata_box.delete(0,"end")      
        self.xfdata_box.insert(1, self.xf.get())
        self.xidata_box.config(to = self.xf.get())
        self.xfdata_box.config(to = self.xf.get())
                
        self.yf.set(len(C))
        self.yfdata_box.delete(0, "end") 
        self.yfdata_box.insert(1, self.yf.get())
        self.yidata_box.config(to = self.yf.get())
        self.yfdata_box.config(to = self.yf.get())
        
        self.tip.set(0) 
        
        self.C = C
        self.Cr = C
        return self.C, self.Cr
    
    def graph1(self, data, row_a):
        """Graphs the raw CPD data"""
        plt.close()
        plt.figure(figsize = (3.5,3.5), dpi = 100, frameon=False)
        im = plt.imshow(data, cmap = 'jet')
        plt.colorbar(im, orientation = 'vertical',fraction = 0.046, 
                     pad = 0.05)

        plt.tight_layout()
        plt.savefig(path.join(tmpdir, im_name1), dpi = 100,
                    bbox_inches = 'tight')

        c = Canvas(self, width = 350, height = 350, bg = 'white')
        c.grid(row = row_a, column = 5, columnspan = 2, rowspan = 50,
               sticky = (W,N), padx = 50) 
        c.background = PhotoImage(file = path.join(tmpdir, im_name1))
        c.create_image(1, 20, image = c.background, anchor=NW)
    
    def interpolation(self, C):
        """Performs inverse distance interpolation"""
        temp=[]
        for j in range(int(self.yidata_box.get()),
                       int(self.yfdata_box.get())):
            temp.append([])
            for i in range(int(self.xidata_box.get()),
                           int(self.xfdata_box.get())):
                temp[j - int(self.yidata_box.get())].append(C[j,i])
        # Sub-array of raw data defined by the min and max x/y value
        C = temp
        C = asarray(C)
        
        n = self.n_box.get()
        n = int(n)

        p = self.p_box.get()
        p = int(p)

        x01 = float(self.x01_.get())
        y01 = float(self.y01_.get()) 
        
        dx1 = float(self.dx1_.get())
        dy1 = float(self.dy1_.get())

        y02 = float(self.y02_.get())
        x02 =  float(self.x02_.get())

        dx2 = float(self.dx2_.get())
        dy2 = float(self.dy2_.get())

        np_x2 = int(self.np_x2_.get())
        np_y2 = int(self.np_y2_.get())

        if self.transpose_chk.get() == 1:
            C = transpose(C)
        
        # Define the dimensions of the raw data grid
        np_x1 = len(C[0])       
        np_y1 = len(C)               
         
        NP1 = np_x1 * np_y1

        x1 = ones(np_x1)
        for i in range(0, np_x1):
            x1[i] = x01 + dx1*i

        y1 = ones(np_y1)
        for i in range(0, np_y1):
            y1[i] = y01 + dy1*i

        if self.Refl_x.get() == 1:
            x1_ = ones(np_x1)
            for i in range(0, len(x1)):
                x1_[i] = x1[-1 - i]
            x1 = x1_

        data_1d = C.flatten()

        # Define the dimensions of interpolated data grid
        NP2 = np_x2 * np_y2

        x2 = ones(np_x2)
        for i in range(0, np_x2):
            x2[i] = x02 + dx2*i

        y2 = ones(np_y2)
        for i in range(0, np_y2):
            y2[i] = y02 + dy2*i
        
        # Intperolation calculation
        k2 = 0
        int_data = ones(NP2)
        dist=[0] * (NP1)
        index=[0] * (NP1)
        for j2 in range (0,np_y2):
            for i2 in range(0,np_x2):
                dist_2d = [] 
                dist_1d = [] 
                for j in range (0,np_y1):
                    dist_2d.append([])
                    for i in range(0,np_x1):
                        dist_2d[j].append(
                            sqrt((x2[i2] - x1[i])**2 + (y2[j2] - y1[j])**2))
                        dist_1d.append(
                            sqrt((x2[i2] - x1[i])**2 + (y2[j2] - y1[j])**2))

                index = sorted(range(len(dist_1d)), key = lambda x: dist_1d[x])
                dist = sorted(dist_1d)
                up = min(n,NP1)
                temp = 0
                lamb = ones(up)
                if dist[0]> 0:
                    for i in range(0,up):
                        lamb[i] = (1/(dist[i]**p))\
                            /sum(ones(up)/power_(dist[0:up],p))
                        temp = temp + data_1d[index[i]]*lamb[i]
                else:
                    temp = data_1d[index[i]]            
                int_data[k2] = temp
                k2 = k2 + 1            

        M1 = []
        N = [] 
        k2 = 0
        for j2 in range (0,np_y2):
            N.append([])
            for i2 in range(0,np_x2):    
                M1.append([x02 + i2*dx2, y02 + j2*dy2, int_data[k2]])
                N[j2].append(int_data[k2])
                k2 = k2 + 1
        M1 = asarray(M1)
        N = asarray(N)

        # Generate a xyz equivalent array
        if np_x2 == np_y2:
            Mt=[]
            for j2 in range (0,np_y2):
                for i2 in range(0,np_x2):    
                    Mt.append([y02 + j2*dy2, x02 + i2*dx2, N[i2,j2]]) 
                    k2 = k2 + 1
            M = append(Mt, M1, axis=1)
        else:
            M = M1
        
        self.N = N
        self.M1 = M1
        self.M = M
        
        self.graph1(self.N, 8)
        
        self.plot_btn.config(state = 'ENABLED')
        return (self.N, self.M1, self.M)
       
    def check_tr(self, C):
        self.info1_Entry.delete(10.0, END)
        if self.Refl_x.get() == 1:
            info_3 = "\nReflect x = %s" %(self.Refl_x.get())
            self.info1_Entry.insert(END, info_3)

        if self.transpose_chk.get() == 1:
            info_4 = "\nTranspose = %s" %(self.transpose_chk.get())
            self.info1_Entry.insert(END, info_4)
        self.interpolation(C)
    
    def plot(self, C, N):
        """Produces images of the raw and interpolated data plots"""
        plt.close('all')
        fig_ext = plt.figure(dpi=200)
        ax1 = fig_ext.add_subplot(121)
        ax1_im = plt.imshow(C, vmin = amin([amin(C),amin(N)]),
                     vmax = amax([amax(C), amax(N)]), cmap = 'jet')
        
        cbar = plt.colorbar(ax1_im, orientation = 'vertical', fraction = 0.046, 
                          pad = 0.05)# aspect='auto') 
        if self.tip.get() == 1:
            cbar.set_label('Work Function (meV)')
        else:
            cbar.set_label('CPD (mV)')
        plt.xlabel('x-axis')
        plt.ylabel('y-axis')
        
        ax2 = fig_ext.add_subplot(122)
        ax2_im = plt.imshow(N, vmin = amin([amin(C), amin(N)]),
                     vmax = amax([amax(C), amax(N)]), cmap = 'jet')
        cbar=plt.colorbar(ax2_im, orientation = 'vertical', fraction = 0.046, 
                          pad = 0.05)#, aspect='auto')
        if self.tip.get() == 1:
            cbar.set_label('Work Function (meV)')
        else:
            cbar.set_label('CPD (mV)')
        plt.xlabel('x-axis')
        plt.ylabel('y-axis')
        
        plt.tight_layout()
        fig_ext = plt.gcf()
        plt.show()

    def N_out(self, N):
        """Exports the data as a csv file via dialogue box"""            
        file_opt = options = {}
        options['filetypes'] = [('all files', '.*'), 
                                ('CSV (Comma Delimited)', '.csv')]
        options['initialfile'] = '.csv'
        options['title'] = 'Save Array as CSV'

        try:
            fo1 = filedialog.asksaveasfile(mode = 'wb', **file_opt)
            savetxt(fo1, N, delimiter = ',', 
                    footer = self.info1_Entry.get(1.0,END))
            fo1.close()
        except IOError:
            messagebox.showinfo(title = "Output file",
                                message = "Output file issue")

    def workfunction(self, C):      
        """Applies linear offset to CPD, useful for correcing 
        for the tip's workfunction and converting CPD into 
        workfunction"""
        if self.tip.get() == 1:
            self.C = C + float(self.tipwf.get())
            self.tip_entry.config(state = 'disabled')     
        else:
            self.C = (C - float(self.tipwf.get()))
            self.tip_entry.config(state = "normal")      
        self.graph1(self.C, 1)
        return self.C

root = Tk()
root.geometry("850x820")
root.title("CPD Visualisation and Interpolation Module")

app = Application(master=root)
app.mainloop()