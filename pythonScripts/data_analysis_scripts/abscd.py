'''A set of classes and functions for reading in, storing, and analyzing Abs/CD Data.'''

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import plotly
import plotly.express as px
import plotly.graph_objs as go
import os

class AbsCD_Data:
    '''Class to hold Abs/CD data from the J-1700.'''
    def __init__(self, data = None, source = None,
                 xlabels = [], ylabels = [], sample = None, 
                 path:'path to folder with info/data csv files' = None) -> None: # type: ignore
        self.data = data
        self.source = source
        self.xlabels = xlabels
        self.ylabels = ylabels
        self.sample = sample
        if path is not None:
            self.parse_data(path)

    def copy(self):
        '''Copy instance of AbsCD_Data Class'''
        new = AbsCD_Data()
        new.data = self.data.copy()
        new.source = self.source
        new.xlabels = self.xlabels.copy()
        new.ylabels = self.ylabels.copy()
        new.sample = self.sample
        return new
    
    def parse_data(self, path):
        #parse through folder of path to get relevant files
       
        #parse files to find log and data
        data_files=[]
        for root,dirs,file in os.walk(top=path):
            if root==path:
                for f in file:
                    if f.endswith('.csv'):
                        if 'log' in f or 'info' in f:
                            log_file=f
                        else:
                            data_files.append(f)
        #use pandas to read in info file
        log = pd.read_csv(path+log_file)
        #use pandas to read in and aggregate data files
        #parse data from each data file
        for i,row in log.iterrows():
            if log.at[i,'File'] in data_files:
                with open(path+log.at[i,'File'], 'r') as f:
                    for line in f:
                        if 'UNITS' in line and 'X' in line:
                            xstr=line.split(',')[-1].replace('\n','')
                            if xstr not in self.xlabels:
                                self.xlabels.append(xstr)
                                log[xstr] = pd.Series(dtype='object')
                        elif 'UNITS' in line and 'Y' in line:
                            ystr=line.split(',')[-1].replace('\n','')
                            if ystr not in self.ylabels:
                                self.ylabels.append(ystr)
                                log[ystr ] = pd.Series(dtype='object')
                        elif 'NPOINTS' in line:
                            npts=int(line.split(',')[-1])
                
                #print(path+log.at[i,'File'])
                data=pd.read_csv(path+log.at[i,'File'], skiprows=22, header=None,nrows=npts)
            
                #write data to xs
                log.at[i, self.xlabels[0]]=pd.Series(data[0]).values
                #write data to ys
                yidx=1
                for lab in self.ylabels:
                    log.at[i, lab]=pd.Series(data[yidx]).values
                    yidx+=1
                
        #Save data to class object
        self.data = log.copy()
        log = None
        if self.source is None:
            self.source = 'Abs/CD Data From ' + path

    def quick_plot(self, x=None, y=None, y_range=None, x_range=None):
    #Makes interactive plotly plots from the data given specific x and y (str or list) parameters as the df column names
        df = self.data
        matplotlib.rcParams['figure.figsize'] = [10, 5]
        #set up y data stucture
        if y==None:
            y=self.ylabels
        if type(y) is str:
            y=[y]
        #set up x data structure
        if x==None:
            x = ('Wavenums' if 'Wavenums' in self.xlabels else 'NANOMETERS')

        #print plot statement for all data
        for yAx in y:
            #Create plotly plot
            fig = go.Figure(layout=go.Layout(
                    width =1000, height=500,title=self.sample, margin=go.layout.Margin(b=50,t=50, l=20)))
            
            fig.update_xaxes(title_text=x)
            fig.update_yaxes(title_text=yAx)

            for i in list(self.data.index):
                fig.add_trace(go.Scatter(x=df.at[i,x],y=df.at[i,yAx],name=str(self.data.at[i,'Id']), mode='lines', line=go.scatter.Line(width=2)))
            
            
            fig.add_hline(y=0)
            if y_range is not None:
                fig.update_layout(yaxis_range=y_range)
            if x_range is not None:
                fig.update_layout(xaxis_range=x_range)
            fig.show()

    def add_wavenums(self):
    #Add an x value of wavenumbers by converting nanometers
        if 'NANOMETERS' in self.xlabels and 'Wavenums' not in self.xlabels:
            self.data['Wavenums'] = np.power(self.data['NANOMETERS'],-1)*10000000
            self.xlabels.append('Wavenums')
        print(self.xlabels)


    def add_deps(self, conc, conc_units='M'):
        '''Add a y value of Delta Epsilon (1/M*cm) by converting from mdeg.
        conc - takes numerical value or name of dataframe column in self.data'''
        if type(conc) is str:
            conc = self.data[conc].to_numpy()
        else:
            try:
                conc = float(conc)
            except:
                print('Please give conc as a number or df column name.')
                return False

        #Convert to M
        if conc_units == 'mM':
            concM = conc / 1000  # Convert mM to M
        elif conc_units == 'uM':
            concM = conc / 1000000  # Convert µM to M
        elif conc_units == 'M':
            concM = conc

        if 'CD/DC [mdeg]' in self.ylabels and 'deps' not in self.ylabels:
            self.data['deps'] = self.data['CD/DC [mdeg]']/(concM * 0.3 * 32980)
            self.ylabels.append('deps')
        print(self.ylabels)

    def add_eps(self, conc, conc_units='M', path_length:'cm'=1): # type: ignore
        '''Add a y value of Epsilon (1/M*cm) by converting from Abs.
        conc - takes numerical value or name of dataframe column in self.data'''
        if type(conc) is str:
            conc = self.data[conc].to_numpy()
        else:
            try:
                conc = float(conc)
            except:
                print('Please give conc as a number or df column name.')
                return False

        #Convert to M
        if conc_units == 'mM':
            concM = conc / 1000  # Convert mM to M
        elif conc_units == 'uM':
            concM = conc / 1000000  # Convert µM to M
        elif conc_units == 'M':
            concM = conc

        if 'ABSORBANCE' in self.ylabels and 'eps' not in self.ylabels:
            self.data['eps'] = self.data['ABSORBANCE']/(concM * path_length)
            self.ylabels.append('eps')
        print(self.ylabels)
            

    
        

    