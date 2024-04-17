import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from IPython.display import Markdown as md
from matplotlib.ticker import FormatStrFormatter
import os
from scipy.optimize import least_squares,minimize,LinearConstraint,NonlinearConstraint,Bounds

class VTVH_Data:
    def __init__(self, path=None):
        self.data = pd.DataFrame()
        self.xlabels = []
        self.ylabels = []
        self.source = None
        self.sample = None
        if path is not None:
            self.parse_vtvh(path)
    
    def copy(self):
        new=VTVH_Data()
        new.data = self.data.copy()
        new.xlabels=self.xlabels.copy()
        new.ylabels=self.ylabels.copy()
        new.source=self.source
        new.sample=self.sample
        return new
        
    #parse data when given path
    def parse_vtvh(self, path, scan_set=0):
        #parse through folder to get relevant files
       
        #parse files to find log and data
        data_files=[]
        for root,dirs,file in os.walk(top=path):
            if root==path:
                for f in file:
                    if f.endswith('.csv'):
                        if 'vtvh_log' in f:
                            log_file=f
                        else:
                            data_files.append(f)
        log=pd.read_csv(path+log_file)
        
        temp_set=scan_set
        idx=0
        num_scans=0
        #parse through data in log file to Only Count Number of scans
        for i,row in log.iterrows():
            #if have reached a new set of scans
            if int(log.at[i,'Scan_Num'])-1 < idx:
                temp_set-=1
            
            idx=int(log.at[i,'Scan_Num'])-1
            #only scan from the set requested
            if temp_set==0:
                if log.at[i,'File'] in data_files:
                    idx=int(log.at[i,'Scan_Num'])-1
                    num_scans+=1 #ADD 1 to number of scans
 
        #intitialize dict that will become raw data df
        vtvh_data={}
        vtvh_data['Id']=['' for n in range(num_scans)]
        vtvh_data['ScanNum']=np.zeros(num_scans, dtype=int)
        vtvh_data['FieldSet']=np.zeros(num_scans)
        vtvh_data['TempSet']=np.zeros(num_scans)
        vtvh_data['FieldAvg']=np.zeros(num_scans)
        vtvh_data['FieldStdev']=np.zeros(num_scans)
        vtvh_data['TempAvg']=np.zeros(num_scans)
        vtvh_data['TempStdev']=np.zeros(num_scans)
        vtvh_data['NumChecks']=np.zeros(num_scans)
        vtvh_data['Path']=['' for n in range(num_scans)]
        
        idx=0
        #parse through data in log file to assign information to data from scans
        for i,row in log.iterrows():
            #if have reached a new set of scans
            if int(log.at[i,'Scan_Num'])-1 < idx:
                scan_set-=1
            idx=int(log.at[i,'Scan_Num'])-1
            #only scan from the set requested
            if scan_set==0:
                #parse data from vtvh log file
                if log.at[i,'File'] in data_files:
                    idx=int(log.at[i,'Scan_Num'])-1
                    id_str='{s}_{t}K_{h}T'.format(s=log.at[i,'Scan_Num'],t=log.at[i,'SampleTemp_SetPt(K)'],h=log.at[i,'Field_SetPoint(T)'])
                    vtvh_data['Id'][idx]=id_str
                    vtvh_data['FieldSet'][idx]=float(log.at[i,'Field_SetPoint(T)'])
                    vtvh_data['TempSet'][idx]=float(log.at[i,'SampleTemp_SetPt(K)'])
                    vtvh_data['ScanNum'][idx]=float(log.at[i,'Scan_Num'])
                    vtvh_data['Path'][idx]=path+log.at[i,'File']
                    if 'Avg_Sample_Temp(K)' in log.columns:
                        vtvh_data['FieldAvg'][idx]=float(log.at[i,'Avg_Magnet_Field(T)'])
                        vtvh_data['TempAvg'][idx]=float(log.at[i,'Avg_Sample_Temp(K)'])
                        vtvh_data['FieldStdev'][idx]=float(log.at[i,'StdDev_Magnet_Field(T)'])
                        vtvh_data['TempStdev'][idx]=float(log.at[i,'StdDev_Sample_Temp(K)'])
                        vtvh_data['NumChecks'][idx]=float(log.at[i,'Num_Temp_Checks'])
                    #parse data from each data file
                    with open(path+log.at[i,'File'], 'r') as f:
                        for line in f:
                            if 'UNITS' in line and 'X' in line:
                                xstr=line.split(',')[-1].replace('\n','')
                                if xstr not in self.xlabels:
                                    self.xlabels.append(xstr)
                                    vtvh_data[xstr]=np.zeros(num_scans, dtype=pd.Series)
                            elif 'UNITS' in line and 'Y' in line:
                                ystr=line.split(',')[-1].replace('\n','')
                                if ystr not in self.ylabels:
                                    self.ylabels.append(ystr)
                                    vtvh_data[ystr]=np.zeros(num_scans, dtype=pd.Series)
                            elif 'NPOINTS' in line:
                                npts=int(line.split(',')[-1])
                    data=pd.read_csv(path+log.at[i,'File'], skiprows=21, header=None,nrows=npts)
                    #write data to xs
                    vtvh_data[self.xlabels[0]][idx]=pd.Series(data[0])
                    #write data to ys
                    yidx=1
                    for lab in self.ylabels:
                        vtvh_data[lab][idx]=pd.Series(data[yidx])
                        yidx+=1
            
        #Save data to class object
        self.data = pd.DataFrame(data=vtvh_data)
        self.source = 'Data From ' + path
        
    #Operation to Subtract Zeroes off of data - returns a new VTVH_Data object    
    def subtract_zeros(self, y='CD/DC [mdeg]', zero=-1, new_id=False):
        subtracted = self.data.loc[self.data['FieldSet']!=0].copy()
        #if want to subtract from 0T right before
        if zero == -1:
            for i,row in subtracted.iterrows():
                z=self.data.loc[(self.data['FieldSet']==0) & (self.data['ScanNum']<row['ScanNum'])].sort_values('ScanNum',ignore_index=True).loc[0]
                subtracted.at[i,y]=np.subtract(row[y],z[y])
                if new_id:
                    subtracted.at[i,'Id']=row['Id'] + ' - ' + z['Id']
        elif zero == 'after':
            for i,row in subtracted.iterrows():
                z=self.data.loc[(self.data['FieldSet']==0) & (self.data['ScanNum']>row['ScanNum'])].sort_values('ScanNum',ascending=True,ignore_index=True).loc[0]
                subtracted.at[i,y]=np.subtract(row[y],z[y])
                if new_id:
                    subtracted.at[i,'Id']=row['Id'] + ' - ' + z['Id']
        else: 
            z=self.data.loc[self.data['Id']==zero]
            for i,row in subtracted.iterrows():
                subtracted.at[i,y]=np.subtract(row[y],z[y])
        sub_data = self.copy()
        sub_data.data = subtracted
        sub_data.source = self.source + ' subtracted off ' + ('each\'s previous 0T' if zero==-1 else zero)
        return sub_data
    
    #Do 0.5*(+T - -T) as a replacement for subtracting off zeroes
    def half_subd_fields(self, field=None):
        y='CD/DC [mdeg]'
        if field is None:
            subtracted = self.data.loc[self.data['FieldSet']>0].copy()
        elif type(field) is float or type(field) is int:
            subtracted = self.data.loc[self.data['FieldSet']==field].copy()
        elif type(field) is list:
            subtracted = self.data.loc[self.data['FieldSet'] in field].copy()
          
        #Find negative field of same magnitude
        for i,row in subtracted.iterrows():
            oppo=self.data.loc[(self.data['FieldSet']==(-1)*row['FieldSet']) & (self.data['TempSet']==row['TempSet'])].reset_index().loc[0]
            subtracted.at[i,y]=0.5*np.subtract(row[y],oppo[y])
            subtracted.at[i,'Id']='0.5*(' + row['Id'] + ' - ' + oppo['Id'] + ')'

        sub_data = self.copy()
        sub_data.data = subtracted
        sub_data.source = self.source + ' checked mirroring by adding opposite signed field (+T + -T) for CD only'
        return sub_data
    
    #Do +T + -T (should be subtracted from 0 already) which should go to 0 to check for mirroring
    def check_mirroring(self, field=None):
        y='CD/DC [mdeg]'
        if field is None:
            subtracted = self.data.loc[self.data['FieldSet']>0].copy()
        elif type(field) is float or type(field) is int:
            subtracted = self.data.loc[self.data['FieldSet']==field].copy()
        elif type(field) is list:
            subtracted = self.data.loc[self.data['FieldSet'] in field].copy()
          
        #Find negative field of same magnitude
        for i,row in subtracted.iterrows():
            oppo=self.data.loc[(self.data['FieldSet']==(-1)*row['FieldSet']) & (self.data['TempSet']==row['TempSet'])].reset_index().loc[0]
            subtracted.at[i,y]=np.add(row[y],oppo[y])
            subtracted.at[i,'Id']=row['Id'] + ' + ' + oppo['Id']

        sub_data = self.copy()
        sub_data.data = subtracted
        sub_data.source = self.source + ' checked mirroring by adding opposite signed field (+T + -T) for CD only'
        return sub_data
    
    #Add an x value of wavenumbers by converting nanometers
    def add_wavenums(self):
        if 'NANOMETERS' in self.xlabels and 'Wavenums' not in self.xlabels:
            #cm_list = []
            #for i, row in self.data.iterrows():
            #    nm = row['NANOMETERS']
            #    cm_list.append(np.power(nm,-1)*10000000)
            #self.data['Wavenums'] = cm_list
            self.data['Wavenums'] = np.power(self.data['NANOMETERS'],-1)*10000000
            self.xlabels.append('Wavenums')
        print(self.xlabels)
        
    #Add a y value of Delta Epsilon by converting from mdeg
    def add_deps(self, concM):
        if 'CD/DC [mdeg]' in self.ylabels and 'deps' not in self.ylabels:
            self.data['deps'] = self.data['CD/DC [mdeg]']/(concM * 0.3 * 32980)
            self.ylabels.append('deps')
        print(self.ylabels)
    
    #prints out commands to make plots
    def plot_prep(self, x=None, y=None, field=None, temp=None, add_sample=False):
        print('df= VARNAME.data')
        print('matplotlib.rcParams[\'figure.figsize\'] = [10, 5]')
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
            for i in list(self.data.index):
                if (field == None and temp == None) or self.data.at[i,'FieldSet'] == field or self.data.at[i,'TempSet'] == temp:
                    print('plt.plot(df.at[{row},\'{xdat}\'],df.at[{row},\'{ydat}\'],label=\'{lab}\')'.format(row=i,xdat=x,ydat=yAx,lab=(self.data.at[i,'Id']+self.sample if add_sample else self.data.at[i,'Id'])))
                
            if yAx in ['CD/DC [mdeg]','ABSORBANCE']:
                print('plt.hlines(0,{minx},{maxx},\'k\')'.format(minx=int(min(self.data.at[i,x])), maxx=int(max(self.data.at[i,x]))))
            print('plt.xlim({minx},{maxx})'.format(minx=int(min(self.data.at[i,x])), maxx=int(max(self.data.at[i,x]))))
            print('plt.xlabel(\'{}\')'.format(x))
            print('plt.ylabel(\'{}\')'.format(yAx))
            print('plt.ticklabel_format(axis=\'x\', scilimits=(3,3), useMathText=True)')
            print('plt.legend()')
            print('plt.show()')
    
    #Just makes the plots based on default parameters in plot_prep
    def quick_plot(self, x=None, y=None, field=None, temp=None):
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
            for i in list(self.data.index):
                if (field == None or self.data.at[i,'FieldSet'] == field) and (temp == None or self.data.at[i,'TempSet'] == temp):
                    plt.plot(df.at[i,x],df.at[i,yAx],label=self.data.at[i,'Id'])
            
            if yAx in ['CD/DC [mdeg]','ABSORBANCE']:
                plt.hlines(0,int(min(self.data.at[i,x])),int(max(self.data.at[i,x])),'k')
            plt.xlim(int(min(self.data.at[i,x])),int(max(self.data.at[i,x])))
            plt.xlabel(x)
            plt.ylabel(yAx)
            plt.legend()
            plt.show()
                
        


class Multi_VTVHs:
    def __init__(self, vtvhs=[]):
        self.vtvhs=vtvhs
        self.cmaps=[plt.cm.Blues_r(np.linspace(0,1,10)), 
                    plt.cm.Reds_r(np.linspace(0,1,10)),
                    plt.cm.Greens_r(np.linspace(0,1,10)),
                    plt.cm.Oranges_r(np.linspace(0,1,10)),
                    plt.cm.Purples_r(np.linspace(0,1,10)),
                    plt.cm.Greys_r(np.linspace(0,1,10))]
        self.features=[]
        
    def copy(self):
        new_set=[]
        for run in self.vtvhs:
            new_set.append(run.copy())
        new_multi=Multi_VTVHs(new_set)
        new_multi.features=self.features.copy()
        return new_multi
    
    def update_cmaps(self, length):
        self.cmaps=[plt.cm.Blues_r(np.linspace(0,1,length+2)), 
                            plt.cm.Reds_r(np.linspace(0,1,length+2)),
                            plt.cm.Greens_r(np.linspace(0,1,length+2)),
                            plt.cm.Oranges_r(np.linspace(0,1,length+2)),
                            plt.cm.Purples_r(np.linspace(0,1,length+2)),
                            plt.cm.Greys_r(np.linspace(0,1,length+2))]
    #Operation to Subtract Zeroes off of data - returns a new VTVH_Data object    
    def subtract_zeros(self, y='CD/DC [mdeg]', zero=-1, new_id=False):
        new=self.copy()
        new.vtvhs=[]
        for run in self.vtvhs:
            new.vtvhs.append(run.subtract_zeros(y,zero,new_id))
        return new

    #Do 0.5*(+T - -T) as a replacement for subtracting off zeroes
    def half_subd_fields(self, field=None):
        new=self.copy()
        new.vtvhs=[]
        for run in self.vtvhs:
            new.vtvhs.append(run.half_subd_fields(field))
        return new

    #Do +T + -T (should be subtracted from 0 already) which should go to 0 to check for mirroring
    def check_mirroring(self, field=None):
        new=self.copy()
        new.vtvhs=[]
        for run in self.vtvhs:
            new.vtvhs.append(run.check_mirroring(field))
        return new

    #Add an x value of wavenumbers by converting nanometers
    def add_wavenums(self):
        new=self.copy()
        new.vtvhs=[]
        for run in self.vtvhs:
            new.vtvhs.append(run.add_wavenums())
        return new

    #prints out commands to make plots
    def plot_prep(self, x=None, y='CD/DC [mdeg]', field=None, temp=None, samp_str=None, full_legend=True):
        print('group= VARNAME')
        print('matplotlib.rcParams[\'figure.figsize\'] = [10, 5]')
        if not full_legend:
            print('full_leg=False')
        
        to_plot=[]
        for idx in range(len(self.vtvhs)):
            if (samp_str in self.vtvhs[idx].sample or samp_str is None) and (field in self.vtvhs[idx].data['FieldSet'].unique() or field is None) and (temp in self.vtvhs[idx].data['TempSet'].unique() or temp is None):
                to_plot.append(idx)
        color_idx=0
        for j in to_plot:
            df=self.vtvhs[j].data
            print('df= group.vtvhs[{}].data'.format(j))

            #set up x data structure
            if x==None:
                x = ('Wavenums' if 'Wavenums' in self.vtvhs[j].xlabels else 'NANOMETERS')
            
            k=0
            #print plot statement for all data
            for i in list(self.vtvhs[j].data.index):
                if (field == None or self.vtvhs[j].data.at[i,'FieldSet'] == field) and (temp == None or self.vtvhs[j].data.at[i,'TempSet'] == temp):
                    #create label
                    lab='\'{}T, {}K {}\''.format(df.at[i,'FieldSet'], df.at[i,'TempSet'], self.vtvhs[j].sample)
                    if k!=0 and not full_legend:
                        lab='\'{}T, {}K {}\' if full_leg else None'.format(df.at[i,'FieldSet'], df.at[i,'TempSet'], self.vtvhs[j].sample)
                    #actually print plot command
                    print('plt.plot(df.at[{row},\'{xdat}\'],df.at[{row},\'{ydat}\'],label={lab}, color=group.cmaps[{color_idx}][{k}])'.format(row=i,xdat=x,ydat=y,
                                                                                                        lab=lab,
                                                                                                        color_idx=color_idx, k=k))
                    #track how many lines are being plotted from this dataset
                    k+=1
            #track which dataset its on        
            color_idx+=1

            if y in ['CD/DC [mdeg]','ABSORBANCE']:
                print('plt.hlines(0,{minx},{maxx},\'k\')'.format(minx=int(min(self.vtvhs[j].data.at[i,x])), maxx=int(max(self.vtvhs[j].data.at[i,x]))))
            print('plt.xlim({minx},{maxx})'.format(minx=int(min(self.vtvhs[j].data.at[i,x])), maxx=int(max(self.vtvhs[j].data.at[i,x]))))
            print('plt.xlabel(\'{}\')'.format(x))
            print('plt.ylabel(\'{}\')'.format(y))
            print('plt.ticklabel_format(axis=\'x\', scilimits=(3,3), useMathText=True)')
            print('plt.legend()')
            print('plt.show()')
    #Just makes the plots based on default parameters in plot_prep
    def quick_plot(self, x=None, y='CD/DC [mdeg]', field=None, temp=None, samp_str=None, together=True, full_legend=True, cmap=True):
        matplotlib.rcParams['figure.figsize'] = [10, 5]
        to_plot=[]
        for run in self.vtvhs:
            if (samp_str in run.sample or samp_str is None) and (field in run.data['FieldSet'].unique() or field is None) and (temp in run.data['TempSet'].unique() or temp is None):
                to_plot.append(run)
        for j in range(len(to_plot)):
            df = to_plot[j].data
            #set up x data structure
            if x==None:
                x = ('Wavenums' if 'Wavenums' in to_plot[j].xlabels else 'NANOMETERS')

            #print plot statement for all data
            k=0
            #fix issue when cmaps are too small
            count=0
            for i in list(df.index):
                if (field == None or df.at[i,'FieldSet'] == field) and (temp == None or df.at[i,'TempSet'] == temp):
                    count+=1
            if count>10:
                self.update_cmaps(length=count)
            #plot all relevant data
            for i in list(df.index):
                if (field == None or df.at[i,'FieldSet'] == field) and (temp == None or df.at[i,'TempSet'] == temp):
                    if full_legend or k==0:
                        lab='{}T, {}K {}'.format(df.at[i,'FieldSet'], df.at[i,'TempSet'], to_plot[j].sample)
                    else:
                        lab=None
                    plt.plot(df.at[i,x],df.at[i,y],label=lab, color=(self.cmaps[j][k] if cmap else None))
                    k+=1
            
            #reset cmaps
            self.update_cmaps(length=8)
            
            if y in ['CD/DC [mdeg]','ABSORBANCE']:
                plt.hlines(0,int(min(df.at[i,x])),int(max(df.at[i,x])),'k')
            plt.xlim(int(min(df.at[i,x])),int(max(df.at[i,x])))
            plt.xlabel(x)
            plt.ylabel(y)
            plt.legend()
            if not together:
                plt.show()

'''Fitting functions not part of the class'''

#Functions for fitting
def gauss(x, center, width):
    return np.exp(-1/2*(x-center)**2/(width)**2)

#with xs and ys both as len(3) arrays for Abs then CD then MCD
def resid_abscdmcd(fitvars, xs, ys):
    #fit vars should be in format [energy1, energy2, ..., width1, w2, ..., scalarAbs1, sA2, ..., scalarCD1, sCD2, ..., sMCD, ..., energyMCD, ..., widthMCD]
    num_gauss = int(len(fitvars)/(3+2+2))
        
    total_resid=np.array([])
    #for each y given, calculate gaussians, fit, and add to list of residuals
    for j in range(len(ys)):
        x=xs[j]
        if j < 2:
            #get a list of the individual gaussian y values
            gauss_list = [fitvars[i+((2+j)*num_gauss)]*gauss(x, fitvars[i], fitvars[i+num_gauss]) for i in range(num_gauss)]
        elif j==2:
            #get a list of the individual gaussian y values for MCD
            gauss_list = [fitvars[i+((2+j)*num_gauss)]*gauss(x, fitvars[i+(5*num_gauss)], fitvars[i+(6*num_gauss)]) for i in range(num_gauss)]
        else:
            print('Error: Too Many Y Data Sets')
            break
        #calculate total for Abs with current params
        total_fit = np.sum(gauss_list, axis=0)
        #calculate total residual and add to the list
        total_resid = np.concatenate((total_resid, ys[j] - total_fit))
        
    return np.sum(np.power(total_resid,2))


def check_plot_abscdmcd(result, xs, ys):
    fitvars=result
    #fit vars should be in format [energy1, energy2, ..., width1, w2, ..., scalarAbs1, sA2, ..., scalarCD1, sCD2, ...]
    num_gauss = int(len(fitvars)/(3+2+2))
    #create dataframe for the results
    fit = pd.DataFrame()
    if len(xs)>10:
        fit['x'] = xs.copy()
    #for each y given, calculate gaussians, fit, and add to list of residuals
    for j in range(len(ys)):
        fit = pd.DataFrame()
        ylab = 'y'+str(j)
        #check to see if multiple x lists
        if len(xs) < 10: #probably wouldnt be fitting more than 10 at once? this would mean its just a list of x values
            x=xs[j].copy()
            fit['x_'+ylab] = xs[j].copy()
        else:
            x=xs
            
        if j < 2:
            #get a list of the individual gaussian y values
            gauss_list = [fitvars[i+((2+j)*num_gauss)]*gauss(x, fitvars[i], fitvars[i+num_gauss]) for i in range(num_gauss)]
        elif j==2:
            #get a list of the individual gaussian y values for MCD
            gauss_list = [fitvars[i+((2+j)*num_gauss)]*gauss(x, fitvars[i+(5*num_gauss)], fitvars[i+(6*num_gauss)]) for i in range(num_gauss)]
        else:
            print('Error: Too Many Y Data Sets')
            break
        #calculate total for Abs with current params
        total_fit = np.sum(gauss_list, axis=0)
        
        #Add expt, total fit, and each gaussian
        fit['expt_'+ylab] = ys[j].copy()
        fit['fit_'+ylab]= total_fit.copy()
        for k in range(len(gauss_list)):
            lab=ylab+'_g'+str(k)
            fit[lab] = gauss_list[k]
        if 'x' in fit.columns:
            fit.plot(x='x', y=[h for h in fit.columns if ylab in h])
        else:
            fit.plot(x=str('x_'+ylab), y=[h for h in fit.columns if (ylab in h and 'x_' not in h)])
    return fit

def fit_abs_cd_mcd(energies, widths, intens, xs, ys, low_bds=None, up_bds=None, same_x=True, de_mcd=300, dw_mcd=1000):
    #make sure all are floats 
    energies = [float(e) for e in energies]
    widths= [float(w) for w in widths]
    #intens = np.concatenate((intens))
    intens= [float(ints) for ints in intens]
    
    num_gauss = len(energies)
    #normalize ys, intensities, and bounds
    areas= np.zeros(len(ys))
    nys = []
    nintens = np.array(intens)
    #handle bounds setup
    if low_bds is not None:
        low_bds = [float(lb) for lb in low_bds]
        nl_bds = np.array(low_bds)
    else:
        nl_bds = -1*np.inf
    if up_bds is not None:
        up_bds = [float(ub) for ub in up_bds]
        nu_bds = np.array(up_bds)
        
    else:
        nu_bds = np.inf
    #iterate through ys, calc and store areas, and normalize
    for k in range(len(ys)):
        areas[k]=np.trapz(abs(ys[k]), x=xs[k])#xs if same_x else xs[k])
        nys.append(np.divide(ys[k], areas[k]))
        nintens[k*num_gauss:k*num_gauss+num_gauss] = np.divide(intens[k*num_gauss:k*num_gauss+num_gauss], areas[k])

        #normalize bounds for intensities
        if low_bds is not None:
            nl_bds[k*num_gauss+(2*num_gauss):k*num_gauss+(3*num_gauss)] = np.divide(low_bds[(2*num_gauss)+k*num_gauss:k*num_gauss+(3*num_gauss)], areas[k])
        if up_bds is not None:
            nu_bds[k*num_gauss+(2*num_gauss):k*num_gauss+(3*num_gauss)] = np.divide(up_bds[(2*num_gauss)+k*num_gauss:k*num_gauss+(3*num_gauss)], areas[k])    

    #Normalize the energies to see if that helps
    xmax = np.max([np.max(x) for x in xs])
    n_energies = np.divide(energies,xmax)
    n_widths = np.divide(widths, xmax)
    nxs = np.divide(xs, xmax)
    #normalize bounds for x
    if low_bds is not None:
        nl_bds[0:num_gauss*2] = np.divide(low_bds[0:num_gauss*2], xmax)
    if up_bds is not None:
        nu_bds[0:num_gauss*2] = np.divide(up_bds[0:num_gauss*2], xmax)    
    
    #create additional parameters for mcd
    e_mcd = n_energies
    w_mcd = np.multiply(n_widths,0.9)
    #give same bounds initally as the rest of the e and w
    if low_bds is not None:
        nl_bds = np.concatenate((nl_bds, nl_bds[0:num_gauss*2]))
    if up_bds is not None:
        nu_bds = np.concatenate((nu_bds, nu_bds[0:num_gauss*2]))

    #prepare input lists
    params = np.concatenate((n_energies, n_widths, nintens, e_mcd, w_mcd))
    bounds = Bounds(ub=nu_bds, lb=nl_bds)
    
    def mcdE_fxn(parm, num):
        ng = int(len(parm)/(4+3))
        e_diff = parm[num] - parm[ng*5 + num]
        return (de_mcd/xmax) - abs(e_diff)
    nlc_mcde = [{'type': 'ineq', 'fun': mcdE_fxn, 'args': (i)} for i in range(num_gauss)]
    
    def mcdW_fxn(parm, num):
        ng = int(len(parm)/(4+3))
        w_diff = parm[num+ng] - parm[ng*6 + num]
        return w_diff 
    nlc_mcdw=[{'type': 'ineq', 'fun': mcdW_fxn, 'args': (i)} for i in range(num_gauss)]
    
    consts = np.concatenate(([],[]))
    
    #Run fit
    fit = minimize(resid_abscdmcd, params, bounds=bounds, constraints=consts, method='SLSQP', args=(nxs, nys), tol=1e-8, options={'ftol':1e-8,'eps':1e-9})
    #print(fit['x'])
    
    
    
    nresults = fit['x']
    
    #undo the normalization of the intensities
    results = np.array(nresults)
    for ai in range(len(areas)):
        results[ai*num_gauss+(2*num_gauss):ai*num_gauss+(3*num_gauss)] = np.multiply(nresults[ai*num_gauss+(2*num_gauss):ai*num_gauss+(3*num_gauss)],areas[ai])
    #undo energy normalization
    results[0:2*num_gauss]*=xmax
    results[5*num_gauss:]*=xmax
    
    details = pd.DataFrame()
    details['Energy'] = results[0:num_gauss]
    details['Widths'] = results[num_gauss:2*num_gauss]
    for i in range(len(ys)):
        ylab = 'Inten_y'+str(i)
        details[ylab] = results[(i+2)*num_gauss:(i+3)*num_gauss]
    #Calc fwhm and oscillator strengths (f)
    fwhms = 2*np.sqrt(2*np.log(2))*details['Widths']
    fs = 4.61e-9*details['Inten_y0']*np.exp(-1/(2*(details['Widths'])**2))*fwhms
    details['fwhm'] = fwhms
    details['f'] = fs
    #add Abs max value    
    #details['y0_max'] = np.multiply(np.exp(-1/2*1/(details['Widths'])**2),details['Inten_y0'])
    
    return results, details, fit