'''
This script is intended to be an interactive data proccessing utility for data
collected by the MPMS 3. Only magnetization and susceptibility data may be analyzed
with this script. The automation assumes that standard data collection parameters
are used.

Usage:
It is highly reccomended that both the MPMS file you wish to process and this 
script are in the same directory, otherwise you must use the absolute file
path which needs to be input as a raw string, i.e. r'lorem ipsum'.

data_object = MPMS_Dataset('MPMS_export_file.csv')
data_object.Phi_Export()
data_object.Plot(data='XT')

The class object retains all attributes calculated in processing the data; XT,
sus_T, mag_dict, etc. This fact can be taken advantage of if it is desirous to
do downstream analysis in python. However, it is highly recommended to use Phi.

This script contains a block which auto-detects and enables simplified usage 
from the command line:

    In the console (powershell, command line, etc.) navigate to the work
    directory containing a copy of this script and the desired MPMS data 
    file. (pasting can be accomplished with a right click) e.g.
        >cd ".\work_dir"
        
    Then invoke python, passing this script then the data file as positional 
    arguments:
        >python MPMS_Wizard.py "2521413413_data_file.csv"
    
    Follow the directions that appear on screen and export it for analysis in 
    Phi.
        
Completed 3/11/26 by Joshua Mengell (mengell@unm.edu), questions are less likely
to be able to be answered as time passes.    
'''
import numpy as np
import sys
import os
from matplotlib import pyplot as plt

class MPMS_Dataset:
    
    def __init__(self, path:str, mass=None, MW=None, Xdia=None):
        self.path=path
        self.data_dict={}
        self.dir=os.path.dirname(os.path.realpath(__file__))
        self.filename=os.path.basename(self.path)
        with open(self.path,'r') as file:
            top=file.readline()
            if ',' in top:
                delim=','
            elif '\t' in top:
                delim='\t'
            else:
                raise ValueError('The delimination of the .csv file could not be determined as either tab or commas.')
                
            titles=top.replace('\n', '').split(sep=delim)
            data=file.readlines()
            data=[ x.replace('\n', '').split(sep=delim) for x in data ]
            for title in titles:
                self.data_dict[title]=[]
            for row_num in range(len(data)):
                for col_num in range(len(titles)):
                    self.data_dict[titles[col_num]].append(float(data[row_num][col_num]))
       
        #Extracted data aliased
        self.M_emu=self.data_dict['DC Moment Fixed Ctr (emu)']
        self.err_M_emu=self.data_dict['DC Moment Err Fixed Ctr (emu)']
        self.H=self.data_dict['Magnetic Field (Oe)']
        self.T=self.data_dict['Temperature (K)']
        self.mass=mass
        self.MW=MW
        self.Xdia=Xdia
        self.Xpara=None
        self.XT=None
        self.sus_T=None #Sliced T range for susceptibility experiment
        self.inverse_Xpara=None
        self.ueff=None
        self.mag_dict={}
        
        
        ####################################################
        ################# Data Entry Block #################
        ####################################################
        if self.mass==None:
            s1="Enter the mass of the sample (mg) into the console below,"
            s2= "if it is not known enter None, (This will disable meaninful analysis)."
            s3='---------------------------------------------'
            print(s1,s2,s3,sep='\n')
            
            
            while True:
                
                var=input()
                
                try:
                    
                    var=float(var)
                    
                    self.mass=var
                    print('Mass Value (mg) accepted',s3,sep='\n')
                    break
                except:
                    if type(var)==str:
                        if var.upper()=='NONE':
                            break
                        elif var.upper()=='EXIT':
                            #exit()
                            break
                        else:
                            print("Unaccepted input, the only options are a numeral mass in mg or None",s3,sep='\n')
                
        if self.MW==None:
            s1="Now, enter the Molecular Weight of the sample (g/mol) into the console below, "
            s2="if it is not known enter None, (This will disable meaninful analysis)."
            s3='---------------------------------------------'
            print(s1,s3,sep='\n')
            
            
            while True:
                
                var=input()
                
                try:
                    
                    var=float(var)
                    
                    self.MW=var
                    print('Molecular Weight (g/mol) accepted',s3,sep='\n')
                    break
                except:
                    if type(var)==str:
                        if var.upper()=='NONE':
                            break
                        elif var.upper()=='EXIT':
                            #exit()
                            break
                        else:
                            print("Unaccepted input, the only options are a numeral mass in g/mol or None",s3,sep='\n')
         
        if self.Xdia==None:
            s1="Finally, enter the diamagnetic suseptibility of the sample (Xe-6 emu/mol) into the console below, "
            s2="if it is not known enter None, it will be assumed to be equal to -MW/2e-6 emu/mol (if provided)."
            s3='---------------------------------------------'
            print(s1,s2,s3,sep='\n')
            
            while True:
                
                var=input()
                
                try:
                    
                    var=float(var)
                    
                    self.Xdia=var
                    print('Diamagnetic suseptibility of the sample (Xe-6 emu/mol) accepted',s3,sep='\n')
                    break
                except:
                    if type(var)==str:
                        if var.upper()=='NONE':
                            if self.MW!=None:
                                self.Xdia=self.MW/-2
                            break
                        elif var.upper()=='EXIT':
                            #exit()
                            break
                        else:
                            print("Unaccepted input, the only options are a numeral mass in g/mol or None",s3,sep='\n')
        
        ####################################################
        ############### Data Auto Sort Block ###############
        ####################################################
        # This is only going to work properly if people use 
        # my standard sequence files.
        data_types=[]
        start_idxs=[]
        stop_idxs=[]
        #Diffrentiate between mag and sus: 
        if len(self.M_emu)!=len(self.T) or len(self.T)!=len(self.H):
            raise ValueError('Problem with MPMS file, incongruency detected in column lengths')
        
        for idx in range(len(self.M_emu)-1):
            delta_H=round(self.H[idx+1]-self.H[idx])
            delta_T=round(self.T[idx+1]-self.T[idx],ndigits=2)
            
            if abs(delta_H)<=5000 and abs(delta_H)>=1:
                #Magnetization detected
                
                #check for new dataset
                if len(start_idxs)==len(stop_idxs):
                    data_types.append(['mag',round(self.T[idx],ndigits=1)])
                    start_idxs.append(idx)
                    print(f'Magnetization Experiment ({self.T[idx]} K) detected.')
            
            if abs(delta_H)>=5000 and len(start_idxs)!=len(stop_idxs):
                #Magnetization reset
                stop_idxs.append(idx)
            
            elif delta_T<0 and len(start_idxs)!=len(stop_idxs):
                #Magnetization reset
                stop_idxs.append(idx)
            
            
            if delta_H==0 and delta_T!=0 and delta_T>0:
                #Susceptibility detected
                
                #check for new dataset
                if len(start_idxs)==len(stop_idxs):
                    print(f'Susceptibility Experiment ({self.H[idx]} Oe) detected.')
                    data_types.append('sus')
                    start_idxs.append(idx)
        #Adding last stop idx
        if len(start_idxs)!=len(stop_idxs):
            stop_idxs.append(len(self.H)-1)
        
        #Reindexing for slicing:
        stop_idxs=[n+1 for n in stop_idxs]
        
        #Time for packaging; we make the bold assumption there is only 1 
        #susceptibility experiment.
        if self.mass!=None and self.MW!=None:
            for idx in range(len(data_types)):
                data_class=data_types[idx]
                if type(data_class)==list:
                    #Magnetization
                    exp, T = data_class
                    H = self.H[start_idxs[idx]:stop_idxs[idx]]
                    M = self.M_emu[start_idxs[idx]:stop_idxs[idx]]
                    NB = [m*self.MW/((self.mass/1000)*5585) for m in M]
                    self.mag_dict[T]=[H,NB]
                elif data_class=='sus':
                    H = self.H[start_idxs[idx]:stop_idxs[idx]]
                    M = self.M_emu[start_idxs[idx]:stop_idxs[idx]]
                    T = self.T[start_idxs[idx]:stop_idxs[idx]]
                    #Converting to ndarray for easy calculation revert later
                    X = [m*self.MW/(self.mass/1000) for m in M]/np.multiply(1,H)
                    X = X - self.Xdia*0.000001
                    self.sus_T=T
                    self.Xpara=X.tolist()
                    self.XT=(X*T).tolist()
                    self.inverse_Xpara=(1/X).tolist()
                    self.ueff=(2.828*np.sqrt(X*T)).tolist()
                else:
                    print("I have now clue how you've gotten here")
        else:
            print("Neither mass or molecular weight have been supplied, units of data remain in emu - no quantitative analysis is possible")
            for idx in range(len(data_types)):
                data_class=data_types[idx]
                if type(data_class)==list:
                    #Magnetization
                    exp, T = data_class
                    H = self.H[start_idxs[idx]:stop_idxs[idx]]
                    M = self.M_emu[start_idxs[idx]:stop_idxs[idx]]
                    NB = M
                    self.mag_dict[T]=[H,NB]
                elif data_class=='sus':
                    H = self.H[start_idxs[idx]:stop_idxs[idx]]
                    M = self.M_emu[start_idxs[idx]:stop_idxs[idx]]
                    T = self.T[start_idxs[idx]:stop_idxs[idx]]
                    #Converting to ndarray for easy calculation revert later
                    X=np.multiply(M,1)/H
                    self.sus_T=T
                    self.Xpara=X.tolist()
                    self.XT=(X*T).tolist()
                    self.inverse_Xpara=(1/X).tolist()
                    self.ueff=(2.828*np.sqrt(X*T)).tolist()
                else:
                    print("I have now clue how you've gotten here")
        print('Data Comprehension Complete.')
        
    def Export(self):
        '''
        Export processed data for generic purposes.
        '''
        
        
        common_name= self.filename.split('.')[0]
        #Susceptibility
        with open(os.path.join(self.dir, common_name+'_sus.txt'),'w') as file:
            file.write('Temperature (K)\tXT(emu·mol-1·K)\tX(emu·mol-1·K)\n')
            for idx in range(len(self.sus_T)):
                file.write(f"{self.sus_T[idx]}\t{self.XT[idx]}\t{self.Xpara[idx]}\n")
        
        with open(os.path.join(self.dir, common_name+'_sus.txt'),'w') as file:
            title_str=''
            title_str+='Field (T)\t'
            for temp in self.mag_dict.keys():
                title_str+=f"{temp}K Magnetization (NB)\t"
            title_str+='\n'
            file.write(title_str)
            
            for idx in range(len(self.mag_dict[list(self.mag_dict.keys())[0]][0])):
                temp_str=''
                #Field in T
                temp_str+=f"{self.mag_dict[list(self.mag_dict.keys())[0]][0][idx]/10000}\t"
                for key in self.mag_dict.keys():
                    temp_str+=f"{self.mag_dict[key][1][idx]}\t"
                temp_str+='\n'
                file.write(temp_str)
    
    def Phi_Export(self):
        '''
        Export processed data for analysis in Phi.
        '''
        
        common_name= self.filename.split('.')[0]
        phi_dir=os.path.join(self.dir, f'Phi_dir_{common_name}')
        if not os.path.isdir(phi_dir):
            os.makedirs(phi_dir)
        #Susceptibility
        with open(os.path.join(phi_dir, common_name+'_sus.exp'),'w') as file:
            for idx in range(len(self.sus_T)):
                file.write(f"{self.sus_T[idx]}\t{self.XT[idx]}\n")
        #Magnetization
        with open(os.path.join(phi_dir, common_name+'_mag.exp'),'w') as file:
            for idx in range(len(self.mag_dict[list(self.mag_dict.keys())[0]][0])):
                temp_str=''
                #Field in T
                temp_str+=f"{self.mag_dict[list(self.mag_dict.keys())[0]][0][idx]/10000}\t"
                for key in self.mag_dict.keys():
                    temp_str+=f"{self.mag_dict[key][1][idx]}\t"
                temp_str+='\n'
                file.write(temp_str)
        
        with open(os.path.join(phi_dir, common_name+'.input'),'w') as file:
            pass
                
        
        
    def Plot(self,data=None):
        plot_types=['XT','Mag','reduced_Mag','X','1/X','ueff']
        plot_types_upper=[x.upper() for x in plot_types]
        
        while data==None:
            s1='No dataset selected, select one of the following by typing it into the console below:'
            s2='---------------------------------------------'
            
            
            print(s1,plot_types,s2,sep='\n')
            temp=input()
            if temp.upper() not in plot_types_upper:
                print(f"{temp} is not an allowed input; {plot_types}")
            if temp.upper() in plot_types_upper:
                data=temp
                break
        
        fig,ax = plt.subplots()
        
        if data.upper()=='XT'.upper():
            ax.scatter(self.sus_T,self.XT)
            ax.set_ylabel('XT (emu·mol$^{-1}$·K)')
            ax.set_xlabel('Temperature (K)')
            ax.set_title(f"{self.filename.split('.')[0]} Susceptibility")
        
        if data.upper()=='Mag'.upper():
            for key in self.mag_dict.keys():
                H,NB = self.mag_dict[key]
                ax.scatter(H,NB,label=f"{key} K")
            ax.set_ylabel('Magnetization (NB·mol$^{-1}$)')
            ax.set_xlabel('Field (Oe)')
            ax.set_title(f"{self.filename.split('.')[0]} Magnetization")
            ax.legend()
        
        if data.upper()=='reduced_Mag'.upper():
            for key in self.mag_dict.keys():
                H,NB = self.mag_dict[key]
                reduced_H = np.multiply(H, 1/float(key)).tolist()
                
                ax.scatter(reduced_H,NB,label=f"{key} K",color='k')
            ax.set_ylabel('Magnetization (NB·mol$^{-1}$)')
            ax.set_xlabel('Reduced Field (Oe·K$^{-1}$)')
            ax.set_title(f"{self.filename.split('.')[0]} Magnetization")
            #ax.legend()
        
        if data.upper()=='X'.upper():
            ax.scatter(self.sus_T,self.Xpara)
            ax.set_ylabel('X (emu·mol$^{-1}$)')
            ax.set_xlabel('Temperature (K)')
            ax.set_title(f"{self.filename.split('.')[0]} Susceptibility")
        
        if data.upper()=='1/X'.upper():
            ax.scatter(self.sus_T,self.inverse_Xpara)
            ax.set_ylabel('1/X (emu$^{-1}$·mol)')
            ax.set_xlabel('Temperature (K)')
            ax.set_title(f"{self.filename.split('.')[0]} Susceptibility")
        
        if data.upper()=='ueff'.upper():
            ax.scatter(self.sus_T,self.ueff)
            ax.set_ylabel('μ$_{eff}$')
            ax.set_xlabel('Temperature (K)')
            ax.set_title(f"{self.filename.split('.')[0]} Susceptibility")
        plt.show()
            
            
###############################################################################
######################### Command line runtime Code ###########################
###############################################################################
############ If something is broken try commenting this block out #############
###############################################################################
if len(sys.argv)>1:
    print('Command line operation detected')
    file_name = sys.argv[1]
    mpms_obj=MPMS_Dataset(file_name)
    while True:
        bar='--------------------------------------------------------'
        lines=['To export for analysis in Phi, type Phi below.',
               'To export for generic plotting purpouses, type Kgraph.',
               'To plot processed data, type Plot.',
               'To exit, type Exit.']
        print(bar,*lines,bar,sep='\n')
        good_inputs=['Phi'.upper(), 
                     'Kgraph'.upper(), 
                     'Plot'.upper(), 
                     'Exit'.upper()]
        
        runtime_var = input()
        if runtime_var.upper() not in good_inputs:
            print("Bad input.")
        
        else:
            if runtime_var.upper()=='Phi'.upper():
                mpms_obj.Phi_Export()
            if runtime_var.upper()=='Kgraph'.upper():
                mpms_obj.Export()
            if runtime_var.upper()=='Plot'.upper():
                mpms_obj.Plot()
            if runtime_var.upper()=='Exit'.upper():
                break       
###############################################################################




