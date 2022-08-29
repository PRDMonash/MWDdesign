from scipy import interpolate

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from tools.basic_version2 import create_hill

class experiment():
    def __init__(self, monomer_molar, syringevolume, syringe_rafteq, reactor_volume, raft_mw = 350.6, monomer_mw = 100.117):
        self.monomer_molar = monomer_molar # Concentration Monomer over two syringes (over volume of both syringes)
        self.syringe_volume = syringevolume # same for both syringes
        self.syringe_rafteq = syringe_rafteq # DP if both flowrates are the same
        self.reactor_volume = reactor_volume #r Reactor volume

        self.total_volume = 2 * syringevolume
        self.raft, self.monomer, self.aibn = raft_mw, monomer_mw, 164.21

    def mol_monomer(self):
        return self.monomer_molar* self.total_volume/ 1000
    
    def mol_RAFT(self):
        return self.mol_monomer() / self.syringe_rafteq
    
    #monomer
    def monomer_molar_syringe(self):
        return self.mol_monomer() / self.syringe_volume *1000
    
    def monomer_molpermL(self):
        return self.monomer_molar_syringe() / 1000
    
    def set_flowrateMonomer(self, flowrate):
        self.flowrateMonomer = flowrate

    def monomer_molpermin(self):
        return self.monomer_molpermL() * self.flowrateMonomer
    
    #RAFT
    def raft_molar_syringe(self):
        return self.mol_RAFT() / self.syringe_volume * 1000
    
    def raft_molpermL(self):
        return self.raft_molar_syringe() / 1000
    
    def set_flowrateRaft(self, flowrate):
        self.flowrateRaft = flowrate

    def raft_molpermin(self):
        return self.raft_molpermL() * self.flowrateRaft

    #DP
    def DP_syringe(self):
        return self.monomer_molpermin() / self.raft_molpermin()
    
    def total_flowrate(self):
        try:
            return self.flowrateMonomer + self.flowrateRaft
        except:
            print("set flowrates")
    
    def tres(self):
        return self.reactor_volume / self.total_flowrate() #in min

    #DP to flowrate
    def DPtoFlowrate(self, DP):
        mol_min_raft = self.monomer_molpermin()/ DP
        return mol_min_raft / self.raft_molpermL()

    def tresDP(self, DP):
        try:
            return self.reactor_volume / (self.flowrateMonomer + self.DPtoFlowrate(DP))
        except:
            print("SET flowrates")
    
    def totalflowrate_DP(self, DP):
        return self.flowrateMonomer + self.DPtoFlowrate(DP)

    def raft_molpermin_DP(self, DP):
        return self.raft_molpermL() * self.DPtoFlowrate(DP)

    def polymer_molpermL(self, DP):
        return self.raft_molpermin_DP(DP) / self.totalflowrate_DP(DP)
    
    def polymer_min(self, DP, n):
        """
        time to pump for n mol of polymer, in min
        """
        return n /self.raft_molpermin_DP(DP) # in min

csv = 'MIXC_Experimental_Calibrated_Distributions.csv'
DF = pd.read_csv(csv)

def normalize(y):
    '''
    Normalize list
    '''
    return [(i-min(y))/(max(y)-min(y)) for i in y]

def plot_MlogM(ax):
    mwd = DF['MWD']
    MlogM = DF['MlogM']

    ax.plot(mwd, MlogM)
    ax2.set_xscale('log')

def plot_number(ax):
    mwd = DF['MWD']
    numberdistr = DF['sum_distribution']

    ax.plot(mwd, numberdistr)


def moving_average(distr_dic:dict, window_size = 5):
    if not isinstance(window_size, int) or (window_size%2 == 0):
        raise ValueError('window_size needs to be an odd number.')
    
    i = 0
    x_moving_averages = []
    moving_averages = []
    print(len(distr_dic.values()))
    while i < len(distr_dic.values()) - window_size:
        window = list(distr_dic.values())[i:i+window_size]
        print(window)
        index_x = int((window_size/2)+0.5+i)
        print(index_x)

        x_moving_averages.append(list(distr_dic.keys())[index_x])
        average = sum(window)/ len(window)
        moving_averages.append(average)
        i +=1
    return x_moving_averages, moving_averages

def add_rectangular(ax, start_x, start_y, width, height, color = 'blue',alpha = 0.2):
    ax.add_patch(
            Rectangle(xy=(start_x,start_y), width=width,
                      height=height, alpha = alpha, facecolor= color, edgecolor='black'))

def plot_number_split(ax):
    mwd = DF['MWD']
    numberdistr = normalize(DF['sum_distribution'])

    # only usefull points in distribution
    for i, number in enumerate(numberdistr):
        if number < 0.001:
            print(i)
            break
    print("Usefull distribution:")

    mwd_x = list(mwd[:i])
    numberdistri_y = list(numberdistr[:i])
    lenght = len(numberdistri_y)
    resolution = AMOUNT_OF_SUBDISTR
    subdistri_width = (max(mwd_x) - min(mwd_x))/resolution
    print('MWD range: {} - {}'.format(min(mwd_x), max(mwd_x)))
    print('Total amount of datapoints: {}'.format(lenght))
    print('Amount of subdistributions: {}'.format(resolution))
    print('Delta MW: {}'.format((max(mwd_x) - min(mwd_x))/resolution))
    print('----------')

    ax.plot(mwd_x, numberdistri_y)

    splits = []
    start_distri = min(mwd_x)
    while start_distri <= max(mwd_x):
        splits.append(round(start_distri,0))
        start_distri += subdistri_width

   
    subdistributions = [split - (subdistri_width/2) for split in splits][1:]

    subdistr_conc_dic = {}
    for dotline in subdistributions:
        y_point = interpolate_data(mwd_x, numberdistri_y, dotline)
        ax.vlines(dotline, 0, y_point, linestyle = 'dashed')
        subdistr_conc_dic.update({dotline:y_point})
        add_rectangular(ax, dotline - subdistri_width/2, 0, subdistri_width,y_point)
    ax.set_title("Number of DPs: {}".format(len(subdistributions)))
    ax.set_xlabel('MW')
    ax.set_ylabel('Number of molecules (n)')
    print("Number of DPs: {}".format(len(subdistributions)))
    
    print(subdistr_conc_dic)

    #f(mwd_x, numberdistri_y, 2000)

    

    return subdistr_conc_dic

def interpolate_data(xs, ys, datapoint_x):
    """
    Linear Interpolation
    """
    dict_data = dict(zip(xs, ys))
    previous = 0
    for i, (key, value) in enumerate(dict_data.items()):
        if previous <datapoint_x <key:
            x1 = (list(dict_data.keys())[i])
            x = (list(dict_data.keys())[i-1])
            y1= (list(dict_data.values())[i])
            y =(list(dict_data.values())[i-1])
            y_new = y1 + (datapoint_x-x1) * ((y - y1)/(x -x1))

            return y_new
        previous = key
    return 0

def toDPs(dictionar:dict):
    dic = {((i-RAFT_MW)/MONOMER_MW):j for i,j in dictionar.items()}
    return dic

def plot_experiment(ax, dp_dic:dict):
    exp = experiment(3, 5, 30, 0.5)
    monomer_flowrate = 0.016
    raft_fr, monomer_fr = [], []
    tress = []
    times = []
    exp.set_flowrateMonomer(monomer_flowrate)
    current_time = 0
    for dp, conc in dp_dic.items():
        flowrate = exp.DPtoFlowrate(dp) #raft flowrate
        raft_fr.append(flowrate)
        
        monomer_fr.append(monomer_flowrate)

        tress.append(exp.tresDP(dp))

        time = (exp.polymer_min(dp, conc))
        time = time / 30000
        times.append(time)
        
        ax.plot([current_time, current_time + time], [flowrate, flowrate], color = 'red')
        current_time += time

    average_tres = sum(tress)/len(tress)
    ax.set_title('Average tres: {} min'.format(round(average_tres, 2)))
    ax.plot([0, current_time], [monomer_flowrate, monomer_flowrate], label = 'Monomer')
    ax.legend()
    ax.set_xlabel("time / min")
    ax.set_ylabel('Flowrate')

    return monomer_fr , raft_fr, times
def plot_experiment_hill(ax, dp_dic:dict):
    print(dp_dic)

    
    dps, concs = list(dp_dic.keys()), list(dp_dic.values())
    dps_hill = create_hill(dps, plot=False)
    dp_dic = {}
    for dp in dps_hill:
        for i, dp_old in enumerate(dps):
            if dp == dp_old:
                conc = concs[i]
                dp_dic.update({dp:conc})

    print(dp_dic)
    exp = experiment(3, 5, 30, 0.5)
    monomer_flowrate = 0.016
    raft_fr, monomer_fr = [], []
    tress = []
    times = []
    exp.set_flowrateMonomer(monomer_flowrate)
    current_time = 0
    for dp, conc in dp_dic.items():
        flowrate = exp.DPtoFlowrate(dp) #raft flowrate
        raft_fr.append(flowrate)
        
        monomer_fr.append(monomer_flowrate)

        tress.append(exp.tresDP(dp))

        time = (exp.polymer_min(dp, conc))
        time = time / 30000
        times.append(time)
        
        ax.plot([current_time, current_time + time], [flowrate, flowrate], color = 'red')
        current_time += time
    average_tres = sum(tress)/len(tress)
    ax.set_title('Average tres: {} min'.format(round(average_tres, 2)))
    ax.plot([0, current_time], [monomer_flowrate, monomer_flowrate], label = 'Monomer')
    ax.legend()
    ax.set_xlabel("time / min")
    ax.set_ylabel('Flowrate')

def save_csv(fr_monomer, fr_raft, times):
    times = [round((time * 60),0) for time in times] #given in minutes, convert to seconds
    d = {'raft': fr_raft, 'monomer':fr_monomer, 'times':times}
    df = pd.DataFrame(d)
    #print(df)
    save = input('Save experiment? y/n')
    if save == 'y':
        namecsv = input('Name of csv file (excluding .csv)\n>> ')
        df.to_csv('{}.csv'.format(namecsv))

def calculate_Mn(subs:dict):
    a = [dp * conc for dp, conc in subs.items()]
    b = sum(list(subs.values()))
    print(sum(a)/b)

ax1 = plt.subplot(2,2,1)
ax2 = plt.subplot(2,2,2)
ax3 = plt.subplot(2,2,3)
ax4 = plt.subplot(2,2,4)

AMOUNT_OF_SUBDISTR = 40


MONOMER_MW = 100.117
RAFT_MW = 350.6

plot_number(ax1)
plot_MlogM(ax2)
subdistri_dict = plot_number_split(ax3)
calculate_Mn(subdistri_dict)
DP_dict = toDPs(subdistri_dict)

#monomer, raft, times = plot_experiment(ax4, DP_dict)
plot_experiment_hill(ax4, DP_dict)

#save_csv(monomer, raft, times)
plt.show()

exit()

exp = experiment(3, 5, 30, 0.5)
print(exp)

exp.set_tres_monomersyringe(20)
exp.set_concentrations(list(DP_dict.values()))
exp.set_DPs(list(DP_dict.keys()))
exp.generate_plugs()
exp.save_experiment("test.png")

