import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

file = "C:/Users/jvan0036/Python/scripts/simulation/Experiments_2022/program/Flowrates_1402_MIXB_FLOWcsv.csv"

class flowrate_csv():
    def __init__(self, file, plugvolume = 0.25) -> None:
        if not os.path.exists(file):
            print(f'{file} does not exist...')

        self.plugvolume = plugvolume
        self.df = pd.read_csv(file)
        self.raft, self.monomer = self.get_flowrates()
        self.totalFRdict = self.create_totalFRDict()
        self.totaltimedict = self.create_timeFRDict()

    def get_flowrates(self, lst = True):
        raft_flowrates = self.df['RAFT']
        monomer_flowrates = self.df['Monomer']
        if lst:
            return list(raft_flowrates), list(monomer_flowrates)
        return raft_flowrates, monomer_flowrates
    
    def create_totalFRDict(self, debug = True):
        totalFRdict = {raft+monomer: [raft, monomer] for raft, monomer in zip(self.raft, self.monomer)}
        if debug:
            print('Times : flowrates')
            for time, flowrates in totalFRdict.items():
                print(f'{time} : {flowrates}')
        return totalFRdict
    
    def create_timeFRDict(self, in_sec = False, debug = True):
        print(self.totalFRdict)
        print('plugvolume : {}'.format(self.plugvolume))
        totaltimedict = {self.plugvolume/total_FR : flowrates for total_FR, flowrates in self.totalFRdict.items()}
        if debug:
            print('Times : flowrates')
            for time, flowrates in totaltimedict.items():
                print(f'{time} : {flowrates}')
        if in_sec:
            totaltimedict = {time_min  * 60 : flowrates for time_min, flowrates in totaltimedict.items()}
        return totaltimedict

    def plot_flowrates(self):
        ax = plt.subplot(1,1,1)
        
        ax.set_xlabel('Reaction Time / miutes', fontsize = 13)
        ax.set_ylabel('flowrate / mL \u00B7 min⁻¹', fontsize = 13)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        ax.xaxis.set_ticks_position('bottom')

        ax.set_xlim(-1,32)
        ax.set_xticks(np.arange(0, 35, 5))

        ax.set_ylim(0,0.02)
        ax.set_yticks(np.arange(0, 0.021, 0.005))

        starttime = 0
        for time, flowrates in self.totaltimedict.items():    
            raft, monomer = flowrates[0], flowrates[1]
            #print(f'monomer {monomer}at and raft at {raft} from {starttime} for {time}')
            if raft == monomer:
                ax.plot([starttime, starttime + time], [raft, raft], c='goldenrod', lw = 6)
            ax.plot([starttime, starttime + time], [raft, raft], c='goldenrod', lw = 2)
            ax.plot([starttime, starttime + time], [monomer, monomer], c='navy', lw = 2)
            starttime += time
        plt.show()

a =flowrate_csv(file, plugvolume=0.25)

a.get_flowrates()


a.plot_flowrates()