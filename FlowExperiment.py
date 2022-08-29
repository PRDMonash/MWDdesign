
class experiment():
    def __init__(self, monomer_molar, syringevolume, syringe_rafteq, reactor_volume, raft_mw = 350.6, monomer_mw = 100.117, monomer_density = 0.94):
        self.monomer_molar = monomer_molar # Concentration Monomer over two syringes (over volume of both syringes)
        self.syringe_volume = syringevolume # same for both syringes
        self.syringe_rafteq = syringe_rafteq # DP if both flowrates are the same
        self.reactor_volume = reactor_volume #r Reactor volume

        self.total_volume = 2 * syringevolume
        self.raft, self.monomer, self.aibn = raft_mw, monomer_mw, 164.21
        self.monomer_density = monomer_density
        self.summary_string = ''

    def mol_monomer(self):
        return self.monomer_molar* self.total_volume/ 1000

    def gram_monomer(self):
        return self.mol_monomer() * self.monomer
    
    def mL_monomer(self, rounding = 2):
        return round(self.gram_monomer() / self.monomer_density, rounding)
    
    def mol_RAFT(self):
        return self.mol_monomer() / self.syringe_rafteq
    
    def gram_raft(self):
        return self.mol_RAFT() * self.raft
    
    #monomer
    def monomer_molar_syringe(self):
        return self.mol_monomer() / self.syringe_volume *1000
    
    def monomer_molpermL(self):
        return self.monomer_molar_syringe() / 1000
    
    def set_flowrateMonomer(self, flowrate):
        self.flowrateMonomer = flowrate

    def monomer_molpermin(self, rounding = 6):
        return round(self.monomer_molpermL() * self.flowrateMonomer, rounding)
    
    #RAFT
    def raft_molar_syringe(self):
        return self.mol_RAFT() / self.syringe_volume * 1000
    
    def raft_molpermL(self):
        return self.raft_molar_syringe() / 1000
    
    def set_flowrateRaft(self, flowrate):
        self.flowrateRaft = flowrate

    def raft_molpermin(self):
        return self.raft_molpermL() * self.flowrateRaft

    def setDP(self, DP):
        self.DP = DP
    #DP
    def DP_syringe(self):
        return self.monomer_molpermin() / self.raft_molpermin()
    
    def total_flowrate(self):
        try:
            return self.flowrateMonomer + self.raft_flowrate()
        except:
            print("set flowrates")
    
    def tres(self, rounding = 2):
        return round(self.reactor_volume / self.total_flowrate(),2) #in min

    #DP to flowrate
    def raft_molpermin_DP(self):
        return self.monomer_molpermin() / self.DP

    def raft_flowrate(self):
        return self.raft_molpermin_DP() / self.raft_molpermL()

    def tresDP(self, DP):
        try:
            return self.reactor_volume / (self.flowrateMonomer + self.DPtoFlowrate(DP))
        except:
            print("SET flowrates")
    
    def totalflowrate_DP(self, DP):
        return self.flowrateMonomer + self.DPtoFlowrate(DP)

    #def raft_molpermin_DP(self, DP):
        return self.raft_molpermL() * self.DPtoFlowrate(DP)

    def polymer_molpermL(self):
        return self.raft_molpermin_DP() / self.total_flowrate()
     
    def setCollectionTime(self, time):
        self.collectionTime = time
    
    def polymer_mol_collected(self):
        """ Polymer collected for one reactor volume"""
        return self.polymer_molpermL() * self.reactor_volume

    def polymer_min(self, DP, n):
        """
        time to pump for n mol of polymer, in min
        """
        return n /self.raft_molpermin_DP(DP) # in min

    def get_dpDictionary(self, dps):
        self.dpsdict = {}
        for dp in dps:
            self.setDP(dp)
            collected = (self.polymer_mol_collected())
            #self.experiment_summary()
            self.dpsdict.update({dp:collected})
        return self.dpsdict

    def simulatedMWD_1reactorVolume(self, dps):
        dps_dict = self.get_dpDictionary(dps)
        dps_dict = {(dp * self.monomer) + self.raft: moles*100000 for dp,moles in dps_dict.items()} # Mn:moles
        print(dps_dict)
    
    def append_string(self,string):
        self.summary_string = self.summary_string + string +'\n'

    def solution_summary(self):
        self.append_string('Monomer Concentration: {} M'.format(self.monomer_molar))
        self.append_string('Syringe Volume: {} mL (total: {} mL)'.format(self.syringe_volume, self.total_volume))
        self.append_string('DP initial: {}'.format(self.syringe_rafteq))
        self.append_string('\tmoles RAFT: {} ({} gram)'.format(self.mol_RAFT(), self.gram_raft()))
        self.append_string('\tmoles monomer: {} ({} gram)'.format(self.mol_monomer(), self.gram_monomer()))
    
    def RAFTsyringe_summary(self):
        self.append_string('\n# RAFT SYRINGE #')
        self.append_string('RAFT: {} gram = {} moles'.format(self.gram_raft(), self.mol_RAFT()))
        self.append_string('Solvent: {} mL'.format(self.syringe_volume))
        self.append_string('--> RAFT: {} M'.format(self.raft_molar_syringe()))

    def Monomersyringe_summary(self):
        self.append_string('\n# MONOMER SYRINGE #')
        self.append_string('Momomer: {} gram = {} mL = {} moles'.format(self.gram_monomer(), self.mL_monomer(), self.mol_monomer()))
        self.append_string('Solvent: {} mL'.format(round(self.syringe_volume- self.mL_monomer(),2)))
        self.append_string('--> Monomer: {} M ({} mol/mL)'.format(self.monomer_molar_syringe(), self.monomer_molpermL()))

    def settingMonomerFlowrate_summary(self):
        self.append_string('\n# Pick constant monomer flowrate #')
        self.append_string('Flowrate Monomer: {} mL/min'.format(self.flowrateMonomer))
        self.append_string('moles flow Monomer: {} mol/min'.format(self.monomer_molpermin()))
        self.append_string('\t\t\t mL/min * mol/mL = mol/min')
    
    def selectDP_summary(self):
        self.append_string('\n# Select DP #')
        self.append_string('DP desired: {}'.format(self.DP))
        self.append_string('moles flow RAFT: {} mol/min'.format(self.raft_molpermin_DP()))
        self.append_string('\t\t\tmol/min monomer /  DP  = moles/min RAFT')
        self.append_string('Flowrate RAFT: {} mL/min'.format(self.raft_flowrate()))
        self.append_string('\t\t\tmol/min / mol/mL = mL/min')
    
    def totalFlowrate_summary(self):
        self.append_string('\n# Total Flowrate #')
        self.append_string('Total Flowrate: {} mL/min'.format(self.total_flowrate()))
        self.append_string('Residence time: {} min'.format(self.tres()))

    def mol_polymer_summary(self):
        self.append_string('\n# moles polymer collected #')
        self.append_string('moles flow polymer: {} mol/mL'.format(self.polymer_molpermL()))
        self.append_string('\t\t\tmol/min (raft)  /  mL/min (total flowrate)  = mol/mL')
        self.append_string('Collect 1 Reactor Volume ({} mL, {} min): {} mol polymer'.format(self.reactor_volume, self.tres(), self.polymer_mol_collected()))
        

    def experiment_summary(self):
        self.append_string("\t\tFLOW EXPERIMENT")
        self.append_string('-----------------')
        self.append_string("MW RAFT: {} g/mol".format(self.raft))
        self.append_string("MW monomer: {} g/mol".format(self.monomer))
        self.append_string('Reactor Volume: {} mL'.format(self.reactor_volume))
        self.append_string('-----------------')

        self.solution_summary()
        self.RAFTsyringe_summary()
        self.Monomersyringe_summary()
        self.settingMonomerFlowrate_summary()
        self.selectDP_summary()
        self.totalFlowrate_summary()
        self.mol_polymer_summary()
        print(self.summary_string)
        pass





