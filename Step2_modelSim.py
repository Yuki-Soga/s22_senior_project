#!/usr/bin/python

# PYTHON packages
from neuron import h
from neuron import gui
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys
import pickle
import time
import random 

# NEURON specific
h.load_file("nrngui.hoc")
h.load_file("import3d.hoc")
h.load_file("L5PCbiophys3.hoc")
h.load_file("TTC.hoc")

scc = sys.argv[1]
ascc = sys.argv[2]
maxTrials = 100

# Sim Time details
h.dt = 0.1
h.steps_per_ms = 1.0/h.dt
h.tstop = 500

# initial membrane voltage (both dendritic and somatic) 
h.v_init = -64

# Specific dendrite to put synapses on (change here if you want to move to another dendrite)
sectionNum = 44;
h("objref L5PC")
h.L5PC = h.TTC("cell1.asc")
h("forall nseg = 1")
h("forall e_pas = -70")
h("access L5PC.apic[" + str(sectionNum) + "]")

# function that places excitatory inputs
eStimlist = []
eSynlist = []
ePreconlist = []

def placeNMDA(location, stimVec):
  eStimlist.append(stimVec)

  eSynlist.append(h.ProbAMPANMDA2_RATIO(location))
  eSynlist[-1].gmax = 0.0004
  eSynlist[-1].mgVoltageCoeff = 0.12 #def: 0.08

  ePreconlist.append(h.NetCon(eStimlist[-1], eSynlist[-1]))
  ePreconlist[-1].weight[0] = 1
  ePreconlist[-1].delay = 0

# function that places inhibitory inputs
iStimlist = []
iSynlista = []
iPreconlist = []

def placeGABA(location, stimVec):
  iStimlist.append(stimVec)

  iSynlista.append(h.ProbUDFsyn2_lark(location))
  iSynlista[-1].tau_r = 0.18
  iSynlista[-1].tau_d = 5
  iSynlista[-1].e = - 80
  iSynlista[-1].Dep = 0
  iSynlista[-1].Fac = 0
  iSynlista[-1].Use = 0.6
  iSynlista[-1].u0 = 0
  iSynlista[-1].gmax = 0.0 #0.0007
  
  iPreconlist.append(h.NetCon(iStimlist[-1], iSynlista[-1]))
  iPreconlist[-1].weight[0] = 1
  iPreconlist[-1].delay = 0

# setup to record dendritic voltage 
dendVoltageVector = h.Vector()
dendVoltageVector.record(h.L5PC.apic[sectionNum](0.5)._ref_v)

# setup to record somatic voltage 
somaVoltageVector = h.Vector()
somaVoltageVector.record(h.L5PC.soma[0](0.5)._ref_v)

# setup to record simulated time 
timeVector = h.Vector()
timeVector.record(h._ref_t)

# Organize Syanpses on the dendrite
h("access L5PC.apic[" + str(sectionNum) + "]")
h("nseg = 50")

# where each group of excitatory synapses are located on the post-synaptic IT neuron dendrite
# (they are placed at complete random)
numGrp = 3
synPerGrp = 2
group_list = np.empty((synPerGrp*numGrp),dtype='int')
group_list[0:synPerGrp-1] = 1
group_list[synPerGrp:2*synPerGrp-1] = 2
group_list[2*synPerGrp:3*synPerGrp-1] = 3
random.shuffle(group_list)
location_list = np.arange(0.2, 0.83, 0.03)
g1_indices = np.where(group_list == 1)[0]
g2_indices = np.where(group_list == 2)[0]
g3_indices = np.where(group_list == 3)[0]

# array that stores where each group of excitatory synapses are located 
g1_location = [location_list[index] for index in g1_indices]
g2_location = [location_list[index] for index in g2_indices]
g3_location = [location_list[index] for index in g3_indices]

# Object that plays the preCell inputs 
spkTms = h.Vector()
EpreStimVec1 = h.VecStim()
EpreStimVec2 = h.VecStim()
EpreStimVec3 = h.VecStim()
IpreStimVec  = h.VecStim()
EpreStimVec1.play(h.Vector([]))
EpreStimVec2.play(h.Vector([]))
EpreStimVec3.play(h.Vector([]))

# place the synapses
for location in g1_location:
  placeNMDA(location, EpreStimVec1)
for location in g2_location:
  placeNMDA(location, EpreStimVec2)
for location in g3_location:
  placeNMDA(location, EpreStimVec3)

# Add IClamp to mimic other depolarizing inputs to the cell
h("access L5PC.soma[0]")
h("nseg = 1")
h_iclamp = h.IClamp(0.5)
h_iclamp.delay = 0
h_iclamp.dur = h.tstop
h_iclamp.amp = .13 # nA

spkTms = {}

numSites = 2 # dend and soma
DEND = 0
SOMA = 1

# Numpy MD array to save the output
VoltageTraces = np.full((numSites,maxTrials,int(h.tstop/h.dt)+1),0,dtype=float)


for featID in list(range(1, 30)):
  #work your way through all synapse groups
  for synGrpID in range(numGrp):
    fn = "Feat{0}_SynGrp{1}_{2}_{3}.pickle".format(featID, synGrpID+1, scc, ascc)
    print(fn)
    fid = open(fn,'rb')
    spt = pickle.load(fid)
    spkTms[synGrpID] = spt
    trials = len(spt)

  # then trials
  for trial in range(maxTrials): 

    # set up preSyn spikes for this trial
    EpreStimVec1.play(h.Vector(spkTms[0][trial]))
    EpreStimVec2.play(h.Vector(spkTms[1][trial]))
    EpreStimVec3.play(h.Vector(spkTms[2][trial]))

    # run simulation
    h.finitialize()
    h.run()
    vec = h.Vector() #why do we need this??
    vec.copy(dendVoltageVector)

    # pause for debugging
    wait = input("Press Enter to continue.") 

    VoltageTraces[DEND,trial,:] = np.array(dendVoltageVector)
    VoltageTraces[SOMA,trial,:] = np.array(somaVoltageVector)

# Save output to file
outFile = "PostSyn_{0}_{1}_{2}.pickle".format(scc, ascc, featID)
pickle.dump(VoltageTraces, open(outFile,'wb'))
print("Output written to {0}".format(outFile))