import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
import pickle
from pint import UnitRegistry
import os
from dotmap import DotMap
import sys  


def make_24Mpouch(**kwargs):

    import numpy as np
    
    def get_nominal(value):
        try:
            nominal = value.n
        except AttributeError:
            nominal = value.magnitude

        return nominal
    
    
    cell = {
        "name": 'name',
        "format": 'single layer pouch',
        "cellstack": 'missing',
        "ecapratio": 'missing', #electrolyte:capacity ratio (~ 1.8 mL/Ah)
        
        "height": 'missing',
        "width": 'missing',
        "nlayers": 'missing', #depth is set by number of layers as pouch cell is unconstrained (1 layer)
        
        "pouchthick": 'missing',
        "pouchdens": 'missing',
        "pouchclearance": 'missing', #clearance on edges for sealing pouch cell (~ 1 mm)
        
        "tabh": 'missing', #tab dimensions
        "tabw": 'missing',
        "tabt": 'missing',
        "tabloc": 'top',
        
        "tabdenspos": 'missing', #2.7 g/cm3 for aluminum
        "tabdensneg": 'missing', #8.9 g/cm3 for nickel/copper
    
        "llifactor": 0.95, #formation losses
        "extramass": 'missing', #unaccounted for mass like tabs etc
    }
    
    for key, value in kwargs.items(): #Load specified properties from arguments
        cell[key] = value
    cell = DotMap(cell)
    unit = cell.unit
    if any(x == 'missing' for x in cell.values()): #Check if any are unspecified
     raise ValueError('Unspecified cell properties.')

    area = cell.nlayers*cell.width*cell.height #jelly roll area
    
    #Capacity
    capacity = min(cell.cellstack.positive.composite.arealcap,cell.cellstack.negative.composite.arealcap) * cell.llifactor * 1*area #double coated
    capacity.ito(unit.A*unit.hr)
    
    #Energy
    energy = capacity*(cell.cellstack.positive.composite.active.avgE-cell.cellstack.negative.composite.active.avgE)
    energy.ito(unit.W*unit.hr)
    
    #Assign component and cell masses
    elytemass = cell.ecapratio*capacity*cell.cellstack.electrolyte.density
    elytemass.ito(unit.g)
    posmass = area*(1*cell.cellstack.positive.composite.thick*cell.cellstack.positive.composite.density) #double coated
    posmass.ito(unit.g)
    posccmass = area*(cell.cellstack.positive.currentcollector.thick*cell.cellstack.positive.currentcollector.density)
    posccmass.ito(unit.g)
    negmass = area*(1*cell.cellstack.negative.composite.thick*cell.cellstack.negative.composite.density) #double coated
    negmass.ito(unit.g)
    negccmass = area*(cell.cellstack.negative.currentcollector.thick*cell.cellstack.negative.currentcollector.density)
    negccmass.ito(unit.g)
    sepmass = area*(1*cell.cellstack.separator.thick*cell.cellstack.separator.density) #2 sided
    sepmass.ito(unit.g)
    jellymass = posmass + posccmass + negmass + negccmass + sepmass + elytemass
    jellymass.ito(unit.g)
    
    pouchmass = cell.pouchdens*cell.pouchthick*(cell.width+cell.pouchclearance)*(cell.height+cell.pouchclearance)*2 #2 layer pouch sealed together
    tabmass = (cell.tabh+cell.pouchclearance)*cell.tabw*cell.tabt*cell.tabdenspos + (cell.tabh+cell.pouchclearance)*cell.tabw*cell.tabt*cell.tabdensneg
    casemass = pouchmass + tabmass + cell.extramass
    casemass.ito(unit.g)
    
    cellmass = casemass + jellymass
    
    #Stackthick
    stackthick = 1*cell.cellstack.positive.composite.thick + cell.cellstack.positive.currentcollector.thick  \
                    + 1*cell.cellstack.negative.composite.thick + cell.cellstack.negative.currentcollector.thick  \
                        + 1*cell.cellstack.separator.thick
    stackthick.ito(unit.cm)
    
    #Cell volume
    depth = cell.nlayers*stackthick + 2*cell.pouchthick
    pouchvol = (cell.width + cell.pouchclearance) * (cell.height + cell.pouchclearance) *  depth
    
    tabvol = (cell.tabh+cell.pouchclearance)*cell.tabw*cell.tabt*2
    volume = pouchvol + tabvol
    volume.ito(unit.cm**3)
    
    #Stack thickness
    stackthick.ito(unit.um)
    
    #NP Ratio
    cell.cellstack.positive.composite.arealcap.ito(unit.mA*unit.hr/unit.cm**2)
    cell.cellstack.negative.composite.arealcap.ito(unit.mA*unit.hr/unit.cm**2)
    NPratio = get_nominal(cell.cellstack.negative.composite.arealcap)/get_nominal(cell.cellstack.positive.composite.arealcap)
    
    #Assign dict
    cell.jlarea = area
    cell.capacity = capacity
    cell.energy = energy
    cell.volume = volume
    cell.stackthick = stackthick
    cell.NPratio = NPratio
    cell.depth = depth
    
    cell.mass.total = cellmass
    cell.mass.jellyroll = jellymass
    cell.mass.case = casemass
    cell.mass.electrolyte = elytemass
    cell.mass.positive = posmass
    cell.mass.positivecc = posccmass
    cell.mass.negative = negmass
    cell.mass.negativecc = negccmass
    cell.mass.separator = sepmass
    
    return cell



#PLOT ELECTRODE STACK THICKNESSES 
def plot_singlelayer_thickbreakdown(cell):
    import numpy as np 
    import matplotlib.pyplot as plt
    figsz = (11,1)
    #     figsz = (20,70)
    dpi = 100

    def get_nominal(value):
        try:
            nominal = value.n
        except AttributeError:
            nominal = value.magnitude
        return nominal

    unit = cell.unit
    cellstack = cell.cellstack

    cellstack.separator.thick.ito(unit.um)
    sep = get_nominal(cellstack.separator.thick)

    cellstack.positive.composite.thick.ito(unit.um)
    pos = get_nominal(cellstack.positive.composite.thick)

    cellstack.positive.currentcollector.thick.ito(unit.um)
    poscc = get_nominal(cellstack.positive.currentcollector.thick)

    cellstack.negative.composite.thick.ito(unit.um)
    neg = get_nominal(cellstack.negative.composite.thick)

    cellstack.negative.currentcollector.thick.ito(unit.um)
    negcc = get_nominal(cellstack.negative.currentcollector.thick)

    thicks = np.array([poscc, pos, sep, neg, negcc])

    cumthick = thicks.cumsum()

    category_colors = plt.get_cmap('RdYlGn')(
            np.linspace(0.15, 0.85, np.size(thicks)))

    labels = ['']
    height = 1

    category_names = ['positive c.c.',
                      'positive electrode','separator','negative electrode',
                      'negative c.c.'] 


    fig, ax = plt.subplots(figsize=figsz, dpi=dpi) 



    i = 0
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='#7B7B7B')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 1
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='#001a33',
            hatch='..',
            edgecolor='k')
    ax.barh(labels,thicc, left=starts, height=height, 
            color='none',edgecolor='k')

    i = 2
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            color='w',
            hatch='xxxx')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 3
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color='k',
            hatch='--',
            edgecolor='#7B7B7B')
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')

    i = 4
    thicc = thicks[i]
    starts = cumthick[i] - thicc
    if cell.cellstack.negative.currentcollector.name=='Cu':
        cc_color = '#D08840'
    else: 
        cc_color = '#7B7B7B'

    ax.barh(labels,thicc, left=starts, height=height,
            label=category_names[i], 
            color=cc_color)
    ax.barh(labels,thicc, left=starts, height=height,
            color='none',edgecolor='k')


    ax.set_ylim((-0.5,0.5))
    ax.set_xlim((0,cumthick[-1]))
    ax.set_xlabel('μm')
 
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.8),
              fancybox=True, ncol=5, frameon=False)  