import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
import pickle
from pint import UnitRegistry
import os
from dotmap import DotMap
import sys  


# MAKE CYLINDRICAL CELL
def make_cylindrical(**kwargs):

    import numpy as np
    
    def get_nominal(value):
        try:
            nominal = value.n
        except AttributeError:
            nominal = value.magnitude

        return nominal
    
    
    cell = {
        "name": 'name',
        "format": 'cylindrical',
        "cellstack": 'missing',
        "ecapratio": 'missing', #electrolyte:capacity ratio (mL/Ah) ~1.8 mL/Ah
        
        "diameter": 'missing',
        "height": 'missing',
        
        "canthick": 'missing',
        "candens": 'missing',
        
        "mandreldiam": 'missing', #void diameter
        "headspace": 'missing',
        
        "llifactor": 0.95, #formation losses
        "extramass": 'missing', #unaccounted for mass like tabs etc
    }
    
    for key, value in kwargs.items(): #Load specified properties from arguments
        cell[key] = value
    cell = DotMap(cell)
    unit = cell.unit
    if any(x == 'missing' for x in cell.values()): #Check if any are unspecified
     raise ValueError('Unspecified cell properties.')
    
    #Calculate: jellyroll area, cell mass, cell capacity, cell energy, cell volume,
    
    stackthick = 2*cell.cellstack.positive.composite.thick + cell.cellstack.positive.currentcollector.thick  \
                    + 2*cell.cellstack.negative.composite.thick + cell.cellstack.negative.currentcollector.thick  \
                        + 2*cell.cellstack.separator.thick
    #Component geometries
    stackthick.ito(unit.cm)
    cell.height.ito(unit.cm)
    h_cell = cell.height
    cell.mandreldiam.ito(unit.cm)
    d_mandrel = cell.mandreldiam
    cell.headspace.ito(unit.cm)
    headspace = cell.headspace
    cell.canthick.ito(unit.cm)
    canthick = cell.canthick
    cell.diameter.ito(unit.cm)
    d_cell = cell.diameter-2*canthick
    
    #Jelly roll area via Archimedes spiral maths
    a = stackthick/(2*np.pi)
    theta = (d_cell/2)*(2*np.pi)/stackthick
    l_all = (a/2)*(theta*(1+theta**2)**0.5 + np.log(theta.n + (1+theta.n**2)**0.5))
    theta_mandrel = (d_mandrel/2)*(2*np.pi)/stackthick
    l_inner = (a/2)*(theta_mandrel*(1+theta_mandrel**2)**0.5 + np.log(theta_mandrel.n + (1+theta_mandrel.n**2)**0.5))
    l_winding = l_all-l_inner
    h = h_cell-headspace-2*canthick
    area = h*l_winding
    area.ito(unit.cm**2)
    
    
    #Capacity
    capacity = min(cell.cellstack.positive.composite.arealcap,cell.cellstack.negative.composite.arealcap) * cell.llifactor * 2*area #double coat
    capacity.ito(unit.A*unit.hr)
    
    #Energy
    energy = capacity*(cell.cellstack.positive.composite.active.avgE-cell.cellstack.negative.composite.active.avgE)
    energy.ito(unit.W*unit.hr)
    
    #Assign component and cell masses
    elytemass = cell.ecapratio*capacity*cell.cellstack.electrolyte.density
    elytemass.ito(unit.g)
    posmass = area*(2*cell.cellstack.positive.composite.thick*cell.cellstack.positive.composite.density) #double coated
    posmass.ito(unit.g)
    posccmass = area*(cell.cellstack.positive.currentcollector.thick*cell.cellstack.positive.currentcollector.density)
    posccmass.ito(unit.g)
    negmass = area*(2*cell.cellstack.negative.composite.thick*cell.cellstack.negative.composite.density) #double coated
    negmass.ito(unit.g)
    negccmass = area*(cell.cellstack.negative.currentcollector.thick*cell.cellstack.negative.currentcollector.density)
    negccmass.ito(unit.g)
    sepmass = area*(2*cell.cellstack.separator.thick*cell.cellstack.separator.density) #2 sided
    sepmass.ito(unit.g)
    jellymass = posmass + posccmass + negmass + negccmass + sepmass + elytemass
    jellymass.ito(unit.g)
    d_cell = cell.diameter
    canmass = cell.candens*cell.canthick*(np.pi*(d_cell*h_cell) + 2*np.pi*(d_cell/2)**2) + cell.extramass
    canmass.ito(unit.g)
    cellmass = canmass + jellymass
    
    #Cell volume
    volume = (np.pi*(d_cell/2)**2)*h_cell
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
    
    cell.mass.total = cellmass
    cell.mass.jellyroll = jellymass
    cell.mass.case = canmass
    cell.mass.electrolyte = elytemass
    cell.mass.positive = posmass
    cell.mass.positivecc = posccmass
    cell.mass.negative = negmass
    cell.mass.negativecc = negccmass
    cell.mass.separator = sepmass
    cell.avgE = cell.cellstack.positive.composite.active.avgE-cell.cellstack.negative.composite.active.avgE
    
    return cell
    
# MAKE POUCH CELL
def make_pouch(**kwargs):

    import numpy as np
    
    def get_nominal(value):
        try:
            nominal = value.n
        except AttributeError:
            nominal = value.magnitude

        return nominal
    
    
    cell = {
        "name": 'name',
        "format": 'pouch stacked',
        "cellstack": 'missing',
        "ecapratio": 'missing', #electrolyte:capacity ratio (~ 1.8 mL/Ah)
        
        "height": 'missing',
        "width": 'missing',
        "nlayers": 'missing', #depth is set by number of layers as pouch cell is unconstrained
        
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

    area = cell.nlayers*cell.width*cell.height #jelly roll area#

    #Capacity
    capacity = min(cell.cellstack.positive.composite.arealcap,cell.cellstack.negative.composite.arealcap) * cell.llifactor * 2*area #double coated
    capacity.ito(unit.A*unit.hr)
    
    #Energy
    energy = capacity*(cell.cellstack.positive.composite.active.avgE-cell.cellstack.negative.composite.active.avgE)
    energy.ito(unit.W*unit.hr)
    
    #Assign component and cell masses
    elytemass = cell.ecapratio*capacity*cell.cellstack.electrolyte.density
    elytemass.ito(unit.g)
    posmass = area*(2*cell.cellstack.positive.composite.thick*cell.cellstack.positive.composite.density) #double coated
    posmass.ito(unit.g)
    posccmass = area*(cell.cellstack.positive.currentcollector.thick*cell.cellstack.positive.currentcollector.density)
    posccmass.ito(unit.g)
    negmass = area*(2*cell.cellstack.negative.composite.thick*cell.cellstack.negative.composite.density) #double coated
    negmass.ito(unit.g)
    negccmass = area*(cell.cellstack.negative.currentcollector.thick*cell.cellstack.negative.currentcollector.density)
    negccmass.ito(unit.g)
    sepmass = area*(2*cell.cellstack.separator.thick*cell.cellstack.separator.density) #2 sided
    sepmass.ito(unit.g)
    jellymass = posmass + posccmass + negmass + negccmass + sepmass + elytemass
    jellymass.ito(unit.g)
    
    pouchmass = cell.pouchdens*cell.pouchthick*(cell.width+cell.pouchclearance)*(cell.height+cell.pouchclearance)*2 #2 layer pouch sealed together
    tabmass = (cell.tabh+cell.pouchclearance)*cell.tabw*cell.tabt*cell.tabdenspos + (cell.tabh+cell.pouchclearance)*cell.tabw*cell.tabt*cell.tabdensneg
    casemass = pouchmass + tabmass + cell.extramass
    casemass.ito(unit.g)
    
    cellmass = casemass + jellymass
    #Stackthick
    stackthick = 2*cell.cellstack.positive.composite.thick + cell.cellstack.positive.currentcollector.thick  \
                    + 2*cell.cellstack.negative.composite.thick + cell.cellstack.negative.currentcollector.thick  \
                        + 2*cell.cellstack.separator.thick
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
    cell.avgE = cell.cellstack.positive.composite.active.avgE-cell.cellstack.negative.composite.active.avgE
    
    return cell