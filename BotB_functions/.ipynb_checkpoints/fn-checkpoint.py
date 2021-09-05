import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
import pickle
from pint import UnitRegistry
import os
from dotmap import DotMap
import sys  


# PRINT BATTERY PROPERTY DICTIONARY STRUCTURE
def print_struct(data):
    import json
    print(json.dumps(data, indent=5, default=str))
    
# ACTIVE MATERIAL STRUCTURE
def make_active(**kwargs):
    active = {
        "name": 'missing', #string
        "speccap": 'missing', #mAh/g
        "avgE": 'missing', #V
        "density": 'missing' #g/cm3 pint units
    }
    for key, value in kwargs.items(): #Load specified properties from arguments
        active[key] = value
    active = DotMap(active)
    unit = active.unit
    if any(x == 'missing' for x in active.values()): #Check if any are unspecified
     raise ValueError('Unspecified active material properties.')
    del active.unit
    return active



# COMPOSITE ELECTRODE STRUCTURE
def make_composite(**kwargs):
    composite = {
        "active": 'missing', #active dict/dotmap
        "arealcap": 'missing', #mAh/cm2
        "thick": 'missing', #cm
        "arealload": 'missing', #g/cm2
        "activefrac": 0.95, #default 95% active frac - build in binder domain next (avg PVDF + SBR is 1.45 g/cc)
        "porosity": 'missing',
        "density": 'missing' #g/cm3
    }
    for key, value in kwargs.items(): #Load specified properties from arguments
        composite[key] = value  
    composite = DotMap(composite)
    unit = composite.unit
    keylist = list(kwargs)
    composite = complete_composite(composite,keylist,unit)
    
    if any(x == 'missing' for x in composite.values()): #Check if any are unspecified
     raise ValueError('Unspecified electrode composite properties.')
    del composite.unit
    return composite


def complete_composite(composite,keylist,unit): #Checks provided properties and tries to fill the gaps
    
    def average(lst):
        return sum(lst) / len(lst)

    actdens = composite.active.density
    activefrac = composite.activefrac
    speccap = composite.active.speccap
    
    porosity = composite.porosity
    arealcap = composite.arealcap
    arealload = composite.arealload
    thick = composite.thick
    compdens = composite.density
    
    actdens = actdens*activefrac #No binder densities yet so this will suffice
    
    #POROSITY SWEEP ==============================   
    porosity_array = []
    if 'porosity' not in keylist:
        try:
            porosity = 1-(compdens/actdens)
            porosity_array.append(porosity)
        except TypeError:
            pass
        try:
            porosity = 1-(arealcap/speccap/thick/actdens)
            porosity_array.append(porosity)
        except TypeError:
            pass
        try:
            porosity = 1-(arealload/thick/actdens)
            porosity_array.append(porosity)
        except TypeError:
            pass
        array = porosity_array
        witherror = (0*array[0].units).plus_minus(0) #make sure everything has uncertainty
        array = [x+witherror for x in array]
        if len(array)==0: #No value calculated
            raise ValueError('Unspecified electrode composite properties.')
        elif  abs(1-average([x / array[0] for x in array]))>0.001: #Different values calculated
            raise ValueError('Conflicting defined electrode composite properties.')
        else:
            units = array[0].units
            try:
                value = array[0].n
                error = max([x.s for x in array]) #choose the biggest uncertainty from the array
                prop = unit.Measurement(value,error,units)
            except AttributeError:
                prop = unit.Measurement(array[0],units)
            porosity = prop
    
    #AREALLOAD SWEEP ==============================
    arealload_array = []
    if 'arealload' not in keylist:
        try:
            arealload = arealcap/speccap
            arealload.ito(unit.g/unit.cm**2)
            arealload_array.append(arealload)
        except TypeError:
            pass
        try:
            arealload = compdens*thick
            arealload.ito(unit.g/unit.cm**2)
            arealload_array.append(arealload)
        except TypeError:
            pass
        try:
            arealload = actdens*(1-porosity)*thick
            arealload.ito(unit.g/unit.cm**2)
            arealload_array.append(arealload)
        except TypeError:
            pass
        array = arealload_array
        witherror = (0*array[0].units).plus_minus(0) #make sure everything has uncertainty, even if it is 0
        array = [x+witherror for x in array]
        if len(array)==0: #No value calculated
            raise ValueError('Unspecified electrode composite properties.')
        elif  abs(1-average([x / array[0] for x in array]))>0.001: #Different values calculated
            raise ValueError('Conflicting defined electrode composite properties.')
        else:
            units = array[0].units
            try:
                value = array[0].n
                error = max([x.s for x in array]) #choose the biggest uncertainty from the array
                prop = unit.Measurement(value,error,units)
            except AttributeError:
                prop = unit.Measurement(array[0],units)
            arealload = prop
            

    
    #AREALCAP SWEEP ==============================
    arealcap_array = []
    if 'arealcap' not in keylist:
        try:
            arealcap = speccap*arealload
            arealcap.ito(unit.mA*unit.hr/unit.cm**2)
            arealcap_array.append(arealcap)
        except TypeError:
            pass
        try:
            arealcap = speccap*thick*compdens
            arealcap.ito(unit.mA*unit.hr/unit.cm**2)
            arealcap_array.append(arealcap)
        except TypeError:
            pass
        try:
            arealcap = speccap*thick*actdens*(1-porosity)
            arealcap.ito(unit.mA*unit.hr/unit.cm**2)
            arealcap_array.append(arealcap)
        except TypeError:
            pass

        array = arealcap_array
        witherror = (0*array[0].units).plus_minus(0) #make sure everything has uncertainty, even if it is 0
        array = [x+witherror for x in array]
        if len(array)==0: #No value calculated
            raise ValueError('Unspecified electrode composite properties.')
        elif  abs(1-average([x / array[0] for x in array]))>0.001: #Different values calculated
            raise ValueError('Conflicting defined electrode composite properties.')
        else:
            units = array[0].units
            try:
                value = array[0].n
                error = max([x.s for x in array]) #choose the biggest uncertainty from the array
                prop = unit.Measurement(value,error,units)
            except AttributeError:
                prop = unit.Measurement(array[0],units)
            arealcap = prop
            
    #COMPDENS SWEEP ==============================
    compdens_array = []
    if 'density' not in keylist:
        try:
            compdens = (arealcap/speccap)/thick
            compdens.ito(unit.g/unit.cm**3)
            compdens_array.append(compdens)
        except TypeError:
            pass
        try:
            compdens = actdens*(1-porosity)
            compdens.ito(unit.g/unit.cm**3)
            compdens_array.append(compdens)
        except TypeError:
            pass
        try:
            compdens = arealload/thick
            compdens.ito(unit.g/unit.cm**3)
            compdens_array.append(compdens)
        except TypeError:
            pass
        array = compdens_array
        witherror = (0*array[0].units).plus_minus(0) #make sure everything has uncertainty, even if it is 0
        array = [x+witherror for x in array]
        if len(array)==0: #No value calculated
            raise ValueError('Unspecified electrode composite properties.')
        elif  abs(1-average([x / array[0] for x in array]))>0.001: #Different values calculated
            raise ValueError('Conflicting defined electrode composite properties.')
        else:
            units = array[0].units
            try:
                value = array[0].n
                error = max([x.s for x in array]) #choose the biggest uncertainty from the array
                prop = unit.Measurement(value,error,units)
            except AttributeError:
                prop = unit.Measurement(array[0],units)
            compdens = prop
            
            
            
    
    #THICK SWEEP ==============================
    thick_array = []
    if 'thick' not in keylist:
        try:
            thick = arealload/compdens
            thick.ito(unit.um)
            thick_array.append(thick)
        except TypeError:
            pass
        try:
            thick = arealcap/speccap/compdens
            thick.ito(unit.um)
            thick_array.append(thick)
        except TypeError:
            pass
        try:
            thick = arealload/actdens/(1-porosity)
            thick.ito(unit.um)
            thick_array.append(thick)
        except TypeError:
            pass

        array = thick_array
        witherror = (0*array[0].units).plus_minus(0) #make sure everything has uncertainty, even if it is 0
        array = [x+witherror for x in array]
        if len(array)==0: #No value calculated
            raise ValueError('Unspecified electrode composite properties.')
        elif  abs(1-average([x / array[0] for x in array]))>0.001: #Different values calculated
            raise ValueError('Conflicting defined electrode composite properties.')
        else:
            units = array[0].units
            try:
                value = array[0].n
                error = max([x.s for x in array]) #choose the biggest uncertainty from the array
                prop = unit.Measurement(value,error,units)
            except AttributeError:
                prop = unit.Measurement(array[0],units)
            thick = prop
    
    #WRITE RESULTS
    composite.porosity = porosity
    composite.arealcap = arealcap
    composite.arealload = arealload
    composite.thick = thick
    composite.density = compdens
    return composite