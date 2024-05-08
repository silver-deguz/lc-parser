from dataclasses import dataclass
from typing import Optional, Dict
import numpy as np 

InjectionKeys = {
    "Data Vault",
    "Injection",
    "Injection Number",
    "Position",
    "Comment",
    "Processing Method",
    "Instrument Method",
    "Type",
    "Status",
    "Injection Date",
    "Injection Time",
    "Injection Volume (µL)",
    "Dilution Factor",
    "Weight",
}

ChromatogramDataKeys = {
    "Time Min. (min)",
    "Time Max. (min)",
    "Data Points",
    "Detector",
    "Generating Data System",
    "Exporting Data System",
    "Operator",
    "Signal Quantity",
    "Signal Unit",
    "Signal Min.",
    "Signal Max.",
    "Channel",
    "Driver Name",
    "Channel Type",
    "Min. Step (s)",
    "Max. Step (s)",
    "Average Step (s)",
}
    
SignalParameterKeys = {
    "Signal Info",
}

"""
Base parser. Parses lines into their metadata represented as dict

Parameters:
lines (list of strings): block lines from a text file
datatype_keys (set): one of InjectionKeys, ChromatogramDataKeys or SignalParameterKeys

Returns:
metadata(Dict[str, Optional[str]): parsed metadata where datatype_keys are the keys of 
the dict
"""
def parse_datatype(lines, datatype_keys):
    metadata: Dict[str, Optional[str]] = {}
    for line in lines:
        if len(line.strip()) == 0:
            continue
        key, value = line.split("\t")
        key = key.strip()
        value = value.strip()
        
        if key in datatype_keys:
            metadata[key] = value
    return metadata

def parse_injection_info(lines):
    return parse_datatype(lines, InjectionKeys)

def parse_chromatogram_data_info(lines):
    return parse_datatype(lines, ChromatogramDataKeys)

def parse_signal_parameter_info(lines):
    return parse_datatype(lines, SignalParameterKeys)

"""
Parses raw chromatogram data and returns an np.array with the columns:
    [time, step, value]

This is an Nx3 array
"""
def parse_raw_chromatogram_data(rows):
    # Initialize empty np.array
    time_step_value = np.zeros((len(rows), 3))
    
    # Iterate through the rows and parse the data
    for i, row in enumerate(rows):
        cols = row.split('\t')
        
        # Parse the time
        time_str = cols[0].strip()
        if time_str == 'n.a.':
            time_val = np.nan
        else:
            time_val = float(time_str)
        
        # Parse the step
        step_str = cols[1].strip()
        if step_str == 'n.a.':
            step_val = np.nan
        else:
            step_val = float(step_str)

        # Parse the value (EU)
        value_str = cols[2].strip()
        if value_str == 'n.a.':
            value_val = np.nan
        else:
            value_val = float(value_str)
        
        time_step_value[i,:] = time_val, step_val, value_val
    return time_step_value
    
"""
potentially use dataclass to better control value types e.g.
    injection_info = InjectionInfo(
        data_vault=metadata["Data Vault"],
        injection=metadata["Injection"],
        injection_number=int(metadata["Injection Number"]),
        position=metadata.get("Position", None),
        comment=metadata.get("Comment", None),
        processing_method=metadata.get("Processing Method", None),
        instrument_method=metadata.get("Instrument Method", None),
        type=metadata.get("Type", None),
        status=metadata.get("Status", None),
        injection_date=metadata.get("Injection Date", None),
        injection_time=metadata.get("Injection Time", None),
        injection_volume=float(metadata.get("Injection Volume (µL)", None)),
        dilution_factor=float(metadata.get("Dilution Factor", None)),
        weight=float(metadata.get("Weight", None)))
"""
@dataclass
class InjectionInfo:
    data_vault: str
    injection: str
    injection_number: int
    position: str = None
    comment: str = None
    processing_method: str = None
    instrument_method: str = None
    type: str = None 
    status: str = None
    injection_date: str = None 
    injection_time: str = None 
    injection_volume: float = None
    dilution_factor: float = None
    weight: float = None

@dataclass
class ChromatogramDataInfo:
    time_min: float 
    time_max: float 
    data_points: int 
    detector: str = None 
    generating_data_system: str = None
    exporting_data_system: str =  None
    operator: str = None 
    signal_unit: str = None 
    signal_min: float = None 
    signal_max: float = None 
    channel: str = None
    driver_name: str = None 
    channel_type: str = None 
    min_step: float = None 
    max_step: float = None 
    average_step: float = None