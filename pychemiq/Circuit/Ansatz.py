# -*- coding: UTF-8 -*-
 
"""

============================

    @author       : Deping Huang
    @mail address : hdp@originqc.com
    @project      : pychemiq
    @source file  : Ansatz.py

============================
"""

from pychemiq import (
    AnsatzFactory
)

from pychemiq.Transform.Mapping import MappingType
from pychemiq import FermionOperator,PauliOperator
import numpy as np
from pychemiq.Transform.Mapping import bravyi_kitaev
from pychemiq.Utils import get_valid_string
import sys

DEF_RESTRICTED                    = "restricted"
DEF_CUTOFF                        = "cutoff"
DEF_MAPPING_METHOD                = "mapping method"
DEF_REORDER                       = "reorder"
DEF_HAMILTONIAN_SIMULATION_SLICES = "slices"
DEF_EVOLUTION_TIME                = "evolution_time"
DEF_EXCITED_LEVEL                 = "excited_level"
DEF_CIRCUIT                       = "circuit"

DEF_JORDAN_WIGNER                 = "Jordan-Wigner"
DEF_PARITY                        = "Parity"
DEF_BRAVYI_Kitaev                 = "Bravyi-Kitaev"
DEF_SEGMENT_PARITY                = "SegmentParity"

restricted = "T"
cutoff     = "F"
slices     = "1"
evo_time   = 1.0
reorder    = "F"

ansatz_options = {}
ansatz_options[DEF_RESTRICTED] =  restricted
ansatz_options[DEF_CUTOFF] =  cutoff
ansatz_options[DEF_REORDER] =  reorder
ansatz_options[DEF_HAMILTONIAN_SIMULATION_SLICES] =  slices
ansatz_options[DEF_EVOLUTION_TIME] = str(evo_time)


def UCC(ucc_type,n_electrons,mapping_type,chemiq=None):
    """
    ucc_type: 
        only "UCCS", "UCCD", "UCCSD" are available
    n_electrons:
        number of electrons
    chemiq:
        object of class ChemiQ
    """
    if chemiq == None:
        err = "parameter chemiq should be assigned"
        raise ValueError(err)

    mapping = ""

    if mapping_type == MappingType.Bravyi_Kitaev:
        mapping = DEF_BRAVYI_Kitaev
    elif mapping_type == MappingType.Jordan_Wigner:
        mapping = DEF_JORDAN_WIGNER  
    elif mapping_type == MappingType.Parity:
        mapping = DEF_PARITY 
    elif mapping_type == MappingType.SegmentParity:
        mapping = DEF_SEGMENT_PARITY
    else:
        err = "ERROR: no such mapping type!!!"
        raise ValueError(err)

    ansatz_options[DEF_MAPPING_METHOD] =  mapping 
    if ucc_type == "UCCS":
        ansatz_options[DEF_EXCITED_LEVEL] =  "S"
    elif ucc_type == "UCCD":
        ansatz_options[DEF_EXCITED_LEVEL] =  "D"
    elif ucc_type == "UCCSD":
        ansatz_options[DEF_EXCITED_LEVEL] =  "SD"
    else:
        err = "ERROR: invalid ucc type"
        raise ValueError(err)

    ansatz = AnsatzFactory.makeAnsatz("ucc",chemiq.qvec,n_electrons,ansatz_options)
    return ansatz

def HardwareEfficient(n_electrons,chemiq=None):
    if chemiq == None:
        err = "parameter chemiq should be assigned"
        raise ValueError(err)
    
    ansatz = AnsatzFactory.makeAnsatz("hardwareefficient",chemiq.qvec,n_electrons,ansatz_options)
    return ansatz

def SymmetryPreserved(n_electrons,chemiq=None):
    if chemiq == None:
        err = "parameter chemiq should be assigned"
        raise ValueError(err)
    
    ansatz = AnsatzFactory.makeAnsatz("symmetrypreserved",chemiq.qvec,n_electrons,ansatz_options)
    return ansatz

def UserDefine(n_electrons,
        circuit=None,
        fermion=None,
        pauli=None,
        chemiq=None,
        parameter_matrix=None,
        opt=None,
        mapping_type=MappingType.Jordan_Wigner):
    """
    circuit: 
        originir string of the circuit of the ansatz
    n_electrons:
        number of electrons
    n_qubits:
        number of qubits
    chemiq:
        object of class ChemiQ
    """

    if circuit == None and fermion==None and pauli==None:
        err = "no circuit or fermion or pauli is assigned"
        raise ValueError(err)

    #if option != None:
    #    for key in option:
    #        ansatz_options[key] = option[key]
    #    option = ansatz_options
    #else:
    option = ansatz_options
    if circuit != None:
        circuit = get_valid_string(circuit)
        option[DEF_CIRCUIT] =  circuit
    if fermion != None:
        if type(fermion) == type(FermionOperator()):
            fermion = fermion.to_string()
        else:
            fermion = get_valid_string(fermion)
        option["fermion"] =  fermion
    if pauli != None:
        if type(pauli) == type(PauliOperator()):
            pauli = pauli.to_string()
        else:
            pauli = get_valid_string(pauli)
        option["pauli"] =  pauli

    if chemiq == None:
        err = "parameter chemiq should be assigned"
        raise ValueError(err)

    if parameter_matrix != None:
        paras = []
        for line in parameter_matrix:
            item = " ".join(list(map(str,line)))
            paras.append(item)
        option["parameter_matrix"] = "\n".join(paras)

    # mapping type
    map_method = "mapping method"
    if mapping_type == MappingType.Jordan_Wigner:
        option[map_method] = "Jordan-Wigner"
    elif mapping_type == MappingType.Bravyi_Kitaev:
        option[map_method] = "Bravyi-Kitaev"
    elif mapping_type == MappingType.Parity:
        option[map_method] = "Parity"
    elif mapping_type == MappingType.SegmentParity:
        option[map_method] = "SegmentParity"
    else:
        raise RuntimeError("Error mapping type")

    ansatz = AnsatzFactory.makeAnsatz("userdefine",chemiq.qvec,n_electrons,option)
    return ansatz
