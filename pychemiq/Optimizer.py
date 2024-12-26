# -*- coding: UTF-8 -*-
 
"""

============================

    @author       : Deping Huang
    @mail address : hdp@originqc.com
    @date         : 2022-11-01 16:31:28
    @project      : pychemiq
    @source file  : Optimizer.py

============================
"""


from pychemiq import OptimizerFactory
import numpy as np
import sys

DEF_NELDER_MEAD      = "Nelder-Mead"
DEF_POWELL           = "Powell"
DEF_COBYLA           = "COBYLA"
DEF_LBFGSB           = "L-BFGS-B"
DEF_SLSQP            = "SLSQP"
DEF_GRADIENT_DESCENT = "Gradient-Descent"
DEF_LEARNING_RATE    = "learning_rate"


def vqe_solver(
            method=DEF_NELDER_MEAD, 
            ansatz=None,
            pauli=None,
            init_para = None,
            chemiq=None,
            pauli_group="none",
            Learning_rate = 0.1,
            Xatol=1e-4,
            Fatol=1e-4,
            MaxFCalls=200,
            MaxIter=200):
    """
    Docstrings for method vqe_solver
    """

    err = ""
    if chemiq == None:
        err = "ERROR: chemiq is needed!!!"
    if init_para is None:
        err = "ERROR: init_para is needed!!!"
    if pauli== None:
        err = "ERROR: pauli is needed!!!" 
    if ansatz == None:
        err = "ERROR: ansatz is needed!!!"
    if err != "":
        raise ValueError(err)
    
    method = method.upper().replace("_","-")
    if method == "NELDER-MEAD":
        method = DEF_NELDER_MEAD 
    elif method == "POWELL":
        method = DEF_POWELL
    elif method == "GRADIENT-DESCENT":
        method = DEF_GRADIENT_DESCENT

    chemiq.setOptimizerType(method)
    chemiq.setOptimizerIterNum(MaxIter)
    chemiq.setOptimizerFuncCallNum(MaxFCalls)
    chemiq.setOptimizerXatol(Xatol)
    chemiq.setOptimizerFatol(Fatol)
    chemiq.setPauliGroup(pauli_group)

    result = chemiq.getOptimizedData(0,init_para,pauli,chemiq.qvec,ansatz)
    
    return result 
