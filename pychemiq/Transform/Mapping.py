# -*- coding: UTF-8 -*-
 
"""

============================

    @author       : Deping Huang
    @mail address : hdp@originqc.com
    @project      : pychemiq
    @source file  : Mapping.py

============================
"""

from pychemiq import (
    JordanWignerTransform,
    BravyiKitaevTransform,
    ParityTransform,
    MappingType
)

from pychemiq import fermion2pauli as Transform

def jordan_wigner(fermion):
    """
    Jordan-Wigner Transform
    """
    return JordanWignerTransform(fermion)

def bravyi_kitaev(fermion):
    """
    Bravyi-Kitaev Transform
    """
    return BravyiKitaevTransform(fermion)

def parity(fermion):
    """
    Parity Transform
    """
    return ParityTransform(fermion)

def segment_parity(fermion):
    """
    Segment-Parity Transform
    """
    return Transform(fermion,MappingType.SegmentParity)
