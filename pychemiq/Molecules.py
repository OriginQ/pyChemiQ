# -*- coding: UTF-8 -*-
 
"""

============================

    @author       : Deping Huang
    @mail address : hdp@originqc.com
    @date         : 2022-10-20 16:39:17
    @project      : pychemiq
    @source file  : Molecules.py

============================
"""

from pychemiq import (
    Mole,
    set_current_execute_path,
    set_log_levels,
    CCSD
)

import numpy as np
import platform,os

# to find the directory of basis
if "Windows" in platform.platform():
    sep = "\\"
else:
    sep = "/"
path  = os.path.dirname(__file__) + sep 
set_current_execute_path(path)
class Molecules(Mole):
    """
    Docstrings for class Molecules
    geometry: atom types and atom coordinates
        geometry = "H 0 0 0,H 0 0 0.74" or 
        geometry = [
            "H 0 0 0",
            "H 0 0 0.74"
        ]
    """
    def __init__(self, 
        geometry=None,
        basis=None,
        multiplicity=None,
        charge=0,
        active = None,
        nfrozen = None,
        filename="",
        diis = "cdiis",
        diis_n = 8,
        diis_thre = 0.01,
        hamiltonian_thre=1e-8,
        directory=None):

        if type(geometry) == type([]):
            geometry = ",".join(geometry)
        if active == None:
            active = [0,0]
        if nfrozen == None:
            nfrozen = 0
        self.active = active 
        self.nfrozen = nfrozen
        self.geometry     = geometry
        self.basis        = basis
        self.multiplicity = multiplicity
        self.charge       = charge
        if geometry == None:
            self.argPrintError("geometry")
            return 
        if basis == None:
            self.argPrintError("basis")
            return
        if multiplicity == None:
            self.argPrintError("multiplicity")
            return

        self.bohr = False
        self.pure = False
        # set log levels
        set_log_levels("",6,6)

        Mole.__init__(
            self,geometry,basis,charge,multiplicity,
            self.bohr,self.pure,"","rectangular"
        )

        # settings for diis 
        self.setDIIS(diis,diis_n,diis_thre)
        # setting active space
        if self.nfrozen != 0:
            self.active = [0,-1]
            self.setNfrozen(self.nfrozen)
        self.setActiveSpace(self.active)
        self.set_hamiltonian_thre(hamiltonian_thre)
        self.hf_iters = 1000 
        self.hf_threshold = 1e-10
        self.HF(self.hf_iters,self.hf_threshold) 
        self.molecular_hamiltonian = self.getHamiltonian()
        self.getAttributes()  
        self.ccsd = None 
          
        return
    def getAttributes(self):
        """
        Docstrings for method getAttributes
        """

        self.n_atoms = self.getnatom()
        self.n_electrons = self.getElectronNum()
        self.hf_energy = self.getEhf()
        self.nuclear_repulsion = self.getEnuc()
        self.canonical_orbitals = self.getMolecularOrbitals()
        self.n_orbitals = self.getnmo()
        self.n_qubits   = 2*self.getnmo()
        self.orbital_energies = self.getOrbitalEnergies()
        self.overlap_integrals = self.getOverlapMatrix()
        self.one_body_integrals = self.getA()
        self.two_body_integrals = self.getGFromMatrix()
        return

    def getGFromMatrix(self):
        """
        Docstrings for method getG
        """
        G = self.getDoubleIntegrals()
        # number of basis functions
        n = self.one_body_integrals.shape[0]
        G = G.reshape((n,n,n,n))
        return G 
    
    def argPrintError(self,err):
        """
        Docstrings for method printError
        """
        print(f"ERROR: Please input {err}!!!!")
        return
    def get_molecular_hamiltonian(self):
        """
        Docstrings for method get_molecular_hamiltonian
        """
        return self.molecular_hamiltonian
    
    def ccsd_init(self):
        """
        Docstrings for method ccsd_init
        """
        mol = CCSD(
            self.geometry,
            self.basis,
            self.charge,
            self.multiplicity,
            self.bohr,
            self.pure,
            "",
            "rectangular"
        )
        mol.run_scf(self.hf_threshold,self.hf_iters) 
        mol.setActiveSpace(self.active)
        return mol 
    def get_ccsd_energy(self):
        """
        Docstrings for method get_ccsd_energy
        """
        if self.ccsd == None:
            self.ccsd = self.ccsd_init()
            
        return self.ccsd.get_CCSD_energy()
    
    def get_ccsd_init_para(self,excited_level):
        """
        Docstrings for method get_ccsd_init_para
        """
        if excited_level not in ["D","SD"]:
            raise RuntimeError("only D or SD is available")
        if self.ccsd == None:
            self.ccsd = self.ccsd_init()
        return self.ccsd.get_init_para(excited_level,True)