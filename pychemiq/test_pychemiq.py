import pytest
import numpy as np
from Molecules import Molecules
from Optimizer import vqe_solver
from Circuit.Ansatz import (
    UCC,
    HardwareEfficient,
    SymmetryPreserved,
    UserDefine
)
from pychemiq.Utils import (
    get_cc_n_term,
    get_cc
)
from . import QMachineType,ChemiQ
from Transform.Mapping import (
    jordan_wigner,
    bravyi_kitaev,
    parity,
    segment_parity,
    MappingType,
    Transform
)

class Test_Molecules:
    def setup_method(self):
        self.geom = [
            "H 0 0 0",
            "H 0 0 0.74"
        ]
        self.spin = 1
        self.charge = 0
        self.basis =  "sto-3g"

    def testMolecules(self):
        """
        docstring for testMolecules
        """
        mol = Molecules(
            geometry = self.geom,
            basis    = self.basis,
            multiplicity = self.spin,
            charge = self.charge
        )
    
        print(mol.get_molecular_hamiltonian())
        #mol.get_active_space_integrals()
        
        return
    
    def test_geom(self):
        mol = Molecules(
            basis    = self.basis,
            multiplicity = self.spin,
            charge = self.charge
        )
        return
    
    def test_basis(self):
        mol = Molecules(
            geometry = self.geom,
            multiplicity = self.spin,
            charge = self.charge
        )
        return
    
    def test_spin(self):
        mol = Molecules(
            geometry = self.geom,
            basis    = self.basis,
            charge = self.charge
        )
        return
    
    def test_nfrozen(self):
        mol = Molecules(
            geometry = self.geom,
            basis    = self.basis,
            nfrozen = 1,
            charge = self.charge
        )
        mol = Molecules(
            geometry = self.geom,
            basis    = self.basis,
            nfrozen = 0,
            charge = self.charge
        )
        return

class Test_Optimizer:
    def setup_method(self):
        self.geom = [
            "H 0 0 0",
            "H 0 0 0.74"
        ]
        self.spin = 1
        self.charge = 0
        self.basis =  "sto-3g"
        self.mol = Molecules(
            geometry = self.geom,
            basis    = self.basis,
            multiplicity = self.spin,
            charge = self.charge
        )

        fermion = self.mol.get_molecular_hamiltonian()
        self.mapping_type = MappingType.Bravyi_Kitaev
        self.pauli = Transform(fermion,self.mapping_type)
        self.n_qubits = self.pauli.get_max_index() + 1

        self.chemiq = ChemiQ()
        machine_type = QMachineType.CPU_SINGLE_THREAD
        pauli_size = len(self.pauli.data())
        n_elec = self.mol.n_electrons
        self.chemiq.prepare_vqe(machine_type,self.mapping_type,n_elec,pauli_size,self.n_qubits)
        self.ansatz = UCC("UCCSD",self.mol.n_electrons,self.mapping_type,chemiq=self.chemiq)
        self.init_para = np.zeros(self.ansatz.get_para_num())

        print("n_qubits = ",self.n_qubits)

    def testOptimizer(self):

        method = "NELDER_MEAD"
        result = vqe_solver(
                method = method,
                pauli = self.pauli,
                ansatz = self.ansatz,
                chemiq=self.chemiq,
                init_para=self.init_para)
        print(result)
        return

    def test_chemiq(self):

        method = "NELDER_MEAD"
        try:
            result = vqe_solver(
                    method = method,
                    pauli = self.pauli,
                    ansatz = self.ansatz,
                    init_para=self.init_para)
        except ValueError as e:
            print(e)
        return

    def test_init_para(self):

        method = "NELDER_MEAD"
        try:
            init_para = [None]*self.ansatz.get_para_num()
            init_para = np.array(init_para)
            result = vqe_solver(
                    method = method,
                    pauli = self.pauli,
                    ansatz = self.ansatz,
                    chemiq=self.chemiq,
                    init_para=init_para)
        except ValueError as e:
            print(e)
        return

    def test_pauli(self):

        method = "NELDER_MEAD"
        try:
            result = vqe_solver(
                    method = method,
                    ansatz = self.ansatz,
                    chemiq=self.chemiq,
                    init_para=self.init_para)
        except ValueError as e:
            print(e)
        return

    def test_ansatz(self):

        method = "NELDER_MEAD"
        try:
            result = vqe_solver(
                    method = method,
                    pauli = self.pauli,
                    chemiq=self.chemiq,
                    init_para=self.init_para)
        except ValueError as e:
            print(e)
        return

    def loss(self,para,grad,iters,fcalls):
        res = self.chemiq.getLossFuncValue(0,para,grad,iters,fcalls,self.pauli,self.chemiq.qvec,self.ansatz)
        return res

    def test_powell(self):

        method = "POWELL"
        result = vqe_solver(
                method = method,
                pauli = self.pauli,
                ansatz = self.ansatz,
                chemiq=self.chemiq,
                init_para=self.init_para)
        print(result)
        return

    def test_lbfgsb(self):

        method = "L_BFGS_B"
        result = vqe_solver(
                method = method,
                pauli = self.pauli,
                ansatz = self.ansatz,
                chemiq=self.chemiq,
                init_para=self.init_para)
        print(result)
        return

class Test_Ansatz:
    def setup_method(self):
        self.geom = [
            "H 0 0 0",
            "H 0 0 0.74"
        ]
        self.spin = 1
        self.charge = 0
        self.basis =  "sto-3g"
        self.mol = Molecules(
            geometry = self.geom,
            basis    = self.basis,
            multiplicity = self.spin,
            charge = self.charge
        )

        fermion = self.mol.get_molecular_hamiltonian()
        self.mapping_type = MappingType.Bravyi_Kitaev
        self.pauli = Transform(fermion,self.mapping_type)
        self.n_qubits = self.pauli.get_max_index() + 1

        self.chemiq = ChemiQ()
        machine_type = QMachineType.CPU_SINGLE_THREAD
        pauli_size = len(self.pauli.data())
        n_elec = self.mol.n_electrons
        self.chemiq.prepare_vqe(machine_type,self.mapping_type,n_elec,pauli_size,self.n_qubits)

    def test_get_cc_n_term(self):
        for level in ["S","D","SD"]:
            print(get_cc_n_term(4,2,level))
        try:
            print(get_cc_n_term(4,5,"SD"))
        except ValueError as e:
            print(e)
        return

    def test_get_cc(self):

        for level in ["S","D","SD"]:
            n = get_cc_n_term(4,2,level)
            para = [0]*n
            fermion = get_cc(4,2,para,excited_level=level)

        fermion = get_cc(4,4,[])

        try:
            fermion = get_cc(4,5,[])
        except ValueError as e:
            print(e)

        try:
            fermion = get_cc(4,2,[])
        except ValueError as e:
            print(e)

        return 

    def test_UCC(self):
        ansatz = UCC("UCCSD",self.mol.n_electrons,self.mapping_type,chemiq=self.chemiq)
        ansatz = UCC("UCCS",self.mol.n_electrons,self.mapping_type,chemiq=self.chemiq)
        ansatz = UCC("UCCD",self.mol.n_electrons,self.mapping_type,chemiq=self.chemiq)

        self.mapping_type = MappingType.Jordan_Wigner
        ansatz = UCC("UCCSD",self.mol.n_electrons,self.mapping_type,chemiq=self.chemiq)
        self.mapping_type = MappingType.Parity
        ansatz = UCC("UCCSD",self.mol.n_electrons,self.mapping_type,chemiq=self.chemiq)

        try:
            ansatz = UCC("UCCSD",self.mol.n_electrons,self.mapping_type)
        except ValueError as e:
            print(e)

        try:
            ansatz = UCC("test",self.mol.n_electrons,self.mapping_type,chemiq=self.chemiq)
        except ValueError as e:
            print(e)

        try:
            ansatz = UCC("UCCSD",self.mol.n_electrons,"test",chemiq=self.chemiq)
        except ValueError as e:
            print(e)

        return

    def test_HardwareEfficient(self):
        ansatz = HardwareEfficient(2,chemiq=self.chemiq)

        try:
            ansatz = HardwareEfficient(2)
        except ValueError as e:
            print(e)
        return

    def test_SymmetryPreserved(self):
        ansatz = SymmetryPreserved(2,chemiq=self.chemiq)

        try:
            ansatz = SymmetryPreserved(2)
        except ValueError as e:
            print(e)
        return

    def test_UserDefine(self):
        try:
            ansatz = UserDefine(2,chemiq=self.chemiq)
        except ValueError as e:
            print(e)

        try:
            ansatz = UserDefine(2)
        except ValueError as e:
            print(e)

        circuit = "H q[0]"
        try:
            ansatz = UserDefine(2,circuit = circuit)
        except ValueError as e:
            print(e)

        # try:
        #     ansatz = UserDefine(2,circuit = circuit,option=None)
        # except ValueError as e:
        #     print(e)

        fermion = "3+ 2+ 1 0: 0"
        try:
            ansatz = UserDefine(2,fermion = fermion)
        except ValueError as e:
            print(e)
        ansatz = UserDefine(2,circuit = circuit,chemiq=self.chemiq)

        return

class Test_Mapping:
    def setup_method(self):
        self.geom = [
            "H 0 0 0",
            "H 0 0 0.74"
        ]
        self.spin = 1
        self.charge = 0
        self.basis =  "sto-3g"
        self.mol = Molecules(
            geometry = self.geom,
            basis    = self.basis,
            multiplicity = self.spin,
            charge = self.charge
        )

        self.fermion = self.mol.get_molecular_hamiltonian()
    def testMapping(self):
        pauli = jordan_wigner(self.fermion)
        pauli = bravyi_kitaev(self.fermion)
        pauli = parity(self.fermion)
        pauli = segment_parity(self.fermion)


if __name__ == "__main__":
    pytest.main("-s test_pychemiq.py")
