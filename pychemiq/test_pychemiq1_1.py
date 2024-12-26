from pychemiq import *
import pytest
import numpy as np

class Test_all():
    def setup_method(self):
        self.chemiq = ChemiQ()
        self.machine_type = QMachineType.CPU_SINGLE_THREAD
        self.pauli_size = 15
        self.n_elec     = 2
        self.n_qubits   = 4
        self.mapping_type = MappingType.Bravyi_Kitaev

        self.chemiq.prepare_vqe(
            self.machine_type,
            self.mapping_type,
            self.n_elec,
            self.pauli_size, 
            self.n_qubits
        )

        self.ansatz = UCC("UCCD",2,self.mapping_type,chemiq=self.chemiq)
        pauli_BK = [
            ": -0.097067",
            "X0 Z1 X2 : 0.045303",
            "X0 Z1 X2 Z3 : 0.045303",
            "Y0 Y2 : 0.045303",
            "Y0 Y2 Z3 : 0.045303",
            "Z0 : 0.171413",
            "Z0 Z1 : 0.171413",
            "Z0 Z1 Z2 : 0.120625",
            "Z0 Z1 Z2 Z3 : 0.120625",
            "Z0 Z2 : 0.165928",
            "Z0 Z2 Z3 : 0.165928",
            "Z1 : 0.168689",
            "Z1 Z2 : -0.223432",
            "Z1 Z3 : 0.174413",
            "Z2 Z3 : -0.223431"
        ]
        self.hamiltonian = load_hamiltonian_from_string(pauli_BK)

    def test_cdiis(self):
        """
        docstring for test_cdiis
        """
        mol = Molecules(
            geometry = "H 0 0 0,Li 0 0 1",
            basis    = "sto-3g",
            multiplicity = 1,
            charge = 0,
            diis   = "cdiis"
        )
        mol2 = Molecules(
            geometry = "H 0 0 0,Li 0 0 1",
            basis    = "sto-3g",
            multiplicity = 1,
            charge = 0,
            diis   = "none" 
        )
        delta = abs(mol.hf_energy - mol2.hf_energy)
        assert(delta < 1e-8)
        
        return

    def test_pauli_ansatz(self):
        """
        docstring for test_pauli_ansatz
        """
        pauli = [
            "X0 Y2 : -0.125000",
            "X0 Y2 Z3 : -0.125000",
            "X0 Z1 Y2 : -0.125000",
            "X0 Z1 Y2 Z3 : -0.125000",
            "Y0 X2 : 0.125000",
            "Y0 X2 Z3 : 0.125000",
            "Y0 Z1 X2 : 0.125000",
            "Y0 Z1 X2 Z3 : 0.125000",
        ]
        pauli = "\n".join(pauli)

        n_electrons = 2
        ansatz = UserDefine(n_electrons,pauli=pauli,chemiq=self.chemiq)
        
        return          
    def test_hamiltonian(self):
        """
        docstring for test_hamiltonian
        """
        fermion = [
            ": 0.715104",
            "0+ 0 : -1.253310",
            "1+ 0+ 1 0 : -0.674756",
            "1+ 0+ 3 2 : -0.181210",
            "1+ 1 : -1.253310",
            "2+ 0+ 2 0 : -0.482501",
            "2+ 1+ 2 1 : -0.663711",
            "2+ 1+ 3 0 : 0.181210",
            "2+ 2 : -0.475069",
            "3+ 0+ 2 1 : 0.181210",
            "3+ 0+ 3 0 : -0.663711",
            "3+ 1+ 3 1 : -0.482501",
            "3+ 2+ 1 0 : -0.181210",
            "3+ 2+ 3 2 : -0.697652",
            "3+ 3 : -0.475069"
        ]
        fermion = "\n".join(fermion)
        fermion = load_hamiltonian_from_string(fermion)
        pauli_BK = [
            ": -0.097067",
            "X0 Z1 X2 : 0.045303",
            "X0 Z1 X2 Z3 : 0.045303",
            "Y0 Y2 : 0.045303",
            "Y0 Y2 Z3 : 0.045303",
            "Z0 : 0.171413",
            "Z0 Z1 : 0.171413",
            "Z0 Z1 Z2 : 0.120625",
            "Z0 Z1 Z2 Z3 : 0.120625",
            "Z0 Z2 : 0.165928",
            "Z0 Z2 Z3 : 0.165928",
            "Z1 : 0.168689",
            "Z1 Z2 : -0.223432",
            "Z1 Z3 : 0.174413",
            "Z2 Z3 : -0.223431"
        ]
        pauli = load_hamiltonian_from_string(pauli_BK)
        return

    def test_parameter_matrix(self):
        """
        docstring for test_parameter_matrix
        """
        pauli = [
            "X0 Y2 : -0.125000",
            "X0 Y2 Z3 : -0.125000",
            "X0 Z1 Y2 : -0.125000",
            "X0 Z1 Y2 Z3 : -0.125000",
            "Y0 X2 : 0.125000",
            "Y0 X2 Z3 : 0.125000",
            "Y0 Z1 X2 : 0.125000",
            "Y0 Z1 X2 Z3 : 0.125000",
        ]

        n_electrons = 2

        maxtrix = [[-1]*4+[1]*4]
        matrix = [
            [0]*4 + [1]*4,
            [1]*4 + [0]*4
        ]

        ansatz = UserDefine(
            n_electrons,
            pauli=pauli,
            chemiq=self.chemiq,
            parameter_matrix = matrix
        )
        
        return

    def test_get_PES_geoms(self):
        configs = {
            "CH4": [
                "H 0.62527605 0.62527605 0.62527605",
                "C 0.00000000 0.00000000 0.00000000",
                "H -0.62527605 -0.62527605 0.62527605",
                "H -0.62527605 0.62527605 -0.62527605",
                "H 0.62527605 -0.62527605 -0.62527605"

            ],
            "NH3":[
                "H 0.00000000 0.97134700 -0.19948700",
                "N 0.00000000 0.00000000 0.08549400",
                "H -0.84121100 -0.48567300 -0.19948700",
                "H 0.84121100 -0.48567300 -0.19948700"
            ]
        }
        Atoms = [
            [1,2],
            [1,2,3],
            [1,2,3,4]
        ]
        values = [0.8,0.9,1.0,1.1,1.2]
        for mol in configs:
            geom = configs[mol]
            geom = ",".join(geom)
            for atoms in Atoms:
                geoms = get_PES_geoms(geom,atoms,values)
                print(mol,atoms,geoms)

    def test_ccsd(self):
        """
        Docstrings for method test_ccsd
        """
        configs = {
            "bf": [
                "B   -6.8504918   -0.4332537   -2.6667315",
                "F   -5.4062400   -0.4332537   -2.6667315"
            ],
            "c2h2": [
                "C       0.0000071193    -0.0000055702     0.6015557270",
                "C      -0.0000025711    -0.0000491954    -0.6015287172",
                "H       0.0000075581     0.0001505625     1.6676884008",
                "H      -0.0000262625     0.0001237641    -1.6676645547",
            ],
            "h2o": [
                "O 0 0 0.12713100",
                "H 0 0.75801600 -0.50852400",
                "H 0 -0.75801600 -0.50852400"
            ],
            "hcn": [
                "N       0.0000000003     0.0026684533    -0.6028601129",
                "C      -0.0000000000    -0.0000433415     0.5512272572",
                "H      -0.0000000003    -0.0026553507     1.6200066436"
            ],
            "n2": [
                "N 0 0 0",
                "N 0 0 1.1"
            ],
            "nh3": [
                "N 0 0 0.08549400",
                "H 0 0.97134700 -0.19948700",
                "H -0.84121100 -0.48567300 -0.19948700",
                "H 0.84121100 -0.48567300 -0.19948700"
            ],
            "o2": [
                "O 0 0 0",
                "O 0 0 1.2"
            ],
            "p2": [
                "P       0.0000000000     0.0000000000    -1.1010000000",
                "P       0.0000000000     0.0000000000     1.1010000000"
            ],
            "ph3": [
                "P      -0.0001960169    -0.0000074722    -0.3039389321",
                "H       0.5881618309    -1.0422334096     0.4591261393",
                "H       0.6086106494     1.0303352920     0.4591322030",
                "H      -1.1965720738     0.0119054801     0.4590158346"
            ],
            "s2": [
                "S 0 0 0",
                "S 0 0 1.902804136"
            ]
        }
        objects = {
            "bf": [
                -124.36283147749145,
                -0.2744210228435292
            ],
            "c2h2": [
                -77.10272300016271,
                -0.27697445600668547
            ],
            "h2o": [
                -76.23914807275743,
                -0.21600124143365515
            ],
            "hcn": [
                -93.18065900425863,
                -0.2974974944395519
            ],
            "n2": [
                -109.26738926696598,
                -0.31359302607510975
            ],
            "nh3": [
                -56.3972586360081,
                -0.20440198201696097
            ],
            "o2": [
                -149.92593354005515,
                -0.38171879145057697
            ],
            "p2": [
                -681.6881470494035,
                -0.2838355158276682
            ],
            "ph3": [
                -342.64113987505823,
                -0.17063613960831436
            ],
            "s2": [
                -795.3002735548919,
                -0.2928948666486043
            ]
        }

        self.basis =  "cc-pvdz"
        for key in configs:
            self.geom = configs[key]
            self.mol = Molecules(
                geometry = self.geom,
                basis    = self.basis,
                multiplicity = 1,
                charge = 0
            )
            self.mol.pure = True
            res= self.mol.get_ccsd_energy()
            delta = np.array(res) - np.array(objects[key])
            delta = np.average(np.abs(delta))
            assert(delta < 2e-7)
            
        return
    def test_ccsd_init_para(self):
        """
        Docstrings for method test_ccsd_init_para
        """
        self.mol = Molecules(
            geometry = "H 0 0 0,H 0 0 0.74",
            basis    = "sto-3g",
            multiplicity = 1,
            charge = 0
        )
        ccsd_para = self.mol.get_ccsd_init_para("D")
        return
    
    def test_pauli_group(self):
        """
        Docstrings for method test_pauli_group
        """
        method = "SLSQP"
        
        init_para = [0]
        res1 = vqe_solver(
            method = method,
            pauli = self.hamiltonian,
            ansatz = self.ansatz,
            chemiq = self.chemiq,
            init_para=init_para
        )
        res2 = vqe_solver(
            method = method,
            pauli = self.hamiltonian,
            ansatz = self.ansatz,
            chemiq = self.chemiq,
            init_para=init_para,
            pauli_group="naive"
        )
        delta = res1.fun_val - res2.fun_val
        assert(delta < 1e-8)
        return
    def test_real_chip(self):
        """
        Docstrings for method test_real_chip
        """
        set_log_levels("",0,6)
        api_key = "302e020100301006072a8648ce3d020106052b8104001c041730150201010410ce40b64bee8d3bf2ecbeb4e7c46422cd/14820"
        chemiq = ChemiQ()
        circuit = [
            "H q[0]",
            "RX q[2],(1.5707963)",
            "CNOT q[0],q[3]",
            "CNOT q[1],q[3]",
            "CNOT q[2],q[3]",
            "RZ q[3],(0.887217)",
            "CNOT q[0],q[3]",
            "CNOT q[1],q[3]",
            "CNOT q[2],q[3]",
            "DAGGER",
            "H q[0]",
            "RX q[2],(1.5707963)",
            "ENDDAGGER"
        ]
        chemiq.setQCloudMachineApiKey(api_key)
        chemiq.prepare_vqe(
            QMachineType.QCloud,
            self.mapping_type,
            self.n_elec,
            self.pauli_size, 
            self.n_qubits
        )
        ansatz = UserDefine(2,circuit=circuit,chemiq=chemiq)
        res = real_chip_measure(
            ansatz = ansatz,
            mode = "wait",
            pauli = self.hamiltonian,
            chemiq = chemiq,
            api_key = api_key,
            init_para = [-0.2256]
        )
        print(res)
        
        return
      


if __name__ == "__main__":
    pytest.main("-s test_pychemiq1_1.py")