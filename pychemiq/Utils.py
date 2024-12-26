# -*- coding: UTF-8 -*-
 
"""

============================

    @author       : Deping Huang
    @mail address : hdp@originqc.com
    @date         : 2023-02-23 18:28:59
    @project      : pychemiq
    @source file  : Utils.py

============================
"""

from pychemiq import (
    transCC2UCC,
    FermionOperator,
    load_fermion_from_string,
    load_pauli_from_string,
    get_PES_geoms
)

def get_cc_n_term(n_qubits, n_elec,excited_level):   # get_ccsd_n_term接口
    '''
    get the number of coupled cluster single and double(ccsd) terms.
    e.g. 4 qubits, 2 electrons
    then 0 and 1 are occupied,just consider 0->2,0->3,1->2,1->3,01->23
    '''
    if n_elec > n_qubits:
        err = "Qubit num is less than electron num!"
        raise ValueError(err)

    n_diff = n_qubits - n_elec
    res = 0
    if "S" in excited_level:
        res += n_diff * n_elec
    if "D" in excited_level:
        res += n_diff*(n_diff - 1) * n_elec * (n_elec - 1) // 4

    return  res

f_op1 = lambda i,j:str(i) + "+ " + str(j)
f_op2 = lambda i,j,k,l:str(i)+"+ "+str(j)+"+ "+str(k)+" "+str(l)
def get_cc(n_qubits, n_elec, para,excited_level = "SD"):
    '''
    get Coupled cluster single and double terms with parameters(amplitudes before coupled cluster).
    e.g. 4 qubits, 2 electrons
    then 0 and 1 are occupied,just consider 0->2,0->3,1->2,1->3,01->23.
    returned FermionOperator like this:
    { {"2+ 0":para[0]},{"3+ 0":para[1]},{"2+ 1":para[2]},{"3+ 1":para[3]},
    {"3+ 2+ 1 0":para[4]} }
    '''
    if n_elec>n_qubits:
        err = "n_elec is bigger than n_qubits"
        raise ValueError(err)

    if n_elec==n_qubits:
        return FermionOperator()
    if get_cc_n_term(n_qubits, n_elec,excited_level) != len(para):
        err = "parameter number mismatched"
        raise ValueError(err)

    cnt = 0

    fermion_op = FermionOperator()
    if "S" in excited_level:
        for i in range(n_elec):
            for ex in range(n_elec, n_qubits):
                fermion_op += FermionOperator(f_op1(ex,i), para[cnt])
                cnt += 1

    if "D" in excited_level:
        for i in range(n_elec):
            for j in range(i+1,n_elec):
                for ex1 in range(n_elec,n_qubits):
                    for ex2 in range(ex1+1,n_qubits):
                        fermion_op += FermionOperator(f_op2(ex2,ex1,j,i),para[cnt])
                        cnt +=1
    return fermion_op

def get_valid_string(input_str):
    """
    "   abcd" -> "abcd"
    """
    if type(input_str) != type([]):
        input_str = input_str.split("\n")
    lines = []
    for line in input_str:
        line = line.strip()
        if line != "":
            lines.append(line)
    return "\n".join(lines)

def load_hamiltonian_from_string(op_str):
    """
    convert fermion string or pauli string into an operator class
    such as:
        "3+ 2+ 1 0" -> FermionOperator 
        "X0 Y1 Z2 Z3" -> PauliOperator
    """
    op_str = get_valid_string(op_str)
    if "+" in op_str:
        return load_fermion_from_string(op_str)
    else:
        return load_pauli_from_string(op_str)