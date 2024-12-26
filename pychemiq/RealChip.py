# -*- coding: UTF-8 -*-
 
"""

============================

    @author       : Deping Huang
    @mail address : hdp@originqc.com
    @date         : 2024-04-29 15:02:28
    @project      : pychemiq
    @source file  : RealChip.py

============================
"""

def real_chip_measure(
            ansatz=None,
            pauli=None,
            chemiq=None,
            task_id=None,
            api_key=None,
            init_para=None,
            mode = "wait",
            chip_id="72",
            shots=1000,
            cloud_url="https://pyqanda-admin.qpanda.cn",
            amend=True,
            mapping=True,
            circuit_opt=True
            ):
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
    if api_key == None:
        err = "ERROR: api_key is needed!!!"
    if mode == "query" and task_id==None:
        err = "ERROR: task_id is needed for the query mode!!!"
    if mode not in ["submit","query","wait"]:
        err = "ERROR: mode should be submit,query or wait!!!"
    if err != "":
        raise ValueError(err)

    chemiq.set_shots(shots)
    chemiq.setRealChip(mode,chip_id,amend,mapping,circuit_opt)
    #chemiq.init_cloud_machine(cloud_url,api_key)
    
    if mode in ["submit","wait"]:
        task_id = chemiq.submit_chip_task(ansatz,pauli,init_para)
        print("Your task id: ",task_id)
        if mode == "submit":
            return task_id
    if mode in ["query","wait"]:
        energy = chemiq.get_energy_by_wait(pauli,task_id)
        return energy
