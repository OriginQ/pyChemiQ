## pyChemiQ

pyChemiQ是一款由本源量子开发的python软件库, 它基于ChemiQ封装的python接口，可以用于在真实量子计算机或虚拟机上实现费米子体系的模拟与计算。支持用户自定义映射、拟设、优化器，方便用户二次开发。该软件包为量子化学计算和方法开发提供了一个简单、轻量且高效的平台。 pyChemiQ 可用于在量子计算机上使用平均场和后平均场方法模拟分子来解决量子化学问题。 pyChemiQ简化了分子结构输入和量子电路之间的转换，最大限度地减少了进入该领域所需的领域专业知识，让感兴趣的学者更方便解决和研究量子计算机上的电子结构问题。


- pyChemiQ详细使用文档及代码示例详见：[pyChemiQ-tutorial](https://pychemiq-tutorial.readthedocs.io/en/latest/index.html)

- 体验全部最新功能的可视化教学工具ChemiQ，请前往 [官网](https://qcloud.originqc.com.cn/zh/chemistryIntroduce) 下载。

- 若您想在文献中引用pyChemiQ或ChemiQ, 请按照如下格式: Wang Q, Liu H Y, Li Q S, et al. Chemiq: A chemistry simulator for quantum computer[J]. arXiv preprint arXiv:2106.10162, 2021.


## 安装pyChemiQ

我们提供了Linux,macOS,Windows上的python预编译包供安装，目前pyChemiQ支持3.8、3.9、3.10版本的python。

```python
pip install pychemiq
```


## 计算氢分子单电能代码示例

```python
from pychemiq import Molecules,ChemiQ,QMachineType
from pychemiq.Transform.Mapping import (jordan_wigner,MappingType)
from pychemiq.Optimizer import vqe_solver
from pychemiq.Circuit.Ansatz import UCC
import numpy as np

multiplicity = 1
charge = 0
basis =  "sto-3g"
geom = "H 0 0 0,H 0 0 0.74"

mol = Molecules(
    geometry = geom,
    basis    = basis,
    multiplicity = multiplicity,
    charge = charge)

fermion_H2 = mol.get_molecular_hamiltonian()
pauli_H2 = jordan_wigner(fermion_H2)


chemiq = ChemiQ()
machine_type = QMachineType.CPU_SINGLE_THREAD
mapping_type = MappingType.Jordan_Wigner
pauli_size = len(pauli_H2.data())
n_qubits = mol.n_qubits
n_elec = mol.n_electrons
chemiq.prepare_vqe(machine_type,mapping_type,n_elec,pauli_size,n_qubits)
ansatz = UCC("UCCSD",n_elec,mapping_type,chemiq=chemiq)

method = "SLSQP"
init_para = np.zeros(ansatz.get_para_num())
solver = vqe_solver(
        method = method,
        pauli = pauli_H2,
        chemiq = chemiq,
        ansatz = ansatz,
        init_para=init_para)
result = solver.fun_val
n_calls = solver.fcalls
print(result)
```
