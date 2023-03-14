## pyChemiQ简介

pyChemiQ是一款由本源量子开发的python软件库, 它基于ChemiQ封装的python接口，可以用于在真实量子计算机或虚拟机上实现费米子体系的模拟与计算。支持用户自定义映射、拟设、优化器，方便用户二次开发。该软件包为量子化学计算和方法开发提供了一个简单、轻量且高效的平台。 pyChemiQ 可用于在量子计算机上使用平均场和后平均场方法模拟分子来解决量子化学问题。 pyChemiQ简化了分子结构输入和量子电路之间的转换，最大限度地减少了进入该领域所需的领域专业知识，让感兴趣的学者更方便解决和研究量子计算机上的电子结构问题。

目前，pyChemiQ支持输入分子结构得到二次量子化后的Fermion Hamiltonian；映射方面，pyChemiQ支持Jordan-Wigner(JW)变换、Bravyi-Kitaev(BK)变换、Parity变换和Multilayer Segmented Parity(MSP)变换方法将二次量子化的Fermion Hamiltonian算符映射成量子计算机上的pauli Hamiltonian算符；拟设方面，pyChemiQ也支持不同的拟设构造量子电路，例如Unitary Coupled Cluster（UCC）、Hardware-Efficient、symmetry-preserved等拟设构造量子电路；优化器方面，pyChemiQ提供了以下几种经典优化器来进行变分参数的优化：NELDER-MEAD、POWELL、COBYLA、L-BFGS-B、SLSQP和Gradient-Descent。通过自定义映射和拟设方式，用户也可以自行构造或优化量子线路，求得电子结构问题的最优基态能量解。

- 欢迎您在Github上Fork我们的项目或给我们留言：[pyChemiQ 项目](https://github.com/OriginQ/pyChemiQ.git)

- 为特定量子化学问题定制pyChemiQ功能，请联系本源量子陈先生18221003869或给我们发送邮件 [dqa@originqc.com](mailto:dqa@originqc.com)。

- 体验全部最新功能的可视化教学工具ChemiQ，请前往 [官网](https://qcloud.originqc.com.cn/zh/chemistryIntroduce) 下载。

- 若您想在文献中引用pyChemiQ或ChemiQ, 请按照如下格式: Wang Q, Liu H Y, Li Q S, et al. Chemiq: A chemistry simulator for quantum computer[J]. arXiv preprint arXiv:2106.10162, 2021.

  
## 安装pyChemiQ

我们提供了Linux,Windows上的python预编译包供安装，pyChemiQ要求python版本为3.8。除此之外，pyChemiQ仅依赖numpy函数包（版本大于1.19）。若使用Linux版本安装，需保证系统glibc版本大于等于2.29。

如果您已经安装好了python环境和pip工具，pyChemiQ可通过如下命令进行安装：

```python
pip install pychemiq
```

<!--若要运行test5Moles.py，您需要额外安装nlopt，scipy，tqdm 工具包。-->


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
print(result,f"函数共调用{n_calls}次")
```