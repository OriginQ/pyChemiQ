from setuptools import setup, find_packages
import platform
import os 

pd_files = []
if platform.system() == 'Windows':
    pd_files += ['*.pyd', '*.dll', '*.pyi']
else :
    pd_files += ['*.so', '*.pyi']

setup(
    name = 'pychemiq', 
    version='1.1.2', 
    license = "Apache Licence",  
    author = "OriginQ",
    include_package_data = True,  
    packages = [
        "pychemiq",
        "pychemiq.Transform",
        "pychemiq.Circuit",
        "pychemiq.basis",
    ],
    package_dir = {"":os.path.dirname(__file__)},
    package_data = {
        '':pd_files,
        "pychemiq.basis":["*.g94"]
    },
)
