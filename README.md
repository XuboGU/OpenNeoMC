## OpenNeoMC

[OpenMC](https://docs.openmc.org/en/stable/index.html) is a community-developed Monte Carlo neutron and photon transport simulation code for nuclear research. OpenMC was originally developed by members of the [Computational Reactor Physics Group](http://crpg.mit.edu/) at the [Massachusetts Institute of Technology](https://web.mit.edu/) starting in 2011

[NEORL](https://neorl.readthedocs.io/en/latest/index.html) (**N**euro**E**volution **O**ptimization with **R**einforcement **L**earning) is a set of implementations of hybrid algorithms combining neural networks and evolutionary computation based on a wide range of machine learning and evolutionary intelligence architectures. NEORL aims to solve large-scale optimization problems relevant to operation & optimization research, engineering, business, and other disciplines. NEORL was established in MIT back to 2020 with feedback, validation, and usage of different colleagues. 

In this project(OpenNeoMC), we dedicated to combine these two open source tools(OpenMC and NEORL) to empower nuclear physics with the state-of-the-art optimization techniques. We will firstly provide users with easy ways to install the framework that combine NEORL with OpenMC, and a simple example is available to test the framework. Then we offer two practical engineering optimization applications in nuclear physics. More applications that involve both optimization and nuclear physics will be added in the future. We highly welcome users and researchers in nuclear area to contribute OpenNeoMC and solve engineering problems in this framework.

## Installing on Linux/Mac with conda

### Install OpenMC

```shell
conda config --add channels conda-forge
conda search openmc
```

Create a new virtual environment named openneomc

```shell
conda create -n openneomc openmc
```

#### Test OpenMC

Follow with the [official examples]( https://docs.openmc.org/en/stable/examples/index.html) to test the OpenMC

**Cross Section Configuration**

You may encounter the [no cross_sections.xml error]( https://github.com/openmc-dev/openmc/issues/1855) when running OpenMC. This is caused by the missing of nuclear data, you could solve this by [Cross Section Configuration](https://docs.openmc.org/en/stable/usersguide/cross_sections.html)

**Download cross section data**

Various cross section data are available on the [official website](https://openmc.org/data-libraries/), from the OpenMC team, LANL etc. Download one version you need, and then configure it as an environmental variable as follows

**Add environmental variables**

```shell
## in python
import os
os.environ['OPENMC_CROSS_SECTIONS'] = '../cross_sections.xml'

## in shell
export OPENMC_CROSS_SECTIONS=../cross_sections.xml

## in shell you can modify the ~/.bashrc to configure environmental variables automatically
# open ~/.bashrc
vim ~/.bashrc
# add the following command in the end 
export OPENMC_CROSS_SECTIONS=../cross_sections.xml
# update 
source ~/.bashrc

```



### Install NEORL

Install python 3.7 to make sure the stable run of tensorflow-1.14.0 

```shell
conda install python=3.7 
```

```shell
pip install neorl==1.6
```

Change the scikit-learn version to 0.24.2

```shell
pip install scikit-learn==0.24.2
```

test NEORL

```shell
neorl
```

unit test of NEORL

```SHELL
neorl --test
```



## Installing with Docker on Linux/Mac/Windows 



## Contact 

If you have any suggestions or issues, please feel free to contact Xubo Gu(guxubo@alumni.sjtu.edu.cn)

