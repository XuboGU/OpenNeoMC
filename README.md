# OpenNeoMC：an Open-source Tool for Particle Transport Optimization that Combining OpenMC with NEORL

[OpenMC](https://docs.openmc.org/en/stable/index.html) is a community-developed Monte Carlo neutron and photon transport simulation code for particle transport. OpenMC was originally developed by members of the [Computational Reactor Physics Group](http://crpg.mit.edu/) at the [Massachusetts Institute of Technology](https://web.mit.edu/) starting in 2011.

[NEORL](https://neorl.readthedocs.io/en/latest/index.html) (**N**euro**E**volution **O**ptimization with **R**einforcement **L**earning) is a set of implementations of hybrid algorithms combining neural networks and evolutionary computation based on a wide range of machine learning and evolutionary intelligence architectures. NEORL aims to solve large-scale optimization problems relevant to operation & optimization research, engineering, business, and other disciplines. NEORL was established in MIT back in 2020 with feedback, validation, and usage of different colleagues. 

In OpenNeoMC, we combine these two open-source tools to empower particle transport with state-of-the-art optimization techniques. We firstly provide users with easy ways to install the framework that combines NEORL with OpenMC, and a simple example is available to test the framework. Then we offer two practical engineering optimization applications in nuclear physics. More applications that involve both optimization and nuclear physics will be added in the future. We highly welcome users and researchers in the nuclear area to contribute OpenNeoMC and solve engineering problems in this framework.

## Installing OpenNeoMC

### Installation on Linux/Mac with conda

#### Install Conda

Please install conda before proceeding,  it will bring you convenience to install [anaconda](https://www.anaconda.com/products/individual#Downloads) directly, which includes conda and other necessary python  packages.

Create a new virtual environment named openneomc with python3.7 (python3.7 is requied for compatibility purposes)

```shell
conda create -n openneomc python=3.7

# activate the virtual environment 
conda activate openneomc  
```
### Install NEORL

```shell
pip install neorl
```

Check if you have install NEORL successfully by input the following command in the terminal.

```SHELL
neorl
```

If you see the 'NEORL' logo, then you have prepared the OpenNeoMC framework, congratulations! 

Note: you could also have a comprehensive test of NEORL by running the following command (this is time-consuming).

```shell
neorl --test
```

#### Install OpenMC of Version 0.12.2 
Refer to [OpenMC installation](https://docs.openmc.org/en/v0.13.1/quickinstall.html#installing-from-source-on-linux-or-mac-os-x) 

```shell
git clone --branch v0.12.2 https://github.com/openmc-dev/openmc.git
cd openmc
mkdir build && cd build
cmake ..
make
sudo make install
cd .. && pip install .
```

**Cross Section Configuration**
You may encounter the [no cross_sections.xml error]( https://github.com/openmc-dev/openmc/issues/1855) when running OpenMC. This is caused by the missing of nuclear data, you could solve it refer to [Cross Section Configuration](https://docs.openmc.org/en/stable/usersguide/cross_sections.html)

**Download cross section data**
Various cross section data are available on the  [OpenMC official website](https://openmc.org/data-libraries/), from the OpenMC team, LANL, etc. In OpenNeoMC, we use [ENDF/B-VII.1](https://openmc.org/official-data-libraries/) in default. But if you have specific purpose, please feel free to use other nuclear data that you need. 

After downloading the cross-section data file, configure it as an environmental variable as follows. 

**Add environmental variables**

```shell
## Temporary methods
# in python
import os
os.environ['OPENMC_CROSS_SECTIONS'] = '/PATH/cross_sections.xml'
# in shell
export OPENMC_CROSS_SECTIONS=../cross_sections.xml

## Once for all: you can modify the ~/.bashrc to configure environmental variables
# open ~/.bashrc
vim ~/.bashrc
# add the following command in the end 
export OPENMC_CROSS_SECTIONS=/PATH/cross_sections.xml
# update 
source ~/.bashrc
```


#### Test OpenNeoMC

Let's test OpenNeoMC by the  'pin_cell_test.py' example (Remember to configure environmental variables as above). 

```shell
# run 
python pin_cell_test.py
```

If you see the 'NEORL' logo and the log information of OpenMC, then congratulations, you could start your journey of OpenNeoMC! 

### Installing OpenNeoMC with Docker on Linux/Mac/Windows 

Installing OpenNeoMC with docker is highly recommended! In this way, you need not worry about issues like cross-section data and software compatibility, etc. All you need to do are simply pull the image and run it in your own machine with any OS.  

#### Install Docker

Follow the official tutorial to Install docker on your machine:   [get docker](https://docs.docker.com/get-docker/)

#### Install OpenNeoMC

After installing docker, your can easily install use OpenNeoMC framework within only four steps: 

```shell
# Pull docker images from dock hub  
sudo docker pull 489368492/openneomc

# Check the openmc docker images
sudo docker images

# Run the openmc images to create container named `openneomc`
sudo docker run -tid --shm-size=8G --gpus all --name openneomc -v /LocalWorkingDir/:/workspace/ 489368492/openneomc

# Execute the container
sudo docker exec -it openneomc /bin/bash
```

Note: in `docker run` step, the `-v` flag mounts the current working directory into the container, which is very convenient for users. 

Please refer to  [Docker CLI](https://docs.docker.com/engine/reference/commandline/run/) for docker command-line descriptions.

**Other commonly used commands** 

```shell
# Exit the container
exit

# Stop the container
sudo docker stop openneomc

# Start the container
sudo docker start openneomc

# Delete the container
sudo docker rm openneomc

# Delete the image(remove the container first)
sudo docker image rm 489368492/openneomc
```

#### Test OpenNeoMC 

Let's test OpenNeoMC by the  'pin_cell_test.py' example, which can be found at the path '/home'

```shell
# cd /home
cd /home

# run 
python pin_cell_test.py
```

If you see the 'NEORL' logo and the log information of OpenMC, then congratulations, you have successfully run OpenNeoMC in docker way.

The program runs around 3 minutes(may vary depending on your CPU), and the results are like:

```python
------------------------ JAYA Summary --------------------------
Best fitness (y) found: 0.0015497217274231812
Best individual (x) found: [2.01355604]
--------------------------------------------------------------
---JAYA Results---
x: [2.01355604]
y: 0.0015497217274231812
JAYA History:
 [0.018311916874464318, 0.0017114252626817539, 0.0017114252626817539, 0.0017114252626817539, 0.0015497217274231812]
running time:
 155.2281835079193
```

## Reference

OpenMC: https://docs.openmc.org/en/stable

OpenMC docker image: https://hub.docker.com/r/openmc/openmc

NEORL: https://neorl.readthedocs.io/en/latest/

OpenNeoMC docker image: https://hub.docker.com/r/489368492/openneomc 

## Citation
If you use this project in your research, please cite it using the following Bibtex entry:
```bibtex
@article{GU2023109450,
  title = {OpenNeoMC: A framework for design optimization in particle transport simulations based on OpenMC and NEORL},
  journal = {Annals of Nuclear Energy},
  volume = {180},
  pages = {109450},
  year = {2023},
  issn = {0306-4549},
  doi = {https://doi.org/10.1016/j.anucene.2022.109450},
  url = {https://www.sciencedirect.com/science/article/pii/S0306454922004807},
  author = {Xubo Gu and Majdi I. Radaideh and Jingang Liang},
  keywords = {OpenMC, NEORL, Particle transport, Optimization, Criticality search},
}
```

## Contact 

If you have any suggestions or issues, please feel free to contact Xubo Gu(xbgu1024@126.com)

