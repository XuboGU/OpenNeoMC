# Assembly's max k-eff - a discrete problem with constraint 

## Assembly Model

The assembly model refers from the reference paper [*Gradient Informed Design Optimization of Select Nuclear Systems*](https://www.tandfonline.com/doi/abs/10.1080/00295639.2021.1987133)[1].  

**Geometry & Materials**

The assembly model is an 11 X 11 voxelated square grid  with a total width of around 1 meter. Each voxel wad filled with material 1 (void) or material 2 (a fuel-moderator mixture). The fuel-moderator mixture is the traditional UN TRISO fuel embeded in a SiC shell, and the moderator is yttrium hydride[2]. 

The fuel-moderator materials amount is limited: for example, if there only have 57 units of fuel/moderator, then one can set at most 57 voxel as material 2. 

**Boundary condition**

* In X and Y directions, the boundary condition is set as 'void' (i.e.: 'vaccum' in openMC)
* In Z direction, the boundary condition is set as 'reflective' (same in openMC)

**Optimal geometric configuration**

The reason the authors set X and Y directions of the model as void boundary conditions is that the optimal geometric configuration can be predictable - a circle, wherein the center of geometry should be fuel pins while edging area the void pin. This circle geometry can be a benchmark. The following picture is an optimal geometic configuration for 47 fuel pins.

<center>    
    <img style="border-radius: 0.3125em;    
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"     
    src=".\pics\assm_max_paper.png">    
    <br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    
    display: inline-block;    
    color: #999;    
   	padding: 2px;">Optimal geometric configuration (47 fuel pins)[1]</div> 
</center>



## Assembly Model Optimization 

**Problem Definition**

The optimization objective of this case is to find the optimal geometric configuration to maximize k-eff using limited amount of fuel-moderator materials in this 11 X 11 geometry. To be more specific, one needs to determine which material (1 or 2) to fill in every voxel of the geometry. 

In this case, we use up to 61 fuel pins to fill the assembly and need to find the optimal geometric configuration to make the k-eff max. Equation (1) is the math definition of the problem. Where $(p_i = 0)$ represents the void pin and $(p_i = 1)$ represents the fuel pin.

$$
\begin{equation} 
\begin{aligned} 
max_{\vec p}f(\vec p), \\
subject \ to,\\ 
\sum_ip_i\leq61, \\
where \ p_i = 0 \ or \ p_i = 1 
\end{aligned}
\tag 1
\end{equation}
$$



**Possible combinations**

The possible combinations of filling the assembly can be calculated by eqution (2).
$$
\sum_{i=0}^{61} C_{121}^{i} \approx 1.521 \times 10^{36}  \tag 2
$$
This is a quite large number, thus, it's hard to find the optimal configuration. 

### Method

Here we call the [Defferential Evolution](https://neorl.readthedocs.io/en/latest/modules/de.html#differential-evolution-de) algorithm to find the optimal geometric configuration. 

The figure showes the converge curve of Differential Evolution. It showes that the algorithm find the best combination at around generation 380.  The best k-eff found is 1.03796.. 

<img src=".\pics\de_curve.png" alt="converge curve" style="zoom:80%;" />


The geometry of the best solution is shown in the following figure, where a good circle is found. This result coincides with the benchmark.

<center>    
    <img style="border-radius: 0.3125em;    
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"     
    src=".\pics\app2_geo_nonsym.jpg">    
    <br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    
    display: inline-block;    
    color: #999;    
   	padding: 2px;">Optimal geometric configuration found by DE (61 fuel pins)</div> 
</center>



### Leveraging the plane symmetry

The assembly is usually in symmetric structures, of which we can take advantage to reduce the search space. We applied 1/4 symmetry to the assembly and the problem can be scale down to 6 x 6. DE was again used to solve the discrete problem (400 generations with 50 individuals per generation).

<img src=".\pics\de_curve2.png" alt="converge curve" style="zoom:80%;" />

We can note that DE converges much faster in the 1/4 symmetry case . It only takes around 80 genetraions to find the best solution, the k_eff is 1.03092. 

The geometry of the best solution is shown in the following figure,  The result is a good circle and coincides with the benchmark.

<center>    
    <img style="border-radius: 0.3125em;    
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"     
    src=".\pics\app2_geo_14sym.jpg">    
    <br>    
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;    
    display: inline-block;    
    color: #999;    
   	padding: 2px;">Optimal geometric configuration found by DE (1/4 symmetry)</div> 
</center>



## Reference

[1] John Pevey, Briana Hiscox, Austin Williams, Ondřej Chvála, Vladimir Sobes & J. Wesley Hines (2021) Gradient-Informed Design Optimization of Select Nuclear Systems, Nuclear Science and Engineering

[2] Transformational Challenge Reactor Preconceptual Core Design Studie: https://www.osti.gov/servlets/purl/1651338

