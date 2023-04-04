# SANISAND-MS - a cyclic model for sand with ratcheting control

## About
This code is a FORTRAN implementation of the SANISAND-MS model proposed by Liu et al. (2019) - see citation below. 
The core research that led to the development of the model was conducted at TU Delft (Netherlands), in collaboration with Universidad des los Andes (Chile) and University of Bristol (UK).
SANISAND-MS enables advanced modelling of sand's cyclic behaviour, with emphasis on realistic simulation of ratcheting phenomena (cyclic accumulation of permanent strains).
The model was originally conceived to study the tilting response of offshore monopile foundations subjected to lateral cyclic loading, as is testified by its first 3D FE modelling application (Liu et al, 2022).
In its present form, the code may be directly used as UDSM (user-defined soil model) for PLAXIS 3D, for which the code has been extensively tested by the developers.
Important notes:
> - the code has been extensively tested with regard to the 3D FE modelling of cyclically loaded monopiles. Applications to other geotechnical problems are certainly possible, but caution is recommended;
> - compared to the original work of Liu et al. (2019), some model equations have been slightly modified to improve robustness in 3D FE calculations; a complete list of the equations implemented and some information about the UDSM are available in the section ['UDSM overview'](#UDSM-overview) of this  file;
> - users are welcome to report their comments and observation on the performance of the code to dr. Federico Pisanò (federicopisano83@gmail.com).


## Authors
The authors of the code are:
- Federico Pisanò (project leader and scientific supervisor)
- Zheng Li (code developer)
- Haoyuan Liu (code developer)
 

## Citation
If you use our code, please cite us! The main citations are:

> Original model formulation: 
Liu, H. Y., Abell, J. A., Diambra, A., & Pisanò, F. (2019). Modelling the cyclic ratcheting of sands through memory-enhanced bounding surface plasticity. 
Géotechnique, 69(9), 783-800. DOI: [10.1680/jgeot.17.P.307](https://doi.org/10.1680/jgeot.17.P.307)

> First application to 3D FE modelling of cyclically loaded monopiles: 
Liu, H., Kementzetzidis, E., Abell, J. A., & Pisanò, F. (2022). From cyclic sand ratcheting to tilt accumulation of offshore monopiles: 3D FE modelling using SANISAND-MS. 
Géotechnique, 72(9), 753-768. [DOI: 10.1680/jgeot.20.P.029](https://doi.org/10.1680/jgeot.20.P.029)


## License
This code is license under the GPL v3.0, see the LICENSE file.

Bentley Systems Inc. and affiliates were granted the exclusive right to use it under the terms of the [Unlicense](https://unlicense.org/). 


##  UDSM overview

The constitutive equations implemented in the UDSM subroutine (see table below) are largely based on the formulation described by Liu et al. (2019). Slight modifications have been implemented to improve computational robustness in 3D FE analyses. The same table also reports the 16 constitutive parameters introduced by Liu et al. (2019) and some of the variables included in the StVar array. The second table summarizes information related to the integration of the constitutive equations. It should be noted that, in addition to all material parameters, users should also provide the values of initial ($e_{ini}$), maximum ($e_{max}$) and minimum ($e_{min}$) void ratio. A complete description of the Props and StVar arrays is provided in SANISANMS.f90 and USRADDD.f90. More general information about the implementation of UDSMs in PLAXIS is available in the Material models manual. 

<table>
<thead>
<tr>
<td align="center">  </td>
<td align="center"> Implemented equation  </td>
<td align="center"> Material parameters </td>
<td align="center"> State/internal variables </td>
</tr>
</thead>
<tbody>
<tr>
<td> Elasticity </td>
<td>

```math
G = G_0 \ p_{atm} \frac{(2.97-e)^2} {(1+e)} \left( \frac{p} {p_{atm}} \right)^{1/2}
```

</td>
<td>

```math
G_0
```

</td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>

```math
K = \frac{2 \ (1+\nu) \ G} {3 \ (1-2 \ \nu)}
```

</td>
<td>

```math
\nu
```

</td>
<td>  </td>
</tr>

<tr>
<td> Critical state surface </td>
<td>

```math
\boldsymbol{r}_\theta^c = \sqrt{\frac{2}{3}} \ g(\theta) \ M_c \ \boldsymbol{n}
```

</td>
<td>

```math
M_c
```

</td>
<td>

```math
```

</tr>

<tr>
<td>  </td>
<td>

```math
g(\theta) = \frac{\alpha} {(1+\beta \ \cos 3 \theta)^{1/n}}
```

</td>
<td>

```math
c
```

</td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>

```math
\alpha = 2^{1/n} \ (1 + c^{-n})^{-1/n} \ \ \ \text{and} \ \ \ \beta = \frac{1-c^{-n}}{1+c^{-n} \ \ } \ \ \ \text{with} \ \ \ \ n = 4
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td> Critical state line </td>
<td>

```math
e_c = e_0 - \lambda_c \left( \frac{p}{p_{atm}}\right)^{\xi}
```

</td>
<td>

```math
e_0, \ \lambda_c, \ \xi
```

</td>
<td>  </td>
</tr>

<tr>
<td> Yield surface </td>
<td>

```math
f = \sqrt{(\boldsymbol{s} - p \ \boldsymbol{\alpha}):(\boldsymbol{s} - p \ \boldsymbol{\mathrm{\alpha}})} - \sqrt{\frac{2}{3}} \ p \ m = 0
```

</td>
<td>

```math
m
```

</td>
<td>

```math
\boldsymbol{\alpha}
```

</td>
</tr>

<tr>
<td> Memory surface </td>
<td>

```math
f^M = \sqrt{\left(\boldsymbol{s} - p \ \boldsymbol{\alpha}^M\right):\left(\boldsymbol{s} - p \ \boldsymbol{\alpha}^M\right)} - \sqrt{\frac{2}{3}} \ p \ m^M = 0
```

</td>
<td>  </td>
<td>

```math
\boldsymbol{\alpha}^M, \ \ m^M
```

</td>
</tr>

<tr>
<td> Yield surface hardening </td>
<td>

```math
d \boldsymbol{\alpha} = \frac{2}{3} \ \langle L \rangle \ h \ \left(\boldsymbol{\alpha}_\theta^b-\boldsymbol{\alpha}\right)\ \ \ \ \text{where} \ \ \ L: \ \text{plastic multiplier}
```

</td>
<td>

```math
```

</td>
<td>

```math
\boldsymbol{\alpha}
```

</td>
</tr>


<tr>
<td>  </td>
<td>

```math
\boldsymbol{\alpha}_\theta^b = \sqrt{\frac{2}{3}} \ \left[g(\theta) \ M \ \exp \left(-n^b \psi\right) - m \right] \ \boldsymbol{n} = \sqrt{\frac{2}{3}} \ \left[ M^b \ - m \right] \ \boldsymbol{n}
```

</td>
<td>

```math
n^b
```

</td>
<td>

```math
```

</tr>

<tr>
<td>  </td>
<td>

```math
\boldsymbol{n} = \frac{\boldsymbol{r} - \boldsymbol{\alpha}} { \left\| \boldsymbol{r} - \boldsymbol{\alpha} \right\| }
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>

```math
\psi = e - e_{cs}
```

</td>
<td>  </td>
<td>

```math
e, \ \psi
```

</tr>

<tr>
<td>  </td>
<td>

```math
h=\frac{b_0}{\left\|\left(\boldsymbol{\alpha}-{\boldsymbol{\alpha}_{i n}}\right):\boldsymbol{n}\right\|+C_{\varepsilon}}\ \exp\left[\mu_0\left(\frac{p}{p_{atm}}\right)^{0.5}\left(\frac{b^M}{b_{ref}}\right)^2\right]
```

</td>
<td>

```math
\mu_0
```

</td>
<td>

```math
\boldsymbol{\alpha}_{i n}
```

</tr>

<tr>
<td>  </td>
<td>

```math
C_{\varepsilon}\propto h_0/100 \ \ \  C_{\varepsilon} = 0.001 \ \ \ \ \ \ \text{see Boulanger and Ziotopoulou (2017)} 
```

</td>
<td>


</td>
<td>


</tr>

<tr>
<td>  </td>
<td>

```math
b_0=G_0 h_0\left(1-c_h e\right) / \sqrt{\left(p / p_{atm}\right)}
```

</td>
<td>

```math
h_0, c_h
```

</td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>

```math
b_{ref} = \left(\boldsymbol{\alpha}_\theta^b - \boldsymbol{\alpha}_{\theta + \pi}^b \right): \boldsymbol{n}
```

</td>
<td>  </td>
<td>  </td>
</tr>


<tr>
<td>  </td>
<td>

```math
\boldsymbol{\alpha}_{\theta + \pi}^b = \sqrt{\frac{2}{3}} \ \left[g(\theta + \pi) \ M \ \exp \left(-n^b \psi\right) - m \right] \ \boldsymbol{n} 
```

</td>
<td>


</td>
<td>


</tr>


<tr>
<td>  </td>
<td>

```math
b^M = \left({\boldsymbol{r}_\alpha^M} - \boldsymbol{\alpha}\right): \boldsymbol{n}
```

</td>
<td>

```math
```
</td>
<td>

```math
\boldsymbol{\alpha}
```
</tr>


<tr>
<td>  </td>
<td>

```math
{\boldsymbol{r}_\alpha^M} = \boldsymbol{\alpha}^M + \sqrt{\frac{2}{3}} \ \left(m^M-m\right) \ \boldsymbol{n}
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td> Volumetric plastic flow </td>
<td>

```math
d \varepsilon_{vol}^p = \left\langle L \right\rangle \ D \ \ \ \ \text{where} \ \ \ L: \ \text{plastic multiplier}
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>


```math
D = A_0 \ \exp \left( \frac{\beta \ \left\langle {\tilde{b}_d^M} \right\rangle} {b_{ref}} \right) \left( {\boldsymbol{\alpha}_\theta^d} - \boldsymbol{\alpha} \right) : \boldsymbol{n}
```

</td>
<td>

```math
A_0, \ \beta
```
</td>
<td>

```math
\boldsymbol{\alpha}
```
</tr>



<tr>
<td>  </td>
<td>

```math
\boldsymbol{\alpha}_\theta^d=\sqrt{2/3} \ \left[g(\theta) \ M \ \exp \left(n^d \psi\right) -m \right] \ \boldsymbol{n} = \sqrt{2/3} \ \left[M^d -m \right] \ \boldsymbol{n}
```

</td>
<td>

```math
n^d
```
</td>
<td>

```math
```
</tr>

<tr>
<td>  </td>
<td>

```math
{\tilde{b}_d^M} = \left( \tilde{{\boldsymbol{\alpha}}}_{\theta}^{\mathrm{d}}-\boldsymbol{\alpha}_{\mathrm{in}}\right): \boldsymbol{n}
```

</td>
<td>

```math
```
</td>
<td>

```math
\boldsymbol{\alpha}_{\mathrm{in}}
```
</tr>

<tr>
<td>  </td>
<td>

```math
b_{ref} = \left(\boldsymbol{\alpha}_\theta^b - \boldsymbol{\alpha}_{\theta + \pi}^b \right): \boldsymbol{n}
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td> Deviatoric plastic flow </td>
<td>

```math
d\boldsymbol{e}^p=\langle L\rangle \ \boldsymbol{R}^{\prime}\ \ \ \ \text{where} \ \ \ L: \ \text{plastic multiplier}
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>

```math
\boldsymbol{R}^{\prime} = B \ \boldsymbol{n} - C \left[\boldsymbol{n}^2 - \frac{1}{3} \ \boldsymbol{I}\right]
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>

```math
B = 1 + \frac{3}{2} \ \frac{1-c}{c} \ g(\theta) \ \cos 3 \theta
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>

```math
C = 3 \ \sqrt{\frac{3}{2}} \ \frac{1-c}{c} \ g(\theta)
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td> Memory surface hardening </td>
<td>

```math
d\boldsymbol{\alpha}^M = \frac{2}{3} \ \left\langle L\right\rangle \ h^M \ \left({\boldsymbol{\alpha}_\theta^b} - {\boldsymbol{r}_\alpha^M} \right)\ \ \ \ \text{where} \ \ \ L: \ \text{plastic multiplier}
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>

```math
h^M = \frac{1}{2} \ \left[\frac{b_0} {\left\| \left( {{\boldsymbol{r}_\alpha^M}} - {\boldsymbol{\alpha}_{i n}} \right): \boldsymbol{n} \right\| +C_{\varepsilon} } + \sqrt{\frac{3}{2}} \ \frac{m^M \ f_{shr}\ \langle -D \rangle} {\zeta \ \left( {\boldsymbol{\alpha}_\theta^b} - {\boldsymbol{r}_\alpha^M} \right): \boldsymbol{n}} \right]
```

</td>
<td>

```math
```
</td>
<td>

```math
\boldsymbol{\alpha}_{i n}
```
</tr>


<tr>
<td>  </td>
<td>

```math
\text{or} \ \ h^M = \frac{1}{2} \ \left[ \frac{b_0} {\left\| \left( {{\boldsymbol{r}_\alpha^M}} - {\boldsymbol{\alpha}_{i n}} \right) : \boldsymbol{n} \right\| + C_{\varepsilon}} + \sqrt{\frac{3}{2}} \ \frac{m^M \ f_{shr} \ \langle -D \rangle} {\zeta \ \left( 1 + H  \left( {\boldsymbol{n}^{\mathrm{bM}}} : \boldsymbol{n} \right) \right)} \right] \operatorname{sgn} \left[ \boldsymbol{n}^{\mathrm{bM}} : \boldsymbol{n}\right]
```

</td>
<td>

```math
```
</td>
<td>

```math
\boldsymbol{\alpha}_{i n}
```
</tr>


<tr>
<td>  </td>
<td>

```math
\boldsymbol{n}^{bM} = \frac{ {\boldsymbol{\alpha}_\theta^b} - {\boldsymbol{r}_\alpha^M} } { \left\| {\boldsymbol{\alpha}_\theta^b} - {\boldsymbol{r}_\alpha^M} \right\| } \ ,  \ \  \operatorname{or}  \ \ \boldsymbol{n}^{bM} = \boldsymbol{n}
```

</td>
<td>  </td>
<td>  </td>
</tr>

<tr>
<td>  </td>
<td>

```math
d m^M = \sqrt{\frac{3}{2}} \ d \boldsymbol{\alpha}^M : \boldsymbol{n} - \frac{m^M} {\zeta} \ f_{shr} \ \left\langle-d \varepsilon_{vol}^p\right\rangle
```

</td>
<td>

```math
\zeta
```

</td>
<td>

```math
```

</tr>
</tbody>
</table>



|              | Description           | Parameters          | 
| :----------- | :----------- | :----------- | 
| yield surface tolerance:   | 	The yield surface correction scheme described by Sloan et al. (2001) is adopted. The search of the stress at the yield surface intersection point and the elastoplastic unloading also uses TolF | $$TolF$$
| Error Tolerance:   | The integration of the elastoplastic strain increment is performed with the explicit RKF23 scheme with automatic substepping and error control| $$TolR$$
| Toggle:   |  If =0, the integration of the elastic strain increments is performed using single step Forward Euler scheme. If =1, the elastic trial stress is computed with an explict RKF23 scheme using the same TolR set for the elastoplastic portion. If = 2 all the purely elastic strain increments are computed with RKF23 | $$El Scheme$$
| Toggle:   | If =0, $h^M$ is set to a constant value in the exceptions. If =1, the new $h^M$ equation for tackling corner cases (shown in table 1) is used | $$h^M Exc Eq$$
| Toggle:   | If =0, default criterion is used to limit $h$. If =value, $h$ maximum is defined in input.  | $$h_{Max}$$




> Sloan, S. W., Abbo, A. J. & Sheng, D. (2001). Refined explicit integration of elastoplastic models with automatic error control. 
Engineering computations, Vol. 18 No. 1/2, 121-154. https://doi.org/10.1108/02644400110365842

> Boulanger, R., Ziotopoulou, K. (2017). PM4Sand (Version 3.1). A sand plasticity model for earthquake engineering applications.
Report No. UCD/CGM-17/01. Revised July 2018 University of California, Davis

