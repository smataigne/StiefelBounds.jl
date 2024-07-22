# StiefelBounds.jl
This repository is associated to the publication "Mataigne S., Absil P. -A. Absil, Miolane N., Bounds on the geodesic distances on the Stiefel manifold for a family of Riemannian metrics, 2024". This article contains the theorems and corollary which are mentioned in this file.
We make available the code necessary to produce the numerical experiments presented in the figures 2, 3 and 7 from the article. 

## Requirements
This repository requires a Julia installation with the following packages: `LinearAlgebra`, `SkewLinearAlgebra` for computational aspects and `Plots`,  `LaTeXStrings` for the plots. These packages are easily obtained from the package installation environment as follows. In Julia REPL, press `]` to access the installation environment and for each package, do
```julia
(@v1.6) pkg> add Name_of_Package
```
## Reproduction of the experiments
Before running any other file, run `Manifold.jl` from the `src` folder in local. `src` contains the algorithms to compute minimal geodesics on the Stiefel manifold $\mathrm{St}(n,p)$. 
The files to produce the figures are in the `make_plots` folder under the names `figure2.jl`, `figure3.jl` and `figure7.jl`. In these files, the size of $\mathrm{St}(n,p)$ can be tuned, as well as the number of random samples $N$. In `figure2.jl`, $n\geq2p$ is a requirement of the algorithm to compute minimal geodesics.
## Bilipschitz equivalence of the $\beta$-distances. 
![My nice plot](https://github.com/smataigne/StiefelBounds.jl/blob/main/make_plots/figures/figure_2.png)

Bounds between the canonical and the Euclidean distance on the Stiefel manifold $\mathrm{St}(4,2)$: $\frac{\sqrt{2}}{2} d_{\mathrm{E}}(U,\widetilde{U})\leq d_{\mathrm{c}}(U,\widetilde{U})\leq d_{\mathrm{E}}(U,\widetilde{U})$ for all $U,\widetilde{U}\in \mathrm{St}(4,2)$ (see Corollary 3.4). The plot shows the canonical distances $d_\mathrm{c}$ ($\beta = \frac{1}{2}$) and the Euclidean distances $d_\mathrm{E}$ ($\beta=1$) for 10000 randomly generated pairs $(U,\widetilde{U})$ on $\mathrm{St}(4,2)$.

## Upper and lower bounds on the geodesic distance by the Frobenius distance
![My nice plot](https://github.com/smataigne/StiefelBounds.jl/blob/main/make_plots/figures/figure_3.png)

Bounds between the Frobenius and the Euclidean geodesic distance on the Stiefel manifold $\mathrm{St}(5,4)$: $2\sqrt{p}\arcsin\left(\frac{\Vert U-\widetilde{U}\Vert_\mathrm{F}}{2\sqrt{p}}\right)\leq d_{\mathrm{E}}(U,\widetilde{U})\leq \frac{\pi}{2}\Vert U-\widetilde{U}\Vert_\mathrm{F}$ for all $U,\widetilde{U}\in \mathrm{St}(5,4)$ (see Theorems 5.3 and 6.4). The plot shows the Frobenius distances and the Euclidean distances $d_\mathrm{E}$ ($\beta=1$) for 10000 randomly generated pairs $(U,\widetilde{U})$ on $\mathrm{St}(5.4)$.

## Minimal geodesics reaching the upper bound
![My nice plot](https://github.com/smataigne/StiefelBounds.jl/blob/main/make_plots/figures/figure_10.png)

Evolution of the Euclidean geodesic distance ($\beta=1$) in terms of the Frobenius distance for  minimal geodesics going  from $U$ to $h_k(U)$ for $k=1, 2, 3, 4$ in $\mathrm{St}(5,4)$ (solid black lines). The upper bound, i.e., $d_\mathrm{E}(U,\widetilde{U})\leq \frac{\pi}{2}\Vert U-\widetilde{U}\Vert_\mathrm{F}$ for all $U,\widetilde{U}\in \mathrm{St}(5,4)$ (see Theorem 8.4), is represented by the green dashed line. Notice that the condition $n\geq 2p$ is not satisfied but Theorem 6.4 is verified nonetheless. 
