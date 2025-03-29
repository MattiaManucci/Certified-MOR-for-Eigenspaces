[![arXiv][arxiv-shield]][arxiv-url]
[![DOI][doi-shield]][doi-url]

# [Certified MOR for Eigenspaces][arxiv-url]

 We deal with the efficient and certified numerical approximation of the smallest eigenvalue and the associated eigenspace of a large-scale parametric Hermitian matrix. We rely on projection-based model order reduction (MOR), i.e., we approximate the large-scale problem by projecting it onto a suitable subspace and reducing it to one of a much smaller dimension. 
    Such a subspace is constructed by means of weak greedy-type strategies. We introduce a novel error estimate for the approximation error related to the eigenspace associated with the smallest eigenvalue. 
    Since the difference between the second smallest and the smallest eigenvalue, the so-called spectral gap, is crucial for the reliability of the error estimate, we propose efficiently computable upper and lower bounds for higher eigenvalues and for the spectral gap, which enable the assembly of a subspace for the MOR approximation of the spectral gap. 
    
This folder contains the Matlab codes that reproduce the numerical experiments in the paper [Certified Model Order Reduction for parametric Hermitian eigenproblems][arxiv-url].

##Code info:

* **Test\_Examples**: the main script to run to reproduce the numerical experiments of [preprint][arxiv-url].

<span style="color:blue">**Algorithm 1 Func.:**</span>

* **approx\_smallesteig\_all**: corresponds to Algorithm 1 in [preprint][arxiv-url] when solving the optimization problem with [EigOpt][Ref2];
* **approx\_sGAP\_Dset**: corresponds to Algorithm 1 in [preprint][arxiv-url] when solving the optimization problem for a discrete set $\Xi\subseteq \mathcal{D}$;
* **lambda\_eta\_eps**: function used to verify the subspace exact dimension recovery condition for the spectral gap approximation, see lines 11-16 from Algorithm 1 in [preprint][arxiv-url];
* **approximation\_SG**: function called from **approx\_smallesteig\_all** to evaluate the spectral gap surrogate error and its derivatives required by [EigOpt][Ref2]; see Section 4.1 from [preprint][arxiv-url]; 
* **approximation\_SG\_d**: function called from **approx\_sGAP\_Dset** to evaluate the spectral gap surrogate error; see Section 4.1 from [preprint][arxiv-url]; 


<span style="color:blue">**Algorithm 2 Func.:**</span>

* **App\_SEVES**: corresponds to Algorithm 2 in [preprint][arxiv-url] when solving the optimization problem with [EigOpt][Ref2];
* **subspace\_SCMM**: corresponds to Algorithm 2 in [preprint][arxiv-url] when solving the optimization problem for a discrete set $\Xi\subseteq \mathcal{D}$;
* **lambda\_eta\_eps\_eig**: function used to verify the subspace exact dimension recovery condition for the eigenspace approximation, see lines 11-16 from Algorithm 2 in [preprint][arxiv-url];
* **EVALUATE\_ERROR\_ESTIMATE**: function called from **App\_SEVES** to evaluate the eigenspace error estimate and its derivatives required by [EigOpt][Ref2]; see Section 4.2 from [preprint][arxiv-url]; 
* **Error\_Estimate\_EigVec\_Disc**: function called from **subspace\_SCMM** to evaluate the eigenspace error estimate; see Section 4.2 from [preprint][arxiv-url];
* **EVALUATE\_ERROR\_ESTIMATE\_C**: the same as **Error\_Estimate\_EigVec\_Disc** but return all the error estimate components (residual, eigenvalue error, and spectral gap); 

<span style="color:blue">**eigopt:**</span> the folder containing the functions that perform the optimization over the continuum domain $\mathcal{D}$, the software is developed by the authors of [EigOpt][Ref2];

<span style="color:blue">**Plot\_Functions:**</span> a folder that contains functions to visualise the numerical results. In particular: **plot_Spectral_GAP** reproduces the spectral gap related plots from [preprint][arxiv-url] and **plot_lambdamin_MD** reproduces the eigenspace related plots from [preprint][arxiv-url];

<span style="color:blue">**Data_Matrices:**</span> a folder that contains the matrices of the Quantum Spin Systems test problems used in [preprint][arxiv-url], e.g. test example 2 and 3;

<span style="color:red">**NOTE:**</span> download the test matrices for Example 1 [here][link-drive].


## Citing
If you use this project for academic work, please consider citing our
[publication][arxiv-url]:

    M. Manucci, B. Stamm, Z. Zeng
    Certified Model Order Reduction for parametric Hermitian eigenproblems
    ArXiv e-print ..., 2025.
    
## License
Distributed under the MIT License. See `LICENSE` for more information.

## Other references

* [Subspace Acceleration for Large-Scale Parameter-Dependent Hermitian Eigenproblems][Ref1]
* [Numerical Optimization of Eigenvalues of Hermitian Matrix Functions][Ref2]
* [Uniform Approximation of Eigenproblems of a Large-Scale Parameter-Dependent Hermitian Matrix][Ref3]

## Contacts

* Mattia Manucci - [mattia.manucci@simtech.uni-stuttgart.de](mattia.manucci@simtech.uni-stuttgart.de)



[doi-shield]: https://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.13254480-blue.svg?style=for-the-badge
[doi-url]: ....
[link-drive]: https://drive.google.com/file/d/1Y-RDkTQOvLaeccjgYzoq_qsD4tQLqN5u/view?usp=sharing
[arxiv-shield]: https://img.shields.io/badge/arXiv-2204.13474-b31b1b.svg?style=for-the-badge
[arxiv-url]: ....

[Ref1]: https://doi.org/10.1137/15M1017181
[Ref2]: https://doi.org/10.1137/130933472
[Ref3]: https://arxiv.org/abs/2409.05791
