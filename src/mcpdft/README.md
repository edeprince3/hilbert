<!--<p align="center">
<img src="logo.png" style='height: 30%; width: 50%; object-fit: contain'/> 
<br>
<a href="https://opensource.org/licenses/GPL-2.0"><img src="https://img.shields.io/github/license/edeprince3/v2rdm_casscf.svg" /></a>
</p>-->
<!-- <img src="logo.png" alt="RDM-inoles" align="middle"> -->


# mcpdft

mcpdft is a reduced-density-matrix (RDM)-based electronic structure plugin to Psi4. Given suitable input RDMs (that presumably captures static correlation effects), mcpdft can provide an efficient and accurate description of dynamic correlation effects that may be lacking in method used to generate the input RDMs. This software can also be used to carry out density-based analyses using post-Hartree-Fock methods.

## OVERVIEW

The multiconfiguration pair-density functional theory (MCPDFT) provides an accurate and efficient way of modeling static and dynamic correlation effects at low cost. Both translated and fully-translated versions of Slater and Vosko-Wilk-Nusair random-phase approximation expression III (SVWN3), Perdew-Burke-Ernzerhof (PBE), revised PBE (revPBE), Becke88 exchange and one-parameter correlation functional (BOP) and Becke and Lee-Yang-Parr (BLYP) on-top pair-density exchange-correlation functionals are available at the moment. In addition, the global-, double- and range-separated hybrid multi-configurational OTPDs such as wPBE and LRC-wPBE have also been implemented.

<!--
In summary, RDM-INOLES:

* can provide an interface with any (multiconfigurational) method that is able to provide 1-electron and 2-electron RDMs.
* hosts the variational 2-RDM driven complete active-space self-consistent field (v2RDM-CASSCF) as the reference method [2] by default
* can generate a .wfn file for further analysis of the wavefunction based on the quantum theory of atoms in molecules (QTAIMs)
* uses the reference total density and on-top pair-density (OTPD) functions as the input to build the so-called OTPD exchange-correlation (XC) functionals [1]
* features a double-hybrid MCPDFT method that is based on the linearly-scaled one-parameter double-hybrid (LS1DH) of Toulouse et al. described in Ref [3] 
* will provide and support both scaled and unscaled densities in MCPDFT
* 
-->

## INSTALLATION

To run the Psi4 plugin mcpdft:

* Download Psi4 from github.com: https://github.com/psi4/psi4, check out version 1.9.1 (git checkout f53cdd7), and follow the installation instructions given here: http://psicode.org/psi4manual/master/build_planning.html . Make sure to keep the name of the plugin directory mcpdft .

*  Configure with CMake to generate a Makefile. Run `psi4 --plugin-compile` to get a CMake command. Modify it as needed with `-D` for compiler, libraries, and options.

* Note that, if you configured Psi4 with a fortran compiler, you shouldn't have to specify these things here. If the configure shows no errors, compile the plugin:

  > make

* We recommend the intel compilers (icc, icpc and ifort) for the above procedures.

## INPUT OPTIONS

* **MCPDFT_METHOD** (string):

    The type of the multi-configurational on-top pair-density functional theory adopted for the calculation.
    The legitimate values include MCPDFT, 1H_MCPDFT, 1DH_MCPDFT, RS_MCPDFT, RS1H_MCPDFT, RS1DH_MCPDFT, and LS1DH_MCPDFT. The default
    value is MCPDFT.

* **MCPDFT_REFERENCE** (string):

    The type of reference 1- and 2-RDMs provided for MCPDFT calculations. The allowed values are V2RDM and CI. The
    default is v2RDM.

* **MCPDFT_TYPE** (string):

    The algorithm used for the MCPDFT computation. The valid options are DF and PK. Default is DF.    

* **MCPDFT_FUNCTIONAL** (string):

    The on-top pair-density exchange-correlation functionals used for MCPDFT. The valid options are SVWN, PBE, REVPBE
    BOP, BLYP, WPBE, LRC_WPBE. The default value is SVWN.

* **MCPDFT_TRANSLATION** (string):

    The density transformation scheme for MCPDFT. The valid options are REGULAR and FULL. The default is REGULAR.

* **MCPDFT_OMEGA** (double):

    The range-separation parameter. The default value is 0.0.

* **MCPDFT_LAMBDA** (double):

    The global-hybrid coupling parameter for hybrid functionals. The default value is 0.0.

* **WRITE_QTAIM_WFN** (bool):

    Writes the QTAIM .wfn file (AIMPAC and its successors). The default is false.

<!--
## REFERENCES

[1] M. Mostafanejad and A. E. DePrince III, J. Chem. Theory Comput. 15, 290-302 (2019). "Combining Pair-Density Functional Theory and Variational Two-Electron Reduced-Density Matrix Methods"

[2] J. Fosso-Tande, T.-S. Nguyen, G. Gidofalvi, and A. E. DePrince III, J. Chem. Theory Comput., 12, 2260-2271 (2016). "Large-scale variational two-electron reduced-density-matrix-driven complete active space self-consistent field methods."

[3] J. Toulouse, K. Sharkas, E. Bremond and C. Adamo J. Chem. Phys. 135, 101102 (2011). "Rationale for a new class of double-hybrid approximations in density-functional theory"
-->

