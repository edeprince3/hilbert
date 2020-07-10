#ifndef MISC_H
#define MISC_H

namespace psi{

/// transform three-index integrals to MO basis
void ThreeIndexIntegrals(std::shared_ptr<Wavefunction> ref, long int &nQ, long int memory);

/// after orbital optimization, update Ca/Cb matrices and repack energy-order transformation matrix as pitzer order
void UpdateTransformationMatrix(std::shared_ptr<Wavefunction> ref, std::shared_ptr<Matrix> newMO,
        std::shared_ptr<Matrix> Ca, std::shared_ptr<Matrix> Cb, double * orbopt_transformation_matrix);

}

#endif
