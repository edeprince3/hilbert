#ifndef MISC_H
#define MISC_H

namespace psi{

/// transform three-index integrals to MO basis
void ThreeIndexIntegrals(std::shared_ptr<Wavefunction> ref, long int &nQ, long int memory);

}

#endif
