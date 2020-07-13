/* 
 *  @BEGIN LICENSE
 * 
 *  Hilbert: a space for quantum chemistry plugins to Psi4 
 * 
 *  Copyright (c) 2020 by its authors (LICENSE).
 * 
 *  The copyrights for code used from other parties are included in
 *  the corresponding files.
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 * 
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 *  @END LICENSE
 */

#include <psi4/psi4-dec.h>
#include <psi4/psifiles.h>
#include <psi4/libiwl/iwl.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libpsi4util/PsiOutStream.h>

#include "v2rdm_solver.h"

using namespace psi;



namespace psi{namespace v2rdm_doci{

void v2RDMSolver::ReadAllIntegrals(iwlbuf *Buf) {

  unsigned long int lastbuf;
  Label *lblptr;
  Value *valptr;

  unsigned long int idx, p, q, r, s, pq, rs, pqrs;

  lblptr = Buf->labels;
  valptr = Buf->values;
  lastbuf = Buf->lastbuf;

  outfile->Printf("\n");
  outfile->Printf("        Read integrals......");
  /**
    * first buffer (read in when Buf was initialized)
    */
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = (unsigned long int) lblptr[idx++];
      q = (unsigned long int) lblptr[idx++];
      r = (unsigned long int) lblptr[idx++];
      s = (unsigned long int) lblptr[idx++];
      double val = (double)valptr[Buf->idx];

      // none of this will work with frozen virtuals ...

      int hp = symmetry_full[p];
      int hq = symmetry_full[q];
      int hpq = SymmetryPair(hp,hq);

      pq = ibas_really_full_sym[hpq][p][q];
      rs = ibas_really_full_sym[hpq][r][s];

      long int offset = 0;
      for (int h = 0; h < hpq; h++) {
          offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
      }
      tei_full_sym_[offset + INDEX(pq,rs)] = val;

  }
  /**
    * now do the same for the rest of the buffers
    */
  while(!lastbuf){
      iwl_buf_fetch(Buf);
      lastbuf = Buf->lastbuf;
      for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {

          p = (unsigned long int) lblptr[idx++];
          q = (unsigned long int) lblptr[idx++];
          r = (unsigned long int) lblptr[idx++];
          s = (unsigned long int) lblptr[idx++];
          double val = (double)valptr[Buf->idx];

          int hp = symmetry_full[p];
          int hq = symmetry_full[q];
          int hpq = SymmetryPair(hp,hq);

          pq = ibas_really_full_sym[hpq][p][q];
          rs = ibas_really_full_sym[hpq][r][s];

          long int offset = 0;
          for (int h = 0; h < hpq; h++) {
              offset += (long int)gems_full[h] * ( (long int)gems_full[h] + 1 ) / 2;
          }
          tei_full_sym_[offset + INDEX(pq,rs)] = val;

      }

  }
  outfile->Printf("done.\n\n");
}


void v2RDMSolver::GetTEIFromDisk(){
  struct iwlbuf Buf;
  iwl_buf_init(&Buf,PSIF_MO_TEI,0.0,1,1);
  ReadAllIntegrals(&Buf);
  iwl_buf_close(&Buf,1);
}


}}
