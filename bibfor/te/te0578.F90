! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine te0578(option, nomte)
! person_in_charge: sam.cuvilliez at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/xtelga.h"
#include "asterfort/xthini.h"
    character(len=16) :: option, nomte
!
!-----------------------------------------------------------------------
!
!     BUT: THERMIQUE LINEAIRE / ELEMENTS X-FEM LINEAIRES
!          OPTION : 'TEMP_ELGA'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nfh, nfe, itempn, igeom, nnop, itempg, jpintt
    integer(kind=8) :: jcnset, jheavt, jlonch, jbaslo, jlsn, jlst
    integer(kind=8) :: heavn(27, 5), ino, ig, jheavn, ncompn, jtab(7), iret
    character(len=8) :: elrefp
!
! ----------------------------------------------------------------------
! --- PREALABLES AU CALCUL
! ----------------------------------------------------------------------
!
    option = option
!     ON INTERDIT LES ELTS QUADRATIQUES
    call elref1(elrefp)
    ASSERT(iselli(elrefp))
!
!     CHAMPS IN
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PTEMPER', 'L', itempn)
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PBASLOR', 'L', jbaslo)
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
!     CHAMPS OUT
    call jevech('PTEMP_R', 'E', itempg)
!
!     ELT DE REF PARENT : RECUP NDIM ET NNOP (NOEUDS PARENT)
!     -> RQ : 'RIGI' POUR LA FAMILLE DE PG EST DONC SANS CONSQUENCE
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop)
!
!     RECUP DE NFH (NBRE FCT HEAVISIDE) ET NFE (NBRE FCT SINGULIER)
    call xthini(nomte, nfh, nfe)
!
!   RECUPERATION DE LA DEFINITION DES FONCTIONS HEAVISIDES
    if (nfh .gt. 0 .or. nfe .gt. 0) then
        call jevech('PHEA_NO', 'L', jheavn)
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                    itab=jtab)
        ncompn = jtab(2)/jtab(3)
        ASSERT(ncompn .eq. 5)
        do ino = 1, nnop
            do ig = 1, ncompn
                heavn(ino, ig) = zi(jheavn-1+ncompn*(ino-1)+ig)
            end do
        end do
    end if
!
! ----------------------------------------------------------------------
! --- CALCUL DU CHAMP DE TEMPERATURE AUX POINTS DE GAUSS
! ----------------------------------------------------------------------
!
    call xtelga(ndim, elrefp, nnop, igeom, zr(itempn), &
                zi(jlonch), zi(jcnset), jpintt, zr(jlsn), zr(jlst), &
                heavn, zr(jbaslo), zi(jheavt), nfh, nfe, zr(itempg))
!
end subroutine
