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

subroutine te0571(option, nomte)
! person_in_charge: sam.cuvilliez at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/xrigth.h"
#include "asterfort/xthddl.h"
#include "asterfort/xthini.h"
    character(len=16) :: option, nomte
!
!-----------------------------------------------------------------------
!
!     BUT: CALCUL DES MATRICES DE RIGIDITE ELEMENTAIRES EN THERMIQUE
!          LINEAIRE, ELEMENTS X-FEM LINEAIRES
!
!          OPTION : 'RIGI_THER'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nfh, nfe, igeom, nnop, jpintt, imate, itps, jstno
    integer(kind=8) :: imattt, jcnset, jheavt, jlonch, jbaslo, jlsn, jlst, nddlno
    integer(kind=8) :: heavn(27, 5), ino, ig, jheavn, ncompn, jtab(7), iret
    character(len=8) :: elrefp
!
! ----------------------------------------------------------------------
! --- PREALABLES AU CALCUL DE LA RIGIDITE ELEMENTAIRE
! ----------------------------------------------------------------------
!
!     ON INTERDIT LES ELTS QUADRATIQUES
    call elref1(elrefp)
    ASSERT(iselli(elrefp))
!
!     CHAMPS IN 'CLASSIQUES'
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTR', 'L', itps)
!     CHAMPS IN X-FEM
    call jevech('PSTANO', 'L', jstno)
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PBASLOR', 'L', jbaslo)
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
!     CHAMP OUT
    call jevech('PMATTTR', 'E', imattt)
!
!     ELT DE REF PARENT : RECUP NDIM ET NNOP (NOEUDS PARENT)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop)
!
!     NBRE DE DDLS PAR NOEUD
    call xthini(nomte, nfh, nfe)
    nddlno = 1+nfh+nfe
!
!   RECUPERATION DE LA DEFINITION DES FONCTIONS HEAVISIDES
    if (nfh .gt. 0) then
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
! --- CALCUL DE LA MATRICE DE RIGIDITE ELEMENTAIRE
! ----------------------------------------------------------------------
!
    call xrigth(ndim, elrefp, nnop, imate, itps, &
                igeom, zi(jlonch), zi(jcnset), jpintt, zr(jlsn), &
                zr(jlst), heavn, zr(jbaslo), zi(jheavt), nfh, nfe, &
                zr(imattt))
!
! ----------------------------------------------------------------------
! --- SUPPRESSION DES DDLS SUPERFLUS
! ----------------------------------------------------------------------
!
!     SUPPRESSION DES DDLS SUPERFLUS
    call xthddl(nfh, nddlno, nnop, zi(jstno), option, &
                nomte, mat=zr(imattt))
!
end subroutine
