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

subroutine te0594(option, nomte)
! person_in_charge: sam.cuvilliez at edf.fr
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/xrechp.h"
#include "asterfort/xthddl.h"
#include "asterfort/xthini.h"
    character(len=16) :: option, nomte
!
!-----------------------------------------------------------------------
!
!     BUT: THERMIQUE LINEAIRE / ELEMENTS PRINCIPAUX X-FEM LINEAIRES
!          ECHANGE_PAROI POUR FISSURES X-FEM
!
!          OPTION : 'RIGI_THER_PARO_F' ET 'RIGI_THER_PARO_R'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nfh, nfe, igeom, nnop, jptint, jcface
    integer(kind=8) :: jlonch, jlst, itps, ihechp, jstno, jbasec
    integer(kind=8) :: heavn(27, 5), ino, ig, jheavn, ncompn, jtab(7), iret
    integer(kind=8) :: imattt, nddlno
    character(len=8) :: elrefp
    character(len=4) :: fonree
!
! ----------------------------------------------------------------------
! --- PREALABLES AU CALCUL
! ----------------------------------------------------------------------
!
!     ON INTERDIT LES ELTS QUADRATIQUES
    call elref1(elrefp)
    ASSERT(iselli(elrefp))
!
!     CHAMPS IN CLASSIQUES
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itps)
!
!     SI LE COEFF D'ECHANGE PAROI EST NUL, IL N'Y A RIEN A FAIRE
    if (option .eq. 'RIGI_THER_PARO_R') then
        fonree = 'REEL'
        call jevech('PHECHPR', 'L', ihechp)
        if (abs(zr(ihechp)) .lt. r8prem()) goto 999
    else if (option .eq. 'RIGI_THER_PARO_F') then
        fonree = 'FONC'
        call jevech('PHECHPF', 'L', ihechp)
        if (zk8(ihechp) .eq. '&FOZERO ') goto 999
    else
        ASSERT(.false.)
    end if
!
!     CHAMPS IN X-FEM
    call jevech('PPINTER', 'L', jptint)
    call jevech('PCFACE', 'L', jcface)
    call jevech('PLONGCO', 'L', jlonch)
    call jevech('PLST', 'L', jlst)
    call jevech('PSTANO', 'L', jstno)
    call jevech('PBASECO', 'L', jbasec)
!
!     CHAMPS OUT
    call jevech('PMATTTR', 'E', imattt)
!
!     ELT DE REF PARENT : RECUP NDIM ET NNOP (NOEUDS PARENT)
!     -> RQ : 'RIGI' POUR LA FAMILLE DE PG EST DONC SANS CONSQUENCE
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop)
!
!     RECUP DE NFH (NBRE FCT HEAVISIDE) ET NFE (NBRE FCT SINGULIER)
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
! --- CALCUL DE LA MATRICE DE RIGIDITE ELEMENTAIRE DUE A ECHANGE_PAROI
! ----------------------------------------------------------------------
!
    call xrechp(ndim, elrefp, nnop, igeom, itps, &
                ihechp, jptint, jcface, jlonch, &
                jlst, jbasec, nfh, nfe, fonree, &
                imattt, heavn)
!
! ----------------------------------------------------------------------
! --- SUPPRESSION DES DDLS SUPERFLUS
! ----------------------------------------------------------------------
!
    call xthddl(nfh, nddlno, nnop, zi(jstno), option, &
                nomte, mat=zr(imattt))
!
999 continue
!
end subroutine
