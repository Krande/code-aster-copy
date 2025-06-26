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
!
subroutine gsyste(matr, nchthe, nnoff, gthi, gi)
    implicit none
!
! ......................................................................
!     - FONCTION REALISEE:   FORMATION DU SYSTEME A RESOUDRE
!                            TA A <GI> = TA<G,THETAI>
!
! ENTREE
!
!     MATR         --> VALEURS DE LA MATRICE A
!     NCHTHE       --> NOMBRE DE CHAMPS THETAI
!     GTHI         --> VALEURS DE <G,THETAI>
!     NNOFF        --> NOMBRE DE NOEUDS DU FOND DE FISSURE
!
!  SORTIE
!
!     GI           --> VALEUR DE GI
! ......................................................................
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mgauss.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: istok, nchthe, nnoff
    integer(kind=8) :: i, j, k, iret
!
    real(kind=8) :: gthi(nnoff), gi(nchthe), det
!
    character(len=24) :: matr
!
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iadra1, kk
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(matr, 'L', istok)
    call wkvect('&&GSYSTE.A1', 'V V R8', nchthe*nchthe, iadra1)
!
! INITIALISATION DES VECTEURS ET MATRICES
!
    do i = 1, nchthe
        gi(i) = 0.d0
    end do
!
! CALCUL DU PRODUIT TA*A
!
    do i = 1, nchthe
        do j = 1, nchthe
            do k = 1, nnoff
                kk = iadra1+(i-1)*nchthe+j-1
                zr(kk) = zr(kk)+zr(istok+(k-1)*nchthe+i-1)*zr(istok+(k-1)*nchthe+j-1)
!
            end do
        end do
    end do
!
!  SECOND MEMBRE TAIJ <G,THETHAI>
!
    do i = 1, nchthe
        do j = 1, nnoff
            gi(i) = gi(i)+zr(istok+(j-1)*nchthe+i-1)*gthi(j)
        end do
    end do
!
! RESOLUTION DU SYSTEME LINEAIRE NON SYMETRIQUE PAR GAUSS
!
    call mgauss('NFVP', zr(iadra1), gi, nchthe, nchthe, &
                1, det, iret)
!
    call jedetr('&&GSYSTE.A1')
    call jedema()
end subroutine
