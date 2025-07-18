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
subroutine vp2trd(type, nbvect, alpha, beta, signes, &
                  vecpro, mxiter, nitqr)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/vp2tru.h"
#include "asterfort/vpordo.h"
#include "asterfort/vpqlts.h"
#include "asterfort/wkvect.h"
    character(len=1) :: type
    integer(kind=8) :: nbvect, mxiter, nitqr
    real(kind=8) :: alpha(*), beta(*), signes(*), vecpro(*)
!     ------------------------------------------------------------------
!     RESOLUTION DU SYSTEME TRIDIAGONAL SYMETRIQUE OU NON SYMETRIQUE.
!     ISSU DE LA METHODE DE LANCZOS.
!     ------------------------------------------------------------------
! IN  TYPE   : K1 : TYPE DU PROBLEME
!       'G' -  LA TRIDIAGONALE EST ISSUE D'UN PROBLEME GENERALISE
!       'Q' -  LA TRIDIAGONALE EST ISSUE D'UN PROBLEME QUADRATIQUE
! IN  NBVECT : I : NOMBRE D'EQUATION (==> DE VECTEURS PROPRES)
! VAR ALPHA  : R :
!        EN ENTREE : ALPHA(I) CONTIENT LE I-EME TERME DIAGONAL
!        EN SORTIE : 'G' - CONTIENT LA PULSATION DU PROBLEME INVERSE
!                    'Q' - CONTIENT IM(VAL_PROP_TRI_DIAG)
! VAR BETA   : R  :
!        EN ENTREE : BETA(I) CONTIENT LE TERME SUR-DIAGONAL A(I-1,I)
!                    PAR CONVENTION BETA(1) = 0
!        EN SORTIE : 'G' - CONTIENT L'AMORTISSEMENT DU PROBLEME INVERSE
!                    'Q' - CONTIENT RE(VAL_PROPTRI_DIAG)
! IN  SIGNES
! IN  MXITER : I : NOMBRE D'ITERATION MAXIMUM POUR LA METHODE (QL/QR)
!            : REMARQUE MXITER = 30 EST UN BON CHOIX.
! OUT NITQR : NOMBRE MAXIMAL D'ITERATIONS ATTEINT AVEC LA METHODE QR
!     ------------------------------------------------------------------
!     REMARQUE : 'G' - TRI SUIVANT LES PULSATION CROISSANTES
!                'Q' - PAS DE TRI
!     ------------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
    real(kind=8) :: symet
    character(len=8) :: method
!     ------------------------------------------------------------------
!
!     ---  ON DETECTE LES FREQUENCES INFERIEURES AU SHIFT ---
!-----------------------------------------------------------------------
    integer(kind=8) :: ier, ivec, ivect, ladw1, ladw2, ladwk1, ladwk2
    integer(kind=8) :: ladz1, ladz2, n2
!-----------------------------------------------------------------------
    call jemarq()
    if (type .eq. 'G') then
        symet = signes(1)
        do ivec = 2, nbvect
            symet = min(signes(ivec), symet)
        end do
    else
        symet = -1.d0
    end if
!
    if (symet .gt. 0.d0) then
!
!        --- CAS OU LA TRIDIAGONALE EST SYMETRIQUE ---
        call vpqlts(alpha, beta, nbvect, vecpro, nbvect, &
                    mxiter, ier, nitqr)
        do ivect = 1, nbvect
            beta(ivect) = 0.0d0
        end do
!
    else
!
!        --- CAS OU LA TRIDIAGONALE N'EST PAS SYMETRIQUE ---
        method = 'TRI_DIAG'
        n2 = 2*nbvect
        if (type .eq. 'G') then
!
            call wkvect('&&VP2TRD.ZONE.TRAV', 'V V R', n2, ladw1)
            call wkvect('&&VP2TRD.WK.VPHQRP', 'V V R', n2, ladwk1)
            call wkvect('&&VP2TRD.Z.VPHQRP ', 'V V R', n2*n2, ladz1)
            call vp2tru(method, type, alpha, beta, signes, &
                        vecpro, nbvect, zr(ladw1), zr(ladz1), zr(ladwk1), &
                        mxiter, ier, nitqr)
            call jedetr('&&VP2TRD.ZONE.TRAV')
            call jedetr('&&VP2TRD.WK.VPHQRP')
            call jedetr('&&VP2TRD.Z.VPHQRP ')
        else
            call wkvect('&&VP2TRD.W.VPHQRP ', 'V V R', n2, ladw2)
            call wkvect('&&VP2TRD.A.VPHQRP ', 'V V R', nbvect*nbvect, ladz2)
            call wkvect('&&VP2TRD.WK.VPHQRP', 'V V R', n2, ladwk2)
!
            call vp2tru(method, type, alpha, beta, signes, &
                        zr(ladz2), nbvect, zr(ladw2), vecpro, zr(ladwk2), &
                        mxiter, ier, nitqr)
            call jedetr('&&VP2TRD.W.VPHQRP  ')
            call jedetr('&&VP2TRD.A.VPHQRP  ')
            call jedetr('&&VP2TRD.WK.VPHQRP ')
        end if
    end if
!
    if (nbvect .eq. 1) vecpro(1) = 1.d0
    if (ier .ne. 0) then
        call utmess('F', 'ALGELINE3_55')
    end if
!
!     --- PASSAGE AUX VALEURS PROPRES DU SYSTEME INITIAL ---
    if (type .eq. 'G') then
        do ivect = 1, nbvect
            if (alpha(ivect) .eq. 0.0d0) then
                call utmess('A', 'ALGELINE3_56')
                alpha(ivect) = 1.d+70
            else
                alpha(ivect) = 1.d0/alpha(ivect)
            end if
        end do
!        --- TRI DES ELEMENTS PROPRES PAR ORDRE CROISSANT DES VALEURS
!            ABSOLUES DES VALEURS PROPRES
        call vpordo(1, 0, nbvect, alpha, vecpro, &
                    nbvect)
    else if (type .eq. 'Q') then
!        -- POUR LE PB Q LA RESTORATION DEPEND DE L'APPROCHE (CF:WP2VEC)
    end if
!
    call jedema()
end subroutine
