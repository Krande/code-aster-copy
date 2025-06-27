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
subroutine rc32env(lieu, futotenv, fuseis_in)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rc32env2.h"
#include "asterfort/wkvect.h"
!
    character(len=4) :: lieu
    real(kind=8) :: futotenv, fuseis_in
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE B3200 et ZE200
!     AFFICHAGE DES RESULTATS DANS LA TABLE DE SORTIE
!
!
!  ------------------------------------------------------------------
    integer(kind=8) :: ndim, jfact, num1, num2, k, noccpris, ind1, ind2, jmax
    integer(kind=8) :: i, jfactenv, nb
    real(kind=8) :: fuseism, fuunit, ke, fen, fenel
!
! DEB ------------------------------------------------------------------
!
    call jemarq()
!
!     ----------------------------------------------------------------
! --- CREATION D'UN VECTEUR QUI CONTIENT LES DONNEES ENVIRONNEMENTALES
!     ----------------------------------------------------------------
    ndim = 3*200
    call wkvect('&&RC3200.FACTENV.'//lieu, 'V V R', ndim, jfactenv)
    do i = 1, ndim
        zr(jfactenv+i-1) = 0.d0
    end do
!
!     ----------------------------------------------------------------
! --- BALAYAGE DES COMBINAISONS DE SITUATIONS QUI INTERVIENNENT DANS FU_TOTAL
!     ----------------------------------------------------------------
!
    call jeveuo('&&RC3200.FACT.'//lieu, 'L', jfact)
    call jeveuo('&&RC3200.MAX_RESU.'//lieu, 'L', jmax)
    call getfac('SITUATION', nb)
!
    k = 0
    futotenv = 0.d0
!
555 continue
!
    num1 = zi(jfact+6*k)
    num2 = zi(jfact+6*k+1)
!
    if (num1 .eq. 0) goto 666
    noccpris = zi(jfact+6*k+4)
!
    if (zi(jfact+6*k+5) .eq. 2) then
!-------- le séisme intervient dans cette combinaison
        call jeveuo('&&RC3200.SITUS_RESU.'//lieu, 'L', ind1)
        call jeveuo('&&RC3200.COMBS_RESU.'//lieu, 'L', ind2)
!
        fuseism = fuseis_in
    end if
    if (zi(jfact+6*k+5) .eq. 1) then
!-------- le séisme n'intervient pas dans cette combinaison
        call jeveuo('&&RC3200.SITU_RESU.'//lieu, 'L', ind1)
        call jeveuo('&&RC3200.COMB_RESU.'//lieu, 'L', ind2)
        fuseism = 0.d0
    end if
!
!---- une situation seule a le plus grand fu unitaire
    if (num1 .eq. num2) then
        fuunit = zr(ind1+123*(num1-1)+15)+fuseism
        ke = zr(ind1+123*(num1-1)+18)
        call rc32env2(num1, num2, ke, lieu, fen)
        call rc32env2(num1, num2, 1.0, lieu, fenel)
        zr(jfactenv+3*k) = fen
        zr(jfactenv+3*k+1) = fenel
        zr(jfactenv+3*k+2) = fen*fuunit*noccpris
        futotenv = futotenv+fen*fuunit*noccpris
!
!---- une combinaison de situations a le plus grand fu unitaire
    else
        fuunit = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+17)+fuseism
        ke = zr(ind2+25*nb*(num1-1)+25*(num2-1)-1+20)
        call rc32env2(num1, num2, ke, lieu, fen)
        call rc32env2(num1, num2, 1.0, lieu, fenel)
        zr(jfactenv+3*k) = fen
        zr(jfactenv+3*k+1) = fenel
        zr(jfactenv+3*k+2) = fen*fuunit*noccpris
        futotenv = futotenv+fen*fuunit*noccpris
    end if
!
    k = k+1
    goto 555
!
666 continue
!
    call jedema()
!
end subroutine
