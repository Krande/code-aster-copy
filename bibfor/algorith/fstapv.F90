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
subroutine fstapv(nbpt, fn, t, offset, fnmoyt, &
                  fnmoyc, fnrmst, fnrmsc, fnmax, fnmin, &
                  fmaxmo, fminmo, sfn, sfn2, tchoc, &
                  nbmaxr, nbminr)
!
!       MOYENNAGE STATISTIQUE DES FORCES
!       ALGORITHME TEMPOREL A PAS VARIABLE
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/trapez.h"
#include "asterfort/wkvect.h"
    real(kind=8) :: fn(*), t(*), fnmoyt, fnmoyc, fnrmsc, fnrmst, fnmax, fnmin
    real(kind=8) :: sminr, smaxr, offset, fmaxmo, fminmo, sfn2, sfn, tchoc
!
!
!       ARGUMENTS:
!       ----------------------------------------
!       IN:
!            NBPT         NB DE POINTS DU TABLEAU A ANALYSER
!            FN           TABLEAU A ANALYSER
!            OFFSET       VALEUR DE SEUIL DE DETECTION DES VALEURS
!
!       OUT:
!            FNMOY        VALEUR MOYENNE ( COMPTAGE AU DESSUS DU SEUIL )
!            FNRMS        SQRT DE LA MOYENNE DES CARRES ( RMS POUR DES F
!                         REDRESSEES )
!            FNMAX        VALEUR MAXIMUM ABSOLU DU TABLEAU
!            FNMIN        VALEUR MINIMUM ABSOLU DE LA FONCTION
!            FMAXMO       MOYENNE DES MAXIMAS RELATIFS DE LA FONCTION
!            FMINMO       MOYENNE DES MINIMAS RELATIFS DE LA FONCTION
!
!
!
!       VARIABLES UTILISEES
!       ----------------------------------------
!       SFN SOMME DES FORCES DE CHOC
!       SFN2 SOMME DES CARRES DES FORCES DE CHOC
!       NBMAXR  NOMBRE DE MAXIMAS RELATIFS RENCONTRES
!       NBMINR  NOMBRE DE MINIMAS RELATIFS RENCONTRES
!       SMAXR   SOMME DES MAXIMAS RELATIFS
!       SMINR   SOMME DES MINIMAS RELATIFS
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ift, itct, nbmaxr, nbminr, nbpt
    real(kind=8) :: ttot
!-----------------------------------------------------------------------
    call jemarq()
    sfn = 0.d0
    sfn2 = 0.d0
    fnmax = -10.d20
    fnmin = -fnmax
    nbminr = 0
    nbmaxr = 0
    smaxr = 0.0d0
    sminr = 0.0d0
!
    call wkvect('&&FSTAPV.TEMP.FCNT', 'V V R', nbpt, ift)
    call wkvect('&&FSTAPV.TEMP.FTCT', 'V V R', nbpt, itct)
!
!           RECHERCHE DES EXTREMAS ABSOLUS
!
    do i = 1, nbpt
        if (fn(i) .gt. fnmax) fnmax = fn(i)
        if (fn(i) .lt. fnmin) fnmin = fn(i)
    end do
!
!
!
! --- RECHERCHE DES EXTREMAS RELATIFS
!
    do i = 2, nbpt-1
!
        if ((fn(i) .gt. fn(i-1)) .and. (fn(i) .gt. fn(i+1))) then
            smaxr = smaxr+fn(i)
            nbmaxr = nbmaxr+1
        end if
!
        if ((fn(i) .lt. fn(i-1)) .and. (fn(i) .lt. fn(i+1))) then
            sminr = sminr+fn(i)
            nbminr = nbminr+1
        end if
    end do
!
!
! --- MOYENNE DES EXTREMAS RELATIFS
!
    if (nbminr .ne. 0) then
        fminmo = sminr/dble(nbminr)
    else
        fminmo = 0.d0
    end if
!
    if (nbmaxr .ne. 0) then
        fmaxmo = smaxr/dble(nbmaxr)
    else
        fmaxmo = 0.d0
    end if
!
!
! --- CALCUL DE LA FORCE MOYENNE
!
    do i = 1, nbpt
        if ((abs(fn(i))) .gt. offset) then
            zr(ift+i-1) = fn(i)
            zr(itct+i-1) = 1.d0
        else
            zr(ift+i-1) = 0.d0
            zr(itct+i-1) = 0.d0
        end if
    end do
    call trapez(t, zr(ift), nbpt, sfn)
    call trapez(t, zr(itct), nbpt, tchoc)
    ttot = t(nbpt)-t(1)
    fnmoyt = sfn/ttot
    if (tchoc .gt. 0.d0) then
        fnmoyc = sfn/tchoc
    else
        fnmoyc = 0.d0
    end if
!
!
! --- CALCUL DE LA FORCE QUADRATIQUE MOYENNE
!
    do i = 1, nbpt
        if (abs(fn(i)) .gt. offset) then
            zr(ift+i-1) = fn(i)*fn(i)
        else
            zr(ift+i-1) = 0.d0
        end if
    end do
    call trapez(t, zr(ift), nbpt, sfn2)
    fnrmst = sqrt(sfn2/ttot)
    if (tchoc .gt. 0.d0) then
        fnrmsc = sqrt(sfn2/tchoc)
    else
        fnrmsc = 0.d0
        fnmin = 0.d0
        fnmax = 0.d0
    end if
!
!
    call jedetr('&&FSTAPV.TEMP.FCNT')
    call jedetr('&&FSTAPV.TEMP.FTCT')
    call jedema()
end subroutine
