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
subroutine fstat0(nbpt, fn, offset, fnmoyt, fnmoyc, &
                  fnrmst, fnrmsc, fnmax, fnmin, fmaxmo, &
                  fminmo, nbmaxr, nbminr)
! CETTE ROUTINE EST EN FAIT L'ANCIENNE FSTAT RENOMMEE FSTAT0
!
!       MOYENNAGE STATISTIQUE DES FORCES AMV
!
!
!
    implicit none
    real(kind=8) :: fn(*), fnmoyt, fnmoyc, fnrmsc, fnrmst, fnmax, fnmin
    real(kind=8) :: sminr, smaxr, offset, fmaxmo, fminmo, sfn2, sfn
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
!            FNETYP       ECART TYPE ( VALEUR DITE PARFOIS RMS )
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
!       NBCOUNT NOMBRE DE VALEURS DU TABLEAU > SEUIL
!       NBMAXR  NOMBRE DE MAXIMAS RELATIFS RENCONTRES
!       NBMINR  NOMBRE DE MINIMAS RELATIFS RENCONTRES
!       SMAXR   SOMME DES MAXIMAS RELATIFS
!       SMINR   SOMME DES MINIMAS RELATIFS
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, nbmaxr, nbminr, nbpt, ncount
!-----------------------------------------------------------------------
    sfn = 0.d0
    sfn2 = 0.d0
    fnmax = -10.d20
    fnmin = -fnmax
    ncount = 0
    nbminr = 0
    nbmaxr = 0
    smaxr = 0.0d0
    sminr = 0.0d0
!
    do i = 1, nbpt
!
        if ((abs(fn(i))) .gt. offset) then
            ncount = ncount+1
            sfn = sfn+fn(i)
!
!           RECHERCHE DES EXTREMAS ABSOLUS
!
            if (fn(i) .gt. fnmax) fnmax = fn(i)
            if (fn(i) .lt. fnmin) fnmin = fn(i)
!
        end if
!
    end do
!
    do i = 2, nbpt-1
!
        if ((abs(fn(i))) .gt. offset) then
!
!           RECHERCHE DES EXTREMAS RELATIFS
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
!
        end if
!
    end do
!
!
    if (ncount .ne. 0) then
        fnmoyc = sfn/dble(ncount)
        fnmoyt = sfn/dble(nbpt)
!
    else
        fnmoyc = 0.d0
        fnmoyt = sfn/dble(nbpt)
    end if
!
    if (nbminr .ne. 0) then
        fminmo = sminr/dble(nbminr)
!
    else
        fminmo = 0.d0
    end if
!
    if (nbmaxr .ne. 0) then
        fmaxmo = smaxr/dble(nbmaxr)
!
    else
        fmaxmo = 0.d0
    end if
!
    do i = 1, nbpt
        if (abs(fn(i)) .gt. offset) then
            sfn2 = sfn2+fn(i)**2
        end if
!
    end do
!
    if (ncount .ne. 0) then
        fnrmsc = sqrt(sfn2/dble(ncount))
        fnrmst = sqrt(sfn2/dble(nbpt))
!
    else
        fnrmsc = 0.d0
        fnrmst = sqrt(sfn2/dble(nbpt))
!
        fnmin = 0.d0
        fnmax = 0.d0
    end if
!
    goto 30
!
30  continue
end subroutine
