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

subroutine pmfdgedef(typfib, b, gg, depl, alicom, nbfibr, nbcarm, &
                     vf, nbassepou, maxfipoutre, nbfipoutre, yj, zj, deffibasse, vfv, deffib)
!
!
! -------------------------------------------------------------------------------------------------
!
!         DEFORMATION GENERALISEE PMF ET DEFORMATION AXIALE DES FIBRES
!
! -------------------------------------------------------------------------------------------------
!
!   IN
!       typfib        : type des fibres : 1 ou 2, ou 3 (squelette assemblage)
!       b             : Matrice b pour la position consideree
!       g             : Matrice g pour la position consideree (mode incompatible)
!       depl          : Champ de deplacement sur l'element
!       alicom        : Variable mode incompatible
!       nbfibr        : nombre de fibres
!       nbcarm        : nombre de caractéristiques sur chaque fibre
!       vf(*)         : positions des fibres
!           Types 1 et 2
!               vf(1,*) : Y fibres
!               vf(2,*) : Z fibres
!               vf(3,*) : Aire fibres
!           Types 2
!               vf(4,*) : Yp groupes de fibres
!               vf(5,*) : Zp groupes de fibres
!               vf(6,*) : num du groupe
!       nbassepou     : nombre de sous-poutres si multipoutre
!       maxfipoutre   : nombre maximum de fibres dans les sous-poutres
!       nbfipoutre(*) : nombre de fibres dans les sous-poutres
!       yj(*)         : position Y des sous-poutres
!       zj(*)         : position Z des sous-poutres
!       deffibasse(*) : deformation des fibres sur une sous-poutre
!       vfv(*)        : tableau vide afin de decrire les proprietes sur une sous-poutre
!       dege(6)       : deformations generalisees sur l'element
!
!   OUT
!       defffib(nf) : deformations axiales des fibres
!
! person_in_charge: jean-luc.flejou at edf.fr
! --------------------------------------------------------------------------------------------------
!

    implicit none
#include "asterfort/utmess.h"
#include "asterfort/pmfdge.h"
#include "asterfort/pmfdef.h"
#include "asterfort/r8inir.h"
!

    integer(kind=8) :: typfib, nbfibr, nbcarm, nbassepou, nbfipoutre(*), maxfipoutre
    real(kind=8) :: vf(nbcarm, nbfibr), dege(6), b(4), gg, depl(*)
    real(kind=8) :: alicom, deffib(nbfibr)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ii, i, pos, posfib
    real(kind=8) :: dege2(6), deffibasse(*), vfv(7, *)
    real(kind=8) :: yj(*), zj(*), depl2(12)
!
! --------------------------------------------------------------------------------------------------
!
    if (typfib .eq. 1) then
!       3 caractéristiques utiles par fibre : y z aire
        call pmfdge(b, gg, depl, alicom, dege)
        call pmfdef(typfib, nbfibr, nbcarm, vf, dege, deffib)

    elseif ((typfib .eq. 2) .or. (typfib .eq. 3)) then
        depl2(:) = 0
        dege2(:) = 0
        call r8inir(nbassepou, 0.d0, yj, 1)
        call r8inir(nbassepou, 0.d0, zj, 1)
        pos = 1
        posfib = 0
        call r8inir(maxfipoutre, 0.d0, deffibasse, 1)
        do i = 1, nbassepou
            call r8inir(maxfipoutre*7, 0.d0, vfv, 1)
            yj(i) = vf(4, pos)
            zj(i) = vf(5, pos)
            !Construction champ des DDLs sur poutres
            if (typfib .eq. 2) then
                depl2(1) = depl(1)+depl(5)*zj(i)-depl(6)*yj(i)
                depl2(2) = depl(2)-depl(4)*zj(i)
                depl2(3) = depl(3)+depl(4)*yj(i)
                depl2(4) = depl(4)
                depl2(5) = depl(5)
                depl2(6) = depl(6)
                depl2(7) = depl(7)+depl(11)*zj(i)-depl(12)*yj(i)
                depl2(8) = depl(8)-depl(10)*zj(i)
                depl2(9) = depl(9)+depl(10)*yj(i)
                depl2(10) = depl(10)
                depl2(11) = depl(11)
                depl2(12) = depl(12)
            else
                depl2(1) = depl(1)+depl(8)*zj(i)-depl(9)*yj(i)
                depl2(2) = depl(2)-depl(7)*zj(i)
                depl2(3) = depl(3)+depl(7)*yj(i)
                depl2(4) = depl(4)
                depl2(5) = depl(5)
                depl2(6) = depl(6)
                depl2(7) = depl(10)+depl(17)*zj(i)-depl(18)*yj(i)
                depl2(8) = depl(11)-depl(16)*zj(i)
                depl2(9) = depl(12)+depl(16)*yj(i)
                depl2(10) = depl(13)
                depl2(11) = depl(14)
                depl2(12) = depl(15)
            end if
            do ii = 1, nbfipoutre(i)
                posfib = pos+ii-1
                vfv(1, ii) = vf(1, posfib)-yj(i)
                vfv(2, ii) = vf(2, posfib)-zj(i)
                vfv(3, ii) = vf(3, posfib)
            end do
            call pmfdge(b, gg, depl2, alicom, dege2)
            call pmfdef(1, maxfipoutre, nbcarm, vfv, dege2, deffibasse)
            do ii = 1, nbfipoutre(i)
                deffib(pos+ii-1) = deffibasse(ii)
            end do
            pos = pos+nbfipoutre(i)
        end do

    else
        call utmess('F', 'ELEMENTS2_40', si=typfib)
    end if

end subroutine
