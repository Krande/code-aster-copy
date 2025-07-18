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

subroutine gcsele(motcle, chvolu, ch1d2d, ch2d3d, chpres, &
                  chepsi, chpesa, chrota, lvolu, l1d2d, &
                  l2d3d, lpres, lepsi, lpesa, lrota, lfchar, &
                  lfvolu, lf1d2d, lf2d3d, lfpres, lfepsi, &
                  lfpesa, lfrota, carte0, lpchar, &
                  lccomb)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
! aslint: disable=W1504
!
    character(len=16) :: motcle
    character(len=19) :: carte0
    aster_logical :: lpchar, lccomb
    character(len=19) :: chvolu, ch1d2d, ch2d3d, chpres
    character(len=19) :: chepsi, chpesa, chrota
    aster_logical :: lvolu, l1d2d, l2d3d, lpres, lfchar
    aster_logical :: lepsi, lpesa, lrota
    aster_logical :: lfvolu, lf1d2d, lf2d3d, lfpres
    aster_logical :: lfepsi, lfpesa, lfrota
!
! ----------------------------------------------------------------------
!
! ROUTINE CALC_G
!
! SELECTION DES VARIABLES CORRESPONDANT AU MOT-CLEF ACTIF
!
! ----------------------------------------------------------------------
!
!
! IN  MOTCLE : MOT-CLEF DE LA CHARGE (VOIR LISDEF)
! I/O lvolu  : .TRUE.  ON A UNE CHARGE FORCE_INTERNE
! I/O l1d2d  : .TRUE.  SI ON AU MOINS UNE CHARGE FORCE_CONTOUR
! I/O l2d3d  : .TRUE.  SI ON AU MOINS UNE CHARGE FORCE_FACE
! I/O lpres  : .TRUE.  SI ON AU MOINS UNE CHARGE PRES_REP
! I/O lepsi  : .TRUE.  SI ON AU MOINS UNE CHARGE EPSI_INIT
! I/O lpesa  : .TRUE.  SI ON AU MOINS UNE CHARGE PESANTEUR
! I/O lrota  : .TRUE.  SI ON AU MOINS UNE CHARGE ROTATION
! IN  lfvolu : .TRUE.  SI CHARGE DE TYPE 'FONCTION'
! OUT  lfvolu : .TRUE.  SI CHARGE FORCE_INTERNE DE TYPE 'FONCTION'
! OUT  lf1d2d : .TRUE.  SI CHARGE FORCE_CONTOUR DE TYPE 'FONCTION'
! OUT  lf2d3d : .TRUE.  SI CHARGE FORCE_FACE DE TYPE 'FONCTION'
! OUT  lfpres : .TRUE.  SI CHARGE PRES_REP DE TYPE 'FONCTION'
! OUT  lfepsi : .TRUE.  SI CHARGE EPSI_INIT DE TYPE 'FONCTION'
! OUT  lfpesa : .TRUE.  SI CHARGE PESANTEUR DE TYPE 'FONCTION'
! OUT  lfrota : .TRUE.  SI CHARGE ROTATION DE TYPE 'FONCTION'
! IN  CHVOLU : CARTE POUR FORCE_INTERNE
! IN  ch1d2d : CARTE POUR FORCE_CONTOUR
! IN  ch2d3d : CARTE POUR FORCE_FACE
! IN  chpres : CARTE POUR PRES_REP
! IN  chepsi : CARTE POUR EPSI_INIT
! IN  chpesa : CARTE POUR PESANTEUR
! IN  chrota : CARTE POUR ROTATION
! OUT LPCHAR : .TRUE.  SI C'EST LA PREMIERE FOIS QU'ON A UNE CHARGE DU STYLE COURANT
! OUT CARTE0 : CARTE DE LA CHARGE DU STYLE COURANT
! OUT LFORMU : .TRUE.  SI LE CHARGEMENT 'FONCTION' UTILISE UNE FORMULE
! OUT LCCOMB : .TRUE. SI LE CHARGEMENT EST COMBINABLE
!
! ----------------------------------------------------------------------
!
    lpchar = .false.
    if (motcle .eq. 'FORCE_INTERNE#2D' .or. motcle .eq. 'FORCE_INTERNE#3D') then
        carte0 = chvolu
        if (.not. lvolu) lpchar = .true.
        lvolu = .true.
        lfvolu = lfchar
        lccomb = .true.
    else if (motcle .eq. 'FORCE_CONTOUR') then
        carte0 = ch1d2d
        if (.not. l1d2d) lpchar = .true.
        l1d2d = .true.
        lf1d2d = lfchar
        lccomb = .true.
    else if (motcle .eq. 'FORCE_FACE') then
        carte0 = ch2d3d
        if (.not. l2d3d) lpchar = .true.
        l2d3d = .true.
        lf2d3d = lfchar
        lccomb = .true.
    else if (motcle .eq. 'PRES_REP') then
        carte0 = chpres
        if (.not. lpres) lpchar = .true.
        lpres = .true.
        lfpres = lfchar
        lccomb = .true.
    else if (motcle .eq. 'EPSI_INIT') then
        carte0 = chepsi
        if (.not. lepsi) lpchar = .true.
        lepsi = .true.
        lfepsi = lfchar
        lccomb = .false.
    else if (motcle .eq. 'PESANTEUR') then
        carte0 = chpesa
        if (.not. lpesa) lpchar = .true.
        lpesa = .true.
        lfpesa = lfchar
        lccomb = .false.
    else if (motcle .eq. 'ROTATION') then
        carte0 = chrota
        if (.not. lrota) lpchar = .true.
        lrota = .true.
        lfrota = lfchar
        lccomb = .false.
    else
        write (6, *) 'MOT-CLEF:', motcle
        ASSERT(.false.)
    end if
!
end subroutine
