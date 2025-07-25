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
subroutine dbgcal(optioz, ifm, nbin, lpaiz, lchiz, &
                  nbout, lpaouz, lchouz)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterfort/jelstc.h"
#include "asterfort/utimsd.h"
#include "asterfort/utmess.h"
    character(len=*) :: optioz
    integer(kind=8) :: ifm
    integer(kind=8) :: nbin, nbout
    character(len=*) :: lpaiz(nbin), lpaouz(nbout)
    character(len=*) :: lchiz(nbin), lchouz(nbout)
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE POUR CALCUL
!
! DEBUGAGE DES CHAMPS IN/OUT POUR CALCUL
!
! ----------------------------------------------------------------------
!
!
! IN  OPTION : OPTION CALCULEE
! IN  IFM    : UNITE LOGIQUE D'IMPRESSION
! IN  NBIN   : NOMBRE DE CHAMPS IN
! IN  NBOUT  : NOMBRE DE CHAMPS OUT
! IN  LPAIZ  : NOM DES TYPES DE CHAMP D'ENTREE
! IN  LCHIZ  : NOM DES CHAMPS D'ENTREE
! IN  LPAOUZ : NOM DES TYPES DE CHAMP DE SORTIE
! IN  LCHOUZ : NOM DES CHAMPS DE SORTIE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ich
    character(len=16) :: option
!
    character(len=8) :: k8bid
    integer(kind=8) :: nbval, nbobj
!
! ---------------------------------------------------------------------
!
    option = optioz
!
    write (ifm, *) '***** CALCUL DE L OPTION <', option, '>'
    write (ifm, *) ' ** NBRE CHAMPS IN : ', nbin
    write (ifm, *) ' ** NBRE CHAMPS OUT: ', nbout
!
    write (ifm, *) '***** <CHAMPS_IN>'
    do ich = 1, nbin
        write (ifm, *) ' * CHAMP IN  <', ich, '>'
        write (ifm, *) ' * PARAMETRE <', lpaiz(ich), '>'
        write (ifm, *) ' * CHAMP     <', lchiz(ich), '>'
        if (lpaiz(ich) (1:1) .eq. ' ') then
            call utmess('A', 'PRECALCUL_60', si=ich)
        end if
        if (lchiz(ich) (1:1) .eq. ' ') then
            call utmess('A', 'PRECALCUL_61', si=ich)
        end if
!
        call jelstc(' ', lchiz(ich) (1:19), 1, 0, k8bid, &
                    nbval)
        nbobj = -nbval
        if (nbobj .eq. 0) then
            call jelstc(' ', lchiz(ich), 1, 0, k8bid, &
                        nbval)
            nbobj = -nbval
            if (nbobj .eq. 0) then
                write (ifm, *) ' * SD INTROUVABLE !'
            else
                write (ifm, *) ' * RESUME DE LA SD :'
                call utimsd(ifm, -1, .true._1, .true._1, lchiz(ich), &
                            1, ' ', perm='OUI')
            end if
        else
            write (ifm, *) ' * RESUME DE LA SD :'
            call utimsd(ifm, -1, .true._1, .true._1, lchiz(ich) (1:19), &
                        1, ' ', perm='OUI')
        end if
    end do
!
    write (ifm, *) '***** <CHAMPS_OUT>'
    do ich = 1, nbout
        write (ifm, *) ' * CHAMP OUT <', ich, '>'
        write (ifm, *) ' * PARAMETRE <', lpaouz(ich), '>'
        write (ifm, *) ' * CHAMP     <', lchouz(ich), '>'
        if (lpaouz(ich) (1:1) .eq. ' ') then
            call utmess('A', 'PRECALCUL_62', si=ich)
        end if
        if (lchouz(ich) (1:1) .eq. ' ') then
            call utmess('A', 'PRECALCUL_63', si=ich)
        end if
!
        call jelstc(' ', lchouz(ich) (1:19), 1, 0, k8bid, &
                    nbval)
        nbobj = -nbval
        if (nbobj .eq. 0) then
            call jelstc(' ', lchouz(ich), 1, 0, k8bid, &
                        nbval)
            nbobj = -nbval
            if (nbobj .eq. 0) then
                write (ifm, *) ' * SD INTROUVABLE !'
            else
                write (ifm, *) ' * RESUME DE LA SD :'
                call utimsd(ifm, -1, .true._1, .true._1, lchouz(ich), &
                            1, ' ', perm='OUI')
            end if
        else
            write (ifm, *) ' * RESUME DE LA SD :'
            call utimsd(ifm, -1, .true._1, .true._1, lchouz(ich) (1:19), &
                        1, ' ', perm='OUI')
        end if
    end do
end subroutine
