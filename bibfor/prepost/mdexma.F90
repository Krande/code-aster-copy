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

subroutine mdexma(nofimd, idfimd, nomamd, option, existm, &
                  ndim, codret)
!_____________________________________________________________________
! person_in_charge: nicolas.sellenet at edf.fr
!        FORMAT MED : EXISTENCE D'UN MAILLAGE DANS UN FICHIER
!               - -   --             --
! ______________________________________________________________________
! .        .     .        .                                            .
! .  NOM   . E/S . TAILLE .           DESCRIPTION                      .
! .____________________________________________________________________.
! . NOFIMD .  E  .   1    . NOM DU FICHIER MED                         .
! . NOFIMD .  E  .   1    . OU NUMERO DU FICHIER DEJA OUVERT           .
! . NOMAMD .  E  .   1    . NOM DU MAILLAGE MED VOULU                  .
! . OPTION .  E  .   1    . QUE FAIT-ON SI LE MAILLAGE EST ABSENT :    .
! .        .     .        . 0 : RIEN DE SPECIAL                        .
! .        .     .        . 1 : ON IMPRIME LA LISTE DES MAILLAGES      .
! .        .     .        .     PRESENTS                               .
! . EXISTM .  S  .   1    . .TRUE. OU .FALSE., SELON QUE LE MAILLAGE   .
! .        .     .        . EST PRESENT OU NON                         .
! . NDIM   .  S  .   1    . LA DIMENSION DU MAILLAGE QUAND IL EXISTE   .
! . CODRET .  S  .   1    . CODE DE RETOUR DES MODULES                 .
! ______________________________________________________________________
!
!====
! 0. DECLARATIONS ET DIMENSIONNEMENT
!====
!
    use as_med_module, only: as_med_open
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/as_mmhmii.h"
#include "asterfort/as_mmhnmh.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
    character(len=*) :: nofimd, nomamd
!
    aster_logical :: existm, ficexi, dejouv
!
    integer(kind=8) :: option, ndim, codret
!
! 0.2. ==> COMMUNS
!
! 0.3. ==> VARIABLES LOCALES
!
    integer(kind=8) :: edlect
    parameter(edlect=0)
    integer(kind=8) :: ednstr
    parameter(ednstr=0)
!
!
    integer(kind=8) :: lnomam, nbmaie
    med_idt :: idfimd
    integer(kind=8) :: iaux, jaux, kaux, tyaux
    integer(kind=8) :: vali(2)
!
    character(len=8) :: saux08
    character(len=64) :: noma64
    character(len=64) :: saux64
    character(len=200) :: daux
    character(len=24) :: valk
! ______________________________________________________________________
!
!====
! 1. ON OUVRE LE FICHIER EN LECTURE
!    ON PART DU PRINCIPE QUE SI ON N'A PAS PU OUVRIR, C'EST QUE LE
!    FICHIER N'EXISTE PAS, DONC SANS MAILLAGE A FORTIORI
!====
!
    existm = .false.
    codret = 0
    inquire (file=nofimd, exist=ficexi)
!
    if (.not. ficexi) then
!
        existm = .false.
        codret = 0
!
    else
        if (idfimd .eq. 0) then
            call as_med_open(idfimd, nofimd, edlect, iaux)
            dejouv = .false.
        else
            dejouv = .true.
            iaux = 0
        end if
        if (iaux .eq. 0) then
!
!====
! 2. LE MAILLAGE EST-IL PRESENT ?
!====
!
! 2.1. ==> COMBIEN DE MAILLAGES DANS LE FICHIER
!
            call as_mmhnmh(idfimd, nbmaie, codret)
            if (codret .ne. 0) then
                saux08 = 'mmhnmh'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
!
! 2.2. ==> RECHERCHE DU NUMERO ET DE LA DIMENSION DU MAILLAGE VOULU
!
!
            noma64 = '                                  '//'                    '
            lnomam = lxlgut(nomamd)
            noma64(1:lnomam) = nomamd(1:lnomam)
!
            do iaux = 1, nbmaie
                saux64 = ' '
                call as_mmhmii(idfimd, iaux, saux64, kaux, tyaux, &
                               daux, codret)
                if (codret .ne. 0) then
                    saux08 = 'mmhmii'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
                if (tyaux .ne. ednstr) then
                    call utmess('A', 'MED_79')
                end if
                jaux = lxlgut(saux64)
                if (jaux .eq. lnomam) then
                    if (saux64 .eq. noma64) then
                        ndim = kaux
                        existm = .true.
                        goto 221
                    end if
                end if
            end do
!
            existm = .false.
!
! 2.3. ==> IMPRESSION EVENTUELLE DES MAILLAGES PRESENTS
!
            if (option .ne. 0) then
!
                valk = nofimd
                vali(1) = nbmaie
                call utmess('A+', 'MED_88', sk=valk, si=vali(1))
                do iaux = 1, nbmaie
                    saux64 = ' '
                    call as_mmhmii(idfimd, iaux, saux64, kaux, tyaux, &
                                   daux, codret)
                    jaux = lxlgut(saux64)
                    valk = saux64(1:jaux)
                    call utmess('A+', 'MED_85', sk=valk)
                    if (tyaux .ne. ednstr) then
                        call utmess('A', 'MED_79')
                    end if
                end do
                call utmess('A', 'MED_80', sk=noma64(1:lnomam))
!
            end if
!
221         continue
!
! 2.3. ==> FERMETURE DU FICHIER
!
            if (.not. dejouv) then
                call as_mficlo(idfimd, codret)
                if (codret .ne. 0) then
                    saux08 = 'mficlo'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
                idfimd = 0
            end if
!
        end if
    end if
!
end subroutine
