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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmchai(tychap, tyvarz, vali, tychap_out)
    !
    implicit none
    !
#include "asterc/indik8.h"
#include "asterfort/assert.h"
    !
    character(len=6), intent(in) :: tychap
    character(len=*), intent(in) :: tyvarz
    integer(kind=8) :: vali
    character(len=6), optional, intent(out) :: tychap_out
    !
    ! ----------------------------------------------------------------------
    !
    ! ROUTINE MECA_NON_LINE (CALCUL - UTILITAIRE)
    !
    ! INDEX OU EST STOCKE LE NOM DE LA VARIABLE DANS UNE VARIABLE CHAPEAU
    !
    ! ----------------------------------------------------------------------
    !
    !
    ! IN  TYCHAP : TYPE DE VARIABLE CHAPEAU
    !                MEELEM - NOMS DES MATR_ELEM
    !                MEASSE - NOMS DES MATR_ASSE
    !                VEELEM - NOMS DES VECT_ELEM
    !                VEASSE - NOMS DES VECT_ASSE
    !                SOLALG - NOMS DES CHAM_NO SOLUTIONS
    !                VALINC - VALEURS SOLUTION INCREMENTALE
    ! IN  TYVARI : TYPE DE LA VARIABLE
    !                OU  LONUTI - NOMBRE DE VAR. STOCKEES
    ! OUT VALI   : INDEX OU EST STOCKE LE NOM DE LA VARIABLE DANS UNE
    !              VARIABLE CHAPEAU
    !                OU NOMBRE DE VAR. STOCKEES
    !
    ! ----------------------------------------------------------------------
    !
    integer(kind=8), parameter :: zmeelm = 8
    integer(kind=8), parameter :: zmeass = 4
    integer(kind=8), parameter :: zveelm = 12
    integer(kind=8), parameter :: zveass = 19
    integer(kind=8), parameter :: zsolal = 17
    integer(kind=8), parameter :: zvalin = 28
    !
    character(len=8) :: lmeelm(zmeelm), lmeass(zmeass)
    character(len=8) :: lveelm(zveelm), lveass(zveass)
    character(len=8) :: lsolal(zsolal)
    character(len=8) :: lvalin(zvalin)
    !
    character(len=8) :: tyvari
    !
    data lmeelm/'MEDIRI', 'MEMASS', 'MEAMOR', 'MESUIV',&
     &             'MESSTR', 'MEGEOM', 'MEELTC', 'MEELTF'/
    data lmeass/'MERIGI', 'MEMASS', 'MEAMOR', 'MESSTR'/
    !
    data lveelm/'CNDIRI', 'CNBUDI', 'CNDIDO',&
     &             'CNDIPI', 'CNFEDO', 'CNFEPI', 'CNONDP',&
     &             'CNFSDO', 'CNIMPE', 'CNDIDI', 'CNSSTF',&
     &             'CNREFE'/
    data lveass/'CNDIRI', 'CNBUDI', 'CNDIDO',&
     &             'CNDIPI', 'CNFEDO', 'CNFEPI', 'CNONDP',&
     &             'CNFSDO', 'CNIMPE', 'CNDIDI', 'CNSSTF',&
     &             'CNREFE',&
     &             'CNCINE', 'CNSSTR', 'CNDYNA',&
     &             'CNAMOD', 'CNFEXT',&
     &             'CNVISS', 'CNHYST'/
    !
    data lsolal/'DDEPLA', 'DEPDEL', 'DEPOLD', 'DEPPR1', 'DEPPR2',&
     &             'DVITLA', 'VITDEL', 'VITOLD', 'VITPR1', 'VITPR2',&
     &             'DACCLA', 'ACCDEL', 'ACCOLD', 'ACCPR1', 'ACCPR2',&
     &             'DEPSO1', 'DEPSO2'/
    !
    data lvalin/'DEPMOI', 'SIGMOI', 'VARMOI', 'VITMOI', 'ACCMOI',&
     &             'COMMOI', 'DEPPLU', 'SIGPLU', 'VARPLU', 'VITPLU',&
     &             'ACCPLU', 'COMPLU', 'SIGEXT', 'DEPKM1', 'VITKM1',&
     &             'ACCKM1', 'ROMKM1', 'ROMK', 'STRMOI', 'STRPLU',&
     &             'FEXMOI', 'FEXPLU', 'FAMMOI', 'FAMPLU', 'FLIMOI',&
     &             'FLIPLU', 'FNOMOI', 'FNOPLU'/
    !
    ! ----------------------------------------------------------------------
    !
    tyvari = tyvarz
    if (present(tychap_out)) then
        ASSERT(tychap .eq. 'VEASSE')
        ASSERT(vali .ge. 1)
        ASSERT(vali .le. zveass)
        tychap_out = lveass(vali) (1:6)
        ASSERT(tychap_out .eq. tyvari(1:6))
        goto 99
    end if
    vali = -1
    !
    if (tychap .eq. 'MEELEM') then
        if (tyvari .eq. 'LONMAX') then
            vali = zmeelm
        else
            vali = indik8(lmeelm, tyvari, 1, zmeelm)
        end if
    else if (tychap .eq. 'MEASSE') then
        if (tyvari .eq. 'LONMAX') then
            vali = zmeass
        else
            vali = indik8(lmeass, tyvari, 1, zmeass)
        end if
    else if (tychap .eq. 'VEELEM') then
        if (tyvari .eq. 'LONMAX') then
            vali = zveelm
        else
            vali = indik8(lveelm, tyvari, 1, zveelm)
        end if
    else if (tychap .eq. 'VEASSE') then
        if (tyvari .eq. 'LONMAX') then
            vali = zveass
        else
            vali = indik8(lveass, tyvari, 1, zveass)
        end if
    else if (tychap .eq. 'SOLALG') then
        if (tyvari .eq. 'LONMAX') then
            vali = zsolal
        else
            vali = indik8(lsolal, tyvari, 1, zsolal)
        end if
    else if (tychap .eq. 'VALINC') then
        if (tyvari .eq. 'LONMAX') then
            vali = zvalin
        else
            vali = indik8(lvalin, tyvari, 1, zvalin)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
    !
99  continue
    !
    ASSERT(vali .gt. 0)
    !
end subroutine
