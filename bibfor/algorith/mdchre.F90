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

subroutine mdchre(nlcase, ioc, iliai, mdgene, typnum, &
                  repere, lnoue2)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/orient.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/nlget.h"
#include "asterfort/nlsav.h"
!
    integer(kind=8) :: ioc, iliai
    aster_logical :: lnoue2
    character(len=8) :: repere, sd_nl
    character(len=*) :: nlcase
    character(len=16) :: typnum
    character(len=24) :: mdgene
!
!     ROUTINE APPELEE PAR MDCHOC
!     TRAITEMENT DU REPERE
!
! IN  : MOTFAC : 'CHOC', 'FLAMBAGE', 'ANTI_SISM'
! IN  : IOC    : NUMERO D'OCCURENCE
! IN  : ILIAI  : NUMERO DE LA LIAISON TRAITEE
! IN  : MDGENE : MODELE GENERALISE
! IN  : TYPNUM : TYPE DE LA NUMEROTATION
! OUT : REPERE : REPERE DU NOEUD DE CHOC = 'GLOBAL' OU 'LOCAL'
! IN  : NBNLI  : DIMENSION DES TABLEAUX (NBCHOC+NBSISM+NBFLAM)
! OUT : PARCHO : PARAMETRE DE CHOC:
!                PARCHO(ILIAI,14)= COOR ORIGINE OBSTACLE X REP GLOBAL
!                PARCHO(ILIAI,15)= COOR ORIGINE OBSTACLE Y REP GLOBAL
!                PARCHO(ILIAI,16)= COOR ORIGINE OBSTACLE Z REP GLOBAL
!                PARCHO(ILIAI,45)= NORMALE X
!                PARCHO(ILIAI,46)= NORMALE Y
!                PARCHO(ILIAI,47)= NORMALE Z
! IN  : LNOUE2 : CHOC DEFINIT ENTRE 2 NOEUDS
!     ------------------------------------------------------------------
    integer(kind=8) :: n1, iret, jcoord
    real(kind=8) :: tempo(3), dircho(3), coord(3), txno
    character(len=16) :: motfac, obstyp
    character(len=24) :: mdssno
    real(kind=8), pointer :: coor_no1(:) => null()
    real(kind=8), pointer :: coor_no2(:) => null()
    real(kind=8), pointer :: ob_orig(:) => null()
!     ------------------------------------------------------------------
!
    sd_nl = '&&OP29NL'
    n1 = 0
    repere = '????????'
    motfac = nlcase
!
    if (motfac .eq. 'DIS_CHOC' .or. motfac .eq. 'FLAMBAGE') then
!          ------------------------------------------
        call getvtx('COMPORTEMENT', 'REPERE', iocc=ioc, nbval=0, nbret=n1)
        if (n1 .eq. 0) then
            repere = 'GLOBAL'
        else
            call getvtx('COMPORTEMENT', 'REPERE', iocc=ioc, scal=repere, nbret=n1)
        end if
        call getvr8('COMPORTEMENT', 'ORIG_OBST', iocc=ioc, scal=tempo(1), nbret=n1)
    end if
!
    n1 = -n1
    if (n1 .eq. 3) then
        call getvr8('COMPORTEMENT', 'ORIG_OBST', iocc=ioc, nbval=3, vect=tempo, &
                    nbret=n1)
        if (typnum(1:16) .eq. 'NUME_DDL_SDASTER') then
            call nlsav(sd_nl, _COOR_ORIGIN_OBSTACLE, 3, iocc=iliai, rvect=tempo)
        else
            mdssno = mdgene(1:14)//'.MODG.SSNO'
            if (repere(1:6) .eq. 'GLOBAL') then
                call nlsav(sd_nl, _COOR_ORIGIN_OBSTACLE, 3, iocc=iliai, rvect=tempo)
            else
                call jenonu(jexnom(mdssno, repere), iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH5_39')
                end if
                call wkvect('&&MDCHOC.COORDO', 'V V R', 3, jcoord)
                zr(jcoord) = tempo(1)
                zr(jcoord+1) = tempo(2)
                zr(jcoord+2) = tempo(3)
                call orient(mdgene, repere, jcoord, 1, coord, &
                            1)
                call nlsav(sd_nl, _COOR_ORIGIN_OBSTACLE, 3, iocc=iliai, rvect=coord)
                call jedetr('&&MDCHOC.COORDO')
            end if
        end if
    else
        call nlget(sd_nl, _COOR_NO1, iocc=iliai, vr=coor_no1)
        call nlget(sd_nl, _OBST_TYP, iocc=iliai, kscal=obstyp)
        if (obstyp(1:2) .eq. 'BI') then
            call nlget(sd_nl, _COOR_NO2, iocc=iliai, vr=coor_no2)
            call nlsav(sd_nl, _COOR_ORIGIN_OBSTACLE, 3, iocc=iliai, &
                       rvect=[(coor_no1(1)+coor_no2(1))/2.d0, &
                              (coor_no1(2)+coor_no2(2))/2.d0, &
                              (coor_no1(3)+coor_no2(3))/2.d0])
        else
            call nlsav(sd_nl, _COOR_ORIGIN_OBSTACLE, 3, iocc=iliai, &
                       rvect=[(coor_no1(1))/2.d0, &
                              (coor_no1(2))/2.d0, &
                              (coor_no1(3))/2.d0])
        end if
    end if
!
    if (lnoue2) then
        call nlget(sd_nl, _COOR_NO1, iocc=iliai, vr=coor_no1)
        call nlget(sd_nl, _COOR_NO2, iocc=iliai, vr=coor_no2)
        dircho(1) = coor_no1(1)-coor_no2(1)
        dircho(2) = coor_no1(2)-coor_no2(2)
        dircho(3) = coor_no1(3)-coor_no2(3)
    else
        call nlget(sd_nl, _COOR_NO1, iocc=iliai, vr=coor_no1)
        call nlget(sd_nl, _COOR_ORIGIN_OBSTACLE, iocc=iliai, vr=ob_orig)
        dircho(1) = coor_no1(1)-ob_orig(1)
        dircho(2) = coor_no1(2)-ob_orig(2)
        dircho(3) = coor_no1(3)-ob_orig(3)
    end if
!
    txno = sqrt(dircho(1)**2+dircho(2)**2+dircho(3)**2)
    if (txno .lt. r8prem()) txno = 1.d0
!
! --- DEBUG : UN TRAVAIL DOIT ETRE FAIT SI TXNO = 0.
!
    call nlsav(sd_nl, _NORMAL_VECTOR, 3, iocc=iliai, &
               rvect=[dircho(1)/txno, dircho(2)/txno, dircho(3)/txno])
!
end subroutine
