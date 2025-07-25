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

subroutine mdchan(nlcase, ioc, iliai, mdgene, typnum, &
                  repere, xjeu)
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/angvx.h"
#include "asterfort/angvxy.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedetr.h"
#include "asterfort/normev.h"
#include "asterfort/orient.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/nlget.h"
#include "asterfort/nlsav.h"
    integer(kind=8) :: ioc, iliai
    real(kind=8) :: xjeu
    character(len=8) :: repere, sd_nl
    character(len=*) :: nlcase
    character(len=16) :: typnum
    character(len=24) :: mdgene
!
!     ROUTINE APPELEE PAR MDCHOC
!     RECHERCHE DES ANGLES NAUTIQUES
!
! IN  : MOTFAC : 'DIS_CHOC', 'FLAMBAGE', 'ANTI_SISM'
! IN  : IOC    : NUMERO D'OCCURENCE
! IN  : ILIAI  : NUMERO DE LA LIAISON TRAITEE
! IN  : MDGENE : MODELE GENERALISE
! IN  : TYPNUM : TYPE DE LA NUMEROTATION
! IN  : REPERE : REPERE DU NOEUD DE CHOC = 'GLOBAL' OU 'LOCAL'
! IN  : XJEU   : JEU INITIAL
! IN  : NBNLI  : DIMENSION DES TABLEAUX (NBCHOC+NBSISM+NBFLAM)
! IN  : NOECHO : (ILIAI,9) = TYPE D'OBSTACLE
! OUT : PARCHO : PARAMETRE DE CHOC:
!                PARCHO(ILIAI,17)= SIN A
!                PARCHO(ILIAI,18)= COS A
!                PARCHO(ILIAI,19)= SIN B
!                PARCHO(ILIAI,20)= COS B
!                PARCHO(ILIAI,21)= SIN G
!                PARCHO(ILIAI,22)= COS G
!     ------------------------------------------------------------------
    integer(kind=8) :: n1, jnorm
    real(kind=8) :: txloc(3), tzloc(3), tyloc(3), ang(3), alpha, beta
    real(kind=8) :: normx(3), normy(3), angl, rnorm, rad, eps
    real(kind=8), pointer :: coor_no1(:) => null()
    real(kind=8), pointer :: coor_no2(:) => null()
    character(len=8) :: obst_typ
    character(len=16) :: motfac
!     ------------------------------------------------------------------
!
    rad = r8dgrd()
    eps = r8prem()
    sd_nl = '&&OP29NL'
    motfac = nlcase
!
    if (motfac .eq. 'DIS_CHOC' .or. motfac .eq. 'FLAMBAGE') then
!          ------------------------------------------
        call getvr8('COMPORTEMENT', 'NORM_OBST', iocc=ioc, nbval=3, vect=txloc, &
                    nbret=n1)
        call getvr8('COMPORTEMENT', 'ANGL_VRIL', iocc=ioc, scal=angl, nbret=n1)
        call nlget(sd_nl, _OBST_TYP, iocc=iliai, kscal=obst_typ)
!
        if (n1 .ne. 0) then
            if (typnum .eq. 'NUME_DDL_SDASTER' .or. repere .eq. 'GLOBAL') then
                call angvx(txloc, alpha, beta)
                call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(alpha), cos(alpha)])
                call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(beta), cos(beta)])
            else
                call wkvect('&&MDCHAN.NORM', 'V V R', 3, jnorm)
                zr(jnorm) = txloc(1)
                zr(jnorm+1) = txloc(2)
                zr(jnorm+2) = txloc(3)
                call orient(mdgene, repere, jnorm, 1, normx, &
                            0)
                call angvx(normx, alpha, beta)
                call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(alpha), cos(alpha)])
                call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(beta), cos(beta)])
                call jedetr('&&MDCHAN.NORM')
            end if
            call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=iliai, rvect=[sin(angl*rad), cos(angl*rad)])
!
        else if (obst_typ .eq. 'BI_PLANY') then
            call nlget(sd_nl, _COOR_NO1, iocc=iliai, vr=coor_no1)
            call nlget(sd_nl, _COOR_NO2, iocc=iliai, vr=coor_no2)
            tyloc(1) = (coor_no2(1)-coor_no1(1))
            tyloc(2) = (coor_no2(2)-coor_no1(2))
            tyloc(3) = (coor_no2(3)-coor_no1(3))
            if (typnum .eq. 'NUME_DDL_SDASTER' .or. repere .eq. 'GLOBAL') then
                call angvxy(txloc, tyloc, ang)
                call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(ang(1)), cos(ang(1))])
                call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(ang(2)), cos(ang(2))])
                call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=iliai, rvect=[sin(ang(3)), cos(ang(3))])
            else
                call wkvect('&&MDCHAN.NORM', 'V V R', 3, jnorm)
                zr(jnorm) = txloc(1)
                zr(jnorm+1) = txloc(2)
                zr(jnorm+2) = txloc(3)
                call orient(mdgene, repere, jnorm, 1, normx, &
                            0)
                zr(jnorm) = tyloc(1)
                zr(jnorm+1) = tyloc(2)
                zr(jnorm+2) = tyloc(3)
                call orient(mdgene, repere, jnorm, 1, normy, &
                            0)
                call angvxy(normx, normy, ang)
                call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(ang(1)), cos(ang(1))])
                call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(ang(2)), cos(ang(2))])
                call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=iliai, rvect=[sin(ang(3)), cos(ang(3))])
                call jedetr('&&MDCHAN.NORM')
            end if
!
        else if (obst_typ .eq. 'BI_PLANZ') then
            call nlget(sd_nl, _COOR_NO1, iocc=iliai, vr=coor_no1)
            call nlget(sd_nl, _COOR_NO2, iocc=iliai, vr=coor_no2)
            tzloc(1) = (coor_no2(1)-coor_no1(1))
            tzloc(2) = (coor_no2(2)-coor_no1(2))
            tzloc(3) = (coor_no2(3)-coor_no1(3))
            call provec(tzloc, txloc, tyloc)
            if (typnum .eq. 'NUME_DDL_SDASTER' .or. repere .eq. 'GLOBAL') then
                call angvxy(txloc, tyloc, ang)
                call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(ang(1)), cos(ang(1))])
                call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(ang(2)), cos(ang(2))])
                call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=iliai, rvect=[sin(ang(3)), cos(ang(3))])
            else
                call wkvect('&&MDCHAN.NORM', 'V V R', 3, jnorm)
                zr(jnorm) = txloc(1)
                zr(jnorm+1) = txloc(2)
                zr(jnorm+2) = txloc(3)
                call orient(mdgene, repere, jnorm, 1, normx, &
                            0)
                zr(jnorm) = tzloc(1)
                zr(jnorm+1) = tzloc(2)
                zr(jnorm+2) = tzloc(3)
                call orient(mdgene, repere, jnorm, 1, normy, &
                            0)
                call angvxy(normx, normy, ang)
                call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(ang(1)), cos(ang(1))])
                call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(ang(2)), cos(ang(2))])
                call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=iliai, rvect=[sin(ang(3)), cos(ang(3))])
                call jedetr('&&MDCHAN.NORM')
            end if
!
        else
            angl = 0.d0
            if (typnum .eq. 'NUME_DDL_SDASTER' .or. repere .eq. 'GLOBAL') then
                call angvx(txloc, alpha, beta)
                call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(alpha), cos(alpha)])
                call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(beta), cos(beta)])
            else
                call wkvect('&&MDCHAN.NORM', 'V V R', 3, jnorm)
                zr(jnorm) = txloc(1)
                zr(jnorm+1) = txloc(2)
                zr(jnorm+2) = txloc(3)
                call orient(mdgene, repere, jnorm, 1, normx, &
                            0)
                call angvx(normx, alpha, beta)
                call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(alpha), cos(alpha)])
                call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(beta), cos(beta)])
                call jedetr('&&MDCHAN.NORM')
            end if
            call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=iliai, rvect=[sin(angl*rad), cos(angl*rad)])
        end if
!
    else if (motfac .eq. 'ANTI_SISM') then
!            ---------------------
!
        call nlsav(sd_nl, _DIST_NO1, 1, iocc=iliai, rscal=sqrt(xjeu)/2.d0)
        call nlsav(sd_nl, _DIST_NO2, 1, iocc=iliai, rscal=sqrt(xjeu)/2.d0)
!
! --- VECTEUR NOEUD1 VERS NOEUD2
        call nlget(sd_nl, _COOR_NO1, iocc=iliai, vr=coor_no1)
        call nlget(sd_nl, _COOR_NO2, iocc=iliai, vr=coor_no2)
        tyloc(1) = (coor_no2(1)-coor_no1(1))
        call nlget(sd_nl, _COOR_NO1, iocc=iliai, vr=coor_no1)
        call nlget(sd_nl, _COOR_NO2, iocc=iliai, vr=coor_no2)
        tyloc(2) = (coor_no2(2)-coor_no1(2))
        call nlget(sd_nl, _COOR_NO1, iocc=iliai, vr=coor_no1)
        call nlget(sd_nl, _COOR_NO2, iocc=iliai, vr=coor_no2)
        tyloc(3) = (coor_no2(3)-coor_no1(3))
        call normev(tyloc, rnorm)
        if (abs(rnorm) .le. eps) then
            call utmess('F', 'ALGORITH5_26')
        end if
!
! --- DETERMINATION DES AXES LOCAUX
        if (abs(tyloc(3)) .le. abs(tyloc(1)) .and. abs(tyloc(3)) .le. abs(tyloc(2))) then
            tzloc(1) = -tyloc(2)
            tzloc(2) = tyloc(1)
            tzloc(3) = 0.d0
        elseif (abs(tyloc(2)) .le. abs(tyloc(1)) .and. abs(tyloc(2)) &
                .le. abs(tyloc(3))) then
            tzloc(1) = -tyloc(3)
            tzloc(2) = 0.d0
            tzloc(3) = tyloc(1)
        else
            tzloc(1) = 0.d0
            tzloc(2) = -tyloc(3)
            tzloc(3) = tyloc(2)
        end if
        call provec(tyloc, tzloc, txloc)
        if (typnum .eq. 'NUME_DDL_SDASTER' .or. repere .eq. 'GLOBAL') then
            call angvxy(txloc, tyloc, ang)
            call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(ang(1)), cos(ang(1))])
            call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(ang(2)), cos(ang(2))])
            call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=iliai, rvect=[sin(ang(3)), cos(ang(3))])
        else
            call wkvect('&&MDCHAN.NORM', 'V V R', 3, jnorm)
            zr(jnorm) = txloc(1)
            zr(jnorm+1) = txloc(2)
            zr(jnorm+2) = txloc(3)
            call orient(mdgene, repere, jnorm, 1, normx, &
                        0)
            zr(jnorm) = tyloc(1)
            zr(jnorm+1) = tyloc(2)
            zr(jnorm+2) = tyloc(3)
            call orient(mdgene, repere, jnorm, 1, normy, &
                        0)
            call angvxy(normx, normy, ang)
            call nlsav(sd_nl, _SINCOS_ANGLE_A, 2, iocc=iliai, rvect=[sin(ang(1)), cos(ang(1))])
            call nlsav(sd_nl, _SINCOS_ANGLE_B, 2, iocc=iliai, rvect=[sin(ang(2)), cos(ang(2))])
            call nlsav(sd_nl, _SINCOS_ANGLE_G, 2, iocc=iliai, rvect=[sin(ang(3)), cos(ang(3))])
            call jedetr('&&MDCHAN.NORM')
        end if
!
    end if
!
end subroutine
