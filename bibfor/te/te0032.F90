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
subroutine te0032(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dxqfor.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxroep.h"
#include "asterfort/dxtfor.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
!
    character(len=16) :: option, nomte
!     IN  OPTION : NOM DE L'OPTION A CALCULER
!     IN  NOMTE  : NOM DU TYPE_ELEMENT
!     -----------------------------------------------------------------
!     CALCUL DE PRESSION SUR LES ELEMENTS DKT, DST, DKQ, DSQ ET Q4G
!         OPTIONS TRAITEES   ==>   CHAR_MECA_FRCO3D
!                                  CHAR_MECA_FFCO3D
!                                  CHAR_MECA_PRES_R
!                                  CHAR_MECA_PRES_F
!                                  CHAR_MECA_PESA_R
!     -----------------------------------------------------------------
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano
    integer(kind=8) :: i, j, ier, iplan, jgeom, jcoqu, jvecg, jpres, itemps
    integer(kind=8) :: iadzi, iazk24, lpesa
    real(kind=8) :: pgl(3, 3), xyzl(3, 4), pglo(3), ploc(3)
    real(kind=8) :: vecl(24), for(6, 4), for2(6, 4), rho, epais
    real(kind=8) :: undemi
    real(kind=8) :: valpar(4), dist, excent, pr
    aster_logical :: global, locapr
    character(len=8) :: nompar(4), moplan, nomail
    character(len=24) :: valk
! DEB ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
!
    undemi = 0.5d0
    iplan = 0
!
    call jevech('PGEOMER', 'L', jgeom)
    call jevech('PCACOQU', 'L', jcoqu)
    call jevech('PVECTUR', 'E', jvecg)
!
    if (nno .eq. 3) then
        call dxtpgl(zr(jgeom), pgl)
    else if (nno .eq. 4) then
        call dxqpgl(zr(jgeom), pgl)
    end if
    call utpvgl(nno, 3, pgl, zr(jgeom), xyzl)
!
! --- CAS DES CHARGEMENTS DE FORME REEL
    if (option .eq. 'CHAR_MECA_PRES_R') then
!         ------------------------------
        global = .false.
        call jevech('PPRESSR', 'L', jpres)
        do j = 1, nno
            do i = 1, 6
                for(i, j) = 0.d0
                for2(i, j) = 0.d0
            end do
!----------------------------------------------------------------------
!           LE SIGNE MOINS CORRESPOND A LA CONVENTION :
!              UNE PRESSION POSITIVE PROVOQUE UN GONFLEMENT
!----------------------------------------------------------------------
            for(3, j) = -zr(jpres+j-1)
        end do
!
    else if (option .eq. 'CHAR_MECA_FRCO3D') then
!              ------------------------------
        call jevech('PFRCO3D', 'L', jpres)
        global = abs(zr(jpres+6)) .lt. 1.d-3
        locapr = abs(zr(jpres+6)-3.d0) .lt. 1.d-3
        if (global) then
            call utpvgl(1, 6, pgl, zr(jpres), for(1, 1))
            call utpvgl(1, 6, pgl, zr(jpres+8), for(1, 2))
            call utpvgl(1, 6, pgl, zr(jpres+16), for(1, 3))
            if (nno .eq. 4) then
                call utpvgl(1, 6, pgl, zr(jpres+24), for(1, 4))
            end if
        elseif (locapr) then
            for(:, :) = 0.d0
            do j = 1, nno
                for(3, j) = -zr(jpres-1+8*(j-1)+3)
            end do
        else
            do j = 1, nno
                do i = 1, 5
                    for(i, j) = zr(jpres-1+8*(j-1)+i)
                end do
                for(6, j) = 0.d0
            end do
        end if
        iplan = nint(zr(jpres+7))
!
! --- CAS DES CHARGEMENTS DE FORME FONCTION
!
    else if (option .eq. 'CHAR_MECA_PRES_F') then
!              ------------------------------
        call jevech('PPRESSF', 'L', jpres)
        if (zk8(jpres) .eq. '&FOZERO') goto 999
        call jevech('PINSTR', 'L', itemps)
        valpar(4) = zr(itemps)
        nompar(4) = 'INST'
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        do j = 0, nno-1
            valpar(1) = zr(jgeom+3*j)
            valpar(2) = zr(jgeom+3*j+1)
            valpar(3) = zr(jgeom+3*j+2)
            call fointe('FM', zk8(jpres), 4, nompar, valpar, &
                        pr, ier)
            if (pr .ne. 0.d0) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                valk = nomail
                call utmess('F', 'ELEMENTS4_92', sk=valk)
            end if
        end do
        goto 999
!
    else if (option .eq. 'CHAR_MECA_FFCO3D') then
!              ------------------------------
        call jevech('PFFCO3D', 'L', jpres)
        call jevech('PINSTR', 'L', itemps)
        valpar(4) = zr(itemps)
        nompar(4) = 'INST'
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        global = zk8(jpres+6) .eq. 'GLOBAL'
        locapr = zk8(jpres+6) .eq. 'LOCAL_PR'
        moplan = zk8(jpres+7)
        if (moplan .eq. 'SUP') then
            iplan = 1
        else if (moplan .eq. 'INF') then
            iplan = -1
        else if (moplan .eq. 'MOY') then
            iplan = 2
        end if
!
        if (global) then
!          REPERE GLOBAL
! --       LECTURE DES INTERPOLATIONS DE FX, FY, FZ, MX, MY, MZ
!
            do j = 0, nno-1
                valpar(1) = zr(jgeom+3*j)
                valpar(2) = zr(jgeom+3*j+1)
                valpar(3) = zr(jgeom+3*j+2)
!------------------------------------------------------
!  PAS DE CHANGEMENT DE SIGNE POUR LES FORCES REPARTIES
!------------------------------------------------------
                call fointe('FM', zk8(jpres), 4, nompar, valpar, &
                            for2(1, j+1), ier)
                call fointe('FM', zk8(jpres+1), 4, nompar, valpar, &
                            for2(2, j+1), ier)
                call fointe('FM', zk8(jpres+2), 4, nompar, valpar, &
                            for2(3, j+1), ier)
                call fointe('FM', zk8(jpres+3), 4, nompar, valpar, &
                            for2(4, j+1), ier)
                call fointe('FM', zk8(jpres+4), 4, nompar, valpar, &
                            for2(5, j+1), ier)
                call fointe('FM', zk8(jpres+5), 4, nompar, valpar, &
                            for2(6, j+1), ier)
            end do
!
            call utpvgl(1, 6, pgl, for2(1, 1), for(1, 1))
            call utpvgl(1, 6, pgl, for2(1, 2), for(1, 2))
            call utpvgl(1, 6, pgl, for2(1, 3), for(1, 3))
            if (nno .eq. 4) then
                call utpvgl(1, 6, pgl, for2(1, 4), for(1, 4))
            end if
!
        else if (locapr) then
! --        REPERE LOCAL - CAS D UNE PRESSION
! --        LECTURE DES INTERPOLATIONS DE LA PRESSION PRES
!
            do j = 0, nno-1
                valpar(1) = zr(jgeom+3*j)
                valpar(2) = zr(jgeom+3*j+1)
                valpar(3) = zr(jgeom+3*j+2)
                call fointe('FM', zk8(jpres+2), 4, nompar, valpar, &
                            pr, ier)
!-----------------------------------------------------
!       LE SIGNE MOINS DE FOR(3,J+1) CORRESPOND A LA CONVENTION :
!          UNE PRESSION POSITIVE PROVOQUE UN GONFLEMENT
!-----------------------------------------------------
                for(3, j+1) = -1*pr
                for(1, j+1) = 0.d0
                for(2, j+1) = 0.d0
                for(4, j+1) = 0.d0
                for(5, j+1) = 0.d0
                for(6, j+1) = 0.d0
            end do
!
        else
! --        REPERE LOCAL - CAS DE F1, F2, F3, MF1, MF2
! --        LECTURE DES INTERPOLATIONS DE F1, F2, F3, MF1, MF2
!
            do j = 0, nno-1
                valpar(1) = zr(jgeom+3*j)
                valpar(2) = zr(jgeom+3*j+1)
                valpar(3) = zr(jgeom+3*j+2)
!------------------------------------------------------
!  PAS DE CHANGEMENT DE SIGNE POUR LES FORCES REPARTIES
!------------------------------------------------------
                call fointe('FM', zk8(jpres), 4, nompar, valpar, &
                            for(1, j+1), ier)
                call fointe('FM', zk8(jpres+1), 4, nompar, valpar, &
                            for(2, j+1), ier)
                call fointe('FM', zk8(jpres+2), 4, nompar, valpar, &
                            for(3, j+1), ier)
                call fointe('FM', zk8(jpres+3), 4, nompar, valpar, &
                            for(4, j+1), ier)
                call fointe('FM', zk8(jpres+4), 4, nompar, valpar, &
                            for(5, j+1), ier)
                for(6, j+1) = 0.d0
            end do
        end if
!
    else if (option .eq. 'CHAR_MECA_PESA_R') then
!              ------------------------------
        global = .true.
!
        call dxroep(rho, epais)
        call jevech('PPESANR', 'L', lpesa)
        do i = 1, 3
            pglo(i) = zr(lpesa)*zr(lpesa+i)*rho*epais
        end do
        call utpvgl(1, 3, pgl, pglo, ploc)
        do i = 1, nno
            do j = 1, 3
                for(j, i) = ploc(j)
                for(j+3, i) = 0.d0
            end do
        end do
    end if
!
    if (iplan .ne. 0) then
        epais = zr(jcoqu)
        excent = zr(jcoqu+4)
        if (iplan .eq. 1) then
            dist = excent+undemi*epais
        else if (iplan .eq. -1) then
            dist = excent-undemi*epais
        else if (iplan .eq. 2) then
            dist = excent
        end if
!
        do i = 1, nno
            for(4, i) = for(4, i)-dist*for(2, i)
            for(5, i) = for(5, i)+dist*for(1, i)
        end do
    end if
!
    if (nno .eq. 3) then
        call dxtfor(global, xyzl, pgl, for, vecl)
    else if (nno .eq. 4) then
        call dxqfor(global, xyzl, pgl, for, vecl)
    end if
!
    call utpvlg(nno, 6, pgl, vecl, zr(jvecg))
!
999 continue
end subroutine
