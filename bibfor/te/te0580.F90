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

subroutine te0580(nomopt, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/utmess.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
    character(len=16) :: nomte, nomopt

!-----------------------------------------------------------------------

    character(len=1) :: code
    integer(kind=8) :: jad, itab(8), nbv, iret, k, kpara, mater, icodre(2)
    integer(kind=8) :: imate, idimge, npara, nno, igeom, ndim, ino, ier, ipt, nbpt
    real(kind=8) :: valres(2), valpar(3), vxyz, pr
    character(len=8) :: nompar(3)
    character(len=16) :: nomres(2)
    character(len=24) :: valk(2)
    character(len=32) :: phenom
    character(len=8) :: param
    character(len=8), parameter :: lparam1(2) = (/'PPRESSR', 'PPRESSF'/)
    character(len=8), parameter :: lparam2(6) = ['PPRESSR', 'PPRESSF', 'PFR2D3D', 'PFF2D3D', &
                                                 'PFR1D2D', 'PFF1D2D']
    character(len=8), parameter :: lparam3(2) = (/'PFRCO3D', 'PFFCO3D'/)
!-----------------------------------------------------------------------
! Cette routine realise les calculs elementaires "triviaux" qui ne sont pas
! encore programmes par les elements.
! Par exemple les chargements de Neumann nuls.
!-----------------------------------------------------------------------

    if (nomopt(1:15) .eq. 'CHAR_MECA_PRES_' .or. nomopt(1:15) .eq. 'CHAR_MECA_PRSU_' &
        .or. nomopt(1:15) .eq. 'RIGI_MECA_PRSU_') then
!   ===================================================================================
        do kpara = 1, 2
            param = lparam1(kpara)
            call tecach('NNO', param, 'L', iret, nval=8, itab=itab)
            if (iret .eq. 0) then
                jad = itab(1)
                nbv = itab(2)
                ASSERT(itab(5) .eq. 1 .or. itab(5) .eq. 4)
                if (itab(5) .eq. 1) then
                    do k = 1, nbv
                        if (zr(jad-1+k) .ne. 0.d0) goto 998
                    end do
                else
                    do k = 1, nbv
                        if (zk8(jad-1+k) .ne. '&FOZERO') then
                            call fointe(' ', zk8(jad-1+k), 0, ' ', [0.d0], pr, ier)
                            if (ier .eq. 0 .and. pr .eq. 0.d0) then
                                ! tout va bien ...
                            else
                                goto 998
                            end if
                        end if
                    end do
                end if
            end if
        end do

    elseif (nomopt .eq. 'CHAR_MECA_EPSI_R') then
!   ===================================================================================
        do kpara = 1, 1
            param = 'PEPSINR'
            call tecach('NNO', param, 'L', iret, nval=8, itab=itab)
            if (iret .eq. 0) then
                jad = itab(1)
                nbv = itab(2)
                ASSERT(itab(5) .eq. 1 .or. itab(5) .eq. 4)
                if (itab(5) .eq. 1) then
                    do k = 1, nbv
                        if (zr(jad-1+k) .ne. 0.d0) goto 998
                    end do
                else
                    ASSERT(.false.)
                end if
            end if
        end do

    elseif (nomopt .eq. 'CHAR_MECA_SFCO3D' .or. nomopt .eq. 'CHAR_MECA_SRCO3D' &
            .or. nomopt .eq. 'RIGI_MECA_SFCO3D' .or. nomopt .eq. 'RIGI_MECA_SRCO3D') then
!   ===================================================================================
        do kpara = 1, 2
            param = lparam3(kpara)
            call tecach('NNO', param, 'L', iret, nval=8, itab=itab)
            if (iret .eq. 0) then
                jad = itab(1)
                nbv = itab(2)
                ASSERT(itab(5) .eq. 1 .or. itab(5) .eq. 4)
                ASSERT(mod(nbv, 8) .eq. 0)
                nbpt = nbv/8
                ! on ne conserve que les 6 CMPS FX, FY, ..., MZ
                nbv = 6
                if (itab(5) .eq. 1) then
                    do ipt = 1, nbpt
                        do k = 1, nbv
                            if (zr(jad-1+8*(ipt-1)+k) .ne. 0.d0) goto 998
                        end do
                    end do
                else
                    do ipt = 1, nbpt
                        do k = 1, nbv
                            if (zk8(jad-1+k) .ne. '&FOZERO') then
                                call fointe(' ', zk8(jad-1+8*(ipt-1)+k), 0, ' ', [0.d0], pr, ier)
                                if (ier .eq. 0 .and. pr .eq. 0.d0) then
                                    ! tout va bien ...
                                else
                                    goto 998
                                end if
                            end if
                        end do
                    end do
                end if
            end if
        end do

    elseif (nomopt .eq. 'CALC_G_XFEM' .or. nomopt .eq. 'CALC_G_XFEM_F' &
            .or. nomopt .eq. 'CALC_K_G_XFEM' .or. nomopt .eq. 'CALC_K_G_XFEM_F') then
!   =======================================================================

!       -- le resultat est nul si les forces de bord sont nulles ou
!          si le champ theta est nul.

!       -- on regarde d'abord theta :
        call tecach('OOO', 'PTHETAR', 'L', iret, nval=8, itab=itab)
        ASSERT(iret .eq. 0)
        jad = itab(1)
        nbv = itab(2)
        do k = 1, nbv
            if (zr(jad-1+k) .ne. 0.d0) goto 2
        end do
        goto 999

!       -- on regarde toutes les forces :
2       continue
        do kpara = 1, 6
            param = lparam2(kpara)
            call tecach('NNO', param, 'L', iret, nval=8, itab=itab)
            if (iret .eq. 0) then
                jad = itab(1)
                nbv = itab(2)
                ASSERT(itab(5) .eq. 1 .or. itab(5) .eq. 4)
                if (itab(5) .eq. 1) then
                    do k = 1, nbv
                        if (zr(jad-1+k) .ne. 0.d0) goto 998
                    end do
                else
                    do k = 1, nbv
                        if (zk8(jad-1+k) .ne. '&FOZERO') goto 998
                    end do
                end if
            end if
        end do

    elseif (nomopt .eq. 'AMOR_MECA') then
!   ==========================================================================

!       -- le resultat est nul si les coefficients d'amortissement sont nuls ou absents.

        call jevech('PMATERC', 'L', imate)
        mater = zi(imate)
        call rccoma(mater, 'ELAS', 0, phenom, icodre(1))
        if (icodre(1) .ne. 0) goto 999

        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'

        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno)
        call tecach('ONO', 'PGEOMER', 'L', iret, nval=5, itab=itab)
        igeom = itab(1)
        idimge = itab(2)/nno

        ASSERT(idimge .eq. 2 .or. idimge .eq. 3)

        npara = idimge
        do k = 1, npara
            vxyz = 0.d0
            do ino = 1, nno
                vxyz = vxyz+zr(igeom+idimge*(ino-1)+k-1)
            end do
            valpar(k) = vxyz/nno
        end do

        nomres(1) = 'AMOR_ALPHA'
        nomres(2) = 'AMOR_BETA'
        call rcvalb('RIGI', 1, 1, '+', mater, ' ', phenom, npara, nompar, valpar, 2, &
                    nomres, valres, icodre, 0)

        if (icodre(1) .ne. 0 .and. icodre(2) .ne. 0) goto 999
        if (valres(1) .ne. 0.d0 .or. valres(2) .ne. 0.d0) goto 998

    else
        ASSERT(.false.)
    end if
    goto 999

!   -- erreur :
998 continue
    valk(1) = nomte
    valk(2) = nomopt
    code = 'F'

!   -- le bloc if suivant sera a retirer apres la correction de issue23503
    if (nomopt(1:14) .eq. 'CHAR_MECA_PRES') then
        if (nomte .eq. 'HM_J_AXSE3' .or. nomte .eq. 'HM_J_DPSE3') code = 'F'
    end if

!   -- le bloc if suivant sera a retirer apres la correction de issue23504
    if (nomopt .eq. 'CALC_K_G_XFEM_F') then
        if (nomte .eq. 'MECA_XH_FACE4' .or. nomte .eq. 'MECA_XHT_FACE4' &
            .or. nomte .eq. 'MECA_XT_FACE4' .or. &
            nomte .eq. 'MECA_XH_FACE8' .or. nomte .eq. 'MECA_XHT_FACE8' &
            .or. nomte .eq. 'MECA_XT_FACE8') code = 'A'
    end if

    if (nomopt(1:14) .eq. 'CHAR_MECA_PRES' .and. &
        (lteatt('GRILLE', 'OUI') .or. lteatt('MODELI', 'MMB'))) then
        call utmess(code, 'CALCUL_48')
    else
        call utmess(code, 'CALCUL_42', 2, valk=valk)
    end if

!   -- sortie normale :
999 continue

end subroutine
