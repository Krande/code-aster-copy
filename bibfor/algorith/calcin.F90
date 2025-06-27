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
subroutine calcin(option, max, may, maz, model, &
                  veprj, modx, mody, modz, i, &
                  j, mij)
    implicit none
!
! ROUTINE CALCULANT LA MASSE AJOUTEE SUR LE MODELE THERMIQUE
!  D INTERFACE AINSI QUE LE PREMIER COEFFICIENT D AMORTISSEMENT
!  AJOUTE
! IN :MAX,MAY,MAZ : MATRICES AX ET AY ET AZ
! IN: VEPRJ: PRESSION PROJETEE DUE AU MODE OU CHAMNO J
! IN: MODEL: K2 : CHAINE DISTINGUANT LE TYPE DE MODELISATION
! IN: MODX,MODY,MODZ : DEPLACEMENTS PROJETES
! OUT : MIJ : MASSE AJOUTEE
!-------------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/ddot.h"
    integer(kind=8) :: i, j
    real(kind=8) :: mij
    character(len=*) :: model, option
    character(len=19) :: modx, mody, modz, veprj, max, may, maz
!--------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: imatx, imaty, imatz
    integer(kind=8) :: nbpres
    real(kind=8) :: rx, ry, rz
    real(kind=8), pointer :: vectx(:) => null()
    real(kind=8), pointer :: vecty(:) => null()
    real(kind=8), pointer :: vectz(:) => null()
    real(kind=8), pointer :: vmodx(:) => null()
    real(kind=8), pointer :: vmody(:) => null()
    real(kind=8), pointer :: vmodz(:) => null()
    real(kind=8), pointer :: pres(:) => null()
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(modx//'.VALE', 'L', vr=vmodx)
    call jeveuo(mody//'.VALE', 'L', vr=vmody)
!
    call jeveuo(veprj//'.VALE', 'L', vr=pres)
    call jelira(veprj//'.VALE', 'LONMAX', nbpres)
!
    AS_ALLOCATE(vr=vectx, size=nbpres)
    AS_ALLOCATE(vr=vecty, size=nbpres)
!
! --- RECUPERATION DES DESCRIPTEURS DE MATRICES ASSEMBLEES MAX ET MAY
!
    call mtdscr(max)
    call jeveuo(max(1:19)//'.&INT', 'E', imatx)
    call mtdscr(may)
    call jeveuo(may(1:19)//'.&INT', 'E', imaty)
!
!------MULTIPLICATIONS MATRICE MAX * CHAMNO MODX---------------------
!----------ET MATRICE MAY * CHAMNO MODY------------------------------
!
    call mrmult('ZERO', imatx, vmodx, vectx, 1, &
                .true._1)
    call mrmult('ZERO', imaty, vmody, vecty, 1, &
                .true._1)
!
!--PRODUITS SCALAIRES VECTEURS PRESSION PAR MAX*MODX ET MAY*MODY
!
    b_n = to_blas_int(nbpres)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    rx = ddot(b_n, pres, b_incx, vectx, b_incy)
    b_n = to_blas_int(nbpres)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    ry = ddot(b_n, pres, b_incx, vecty, b_incy)
!
!
!---------------- MENAGE SUR LA VOLATILE ---------------------------
!
!
    AS_DEALLOCATE(vr=vectx)
    AS_DEALLOCATE(vr=vecty)
!
    call detrsd('CHAM_NO', modx)
    call detrsd('CHAM_NO', mody)
!
!
!
!
!
!------ MASSE AJOUTEE = PRESSION*MAX*MODX + PRESSION*MAY*MODY-------
!--------------------------+ PRESSION*MAZ*MODZ  EN 3D---------------
    if (model .eq. '3D') then
!
        call jeveuo(modz//'.VALE', 'L', vr=vmodz)
        AS_ALLOCATE(vr=vectz, size=nbpres)
        call mtdscr(maz)
        call jeveuo(maz(1:19)//'.&INT', 'E', imatz)
        call mrmult('ZERO', imatz, vmodz, vectz, 1, &
                    .true._1)
        b_n = to_blas_int(nbpres)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        rz = ddot(b_n, pres, b_incx, vectz, b_incy)
        AS_DEALLOCATE(vr=vectz)
        call detrsd('CHAM_NO', modz)
        mij = rx+ry+rz
!
    else
        mij = rx+ry
    end if
!
    if ((i .eq. j) .and. (mij .lt. 0) .and. (option .eq. 'MASS_AJOU')) then
        call utmess('A', 'ALGORITH_60')
    end if
    if ((i .eq. j) .and. (mij .lt. 0) .and. (option .eq. 'AMOR_AJOU')) then
        call utmess('A', 'ALGORITH_61')
    end if
!
    call detrsd('CHAM_NO', veprj)
!
!-----------------------------------------------------------------
    call jedema()
end subroutine
