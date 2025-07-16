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
subroutine chor2c(lischa, vecele)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/corich.h"
#include "asterfort/dismoi.h"
#include "asterfort/copisd.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisltc.h"
#include "asterfort/sdchgd.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/dcopy.h"
    character(len=19) :: lischa, vecele
!
! ----------------------------------------------------------------------
!
! TRANSFORMATION D'UNE COLLECTION DE CHAM_NO A VALEURS REELLES EN
! UNE COLLECTION DE CHAM_NO A VALEURS COMPLEXES AVEC PARTIE
! IMAGINAIRE NULLE
!
! ----------------------------------------------------------------------
!
! I/O VECELE : NOM DU VECT_ELEM
!
!
!
!
    character(len=24) :: vachar
    integer(kind=8) :: jvacha
    character(len=19) :: chamno, nume_equa, nume_equa_new, nume_equa_prev
    integer(kind=8) :: jcn, jrefn
    character(len=24) :: resuel, noojb
    character(len=8) :: typech, typsca, nomgd
    integer(kind=8) :: iret, ichar
    integer(kind=8) :: ivec, nbvec, nbvdim, ivale, nbvale
    character(len=4) :: tyresl
    character(len=1) :: typchn
    real(kind=8), pointer :: copie_travail(:) => null()
    character(len=24), pointer :: relr(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- VERIFICATION DU VACHAR
!
    vachar = vecele//'.CHNO'
    call jeexin(vachar, iret)
    ASSERT(iret .ne. 0)
    call jelira(vachar, 'LONMAX', nbvec)
    ASSERT(nbvec .ne. 0)
!
! --- DIMENSIONNEMENT SD DE SAUVEGARDE
!
    call jeveuo(vachar, 'L', jvacha)
    chamno = zk24(jvacha+1-1) (1:19)
    call jeveuo(vecele//'.RELR', 'L', vk24=relr)
    call jelira(chamno//'.VALE', 'LONMAX', nbvdim)
    AS_ALLOCATE(vr=copie_travail, size=nbvdim)
    nume_equa_prev = ' '
!
! --- BOUCLES SUR LES CHAMNO
!
    do ivec = 1, nbvec
!
! ----- NOM DU RESU_ELEM
!
        resuel = relr(ivec)
!
! ----- NOM DU CHAMNO
!
        chamno = zk24(jvacha+ivec-1) (1:19)
!
! ----- RECUPERATION DU NUMERO DE LA CHARGE DU RESU_ELEM
!
        call corich('L', resuel, ichout_=ichar)
        ASSERT((ichar .ne. 0) .and. (ichar .ge. -2))
!
! ----- TYPE DU CHARGEMENT
!
        call dismoi('TYPE_CHAMP', resuel, 'CHAMP', repk=tyresl)
        if (tyresl .eq. 'RESL') then
            call dismoi('TYPE_SCA', resuel, 'RESUELEM', repk=typsca)
            if (typsca .eq. 'R') then
                typchn = 'R'
            else if (typsca .eq. 'C') then
                typchn = 'C'
            else
                ASSERT(.false.)
            end if
        else if (tyresl .eq. 'NOEU') then
            call lisltc(lischa, ichar, typech)
            typchn = 'R'
            if (typech .eq. 'COMP') typchn = 'C'
        else
!
        end if
!
! ----- CONVERSION
!
        if (typchn .eq. 'R') then
            call dismoi("NUME_EQUA", chamno, "CHAM_NO", repk=nume_equa)
            if (nume_equa .ne. nume_equa_prev) then
                noojb = '12345678.NUME000000.PRNO'
                call gnomsd(chamno, noojb, 14, 19)
                noojb(1:8) = chamno(1:8)
                nume_equa_new = noojb(1:19)
                call copisd("NUME_EQUA", "G", nume_equa, nume_equa_new)
                call jeveuo(nume_equa_new//".REFN", "E", jrefn)
                call dismoi("NOM_GD", chamno, "CHAM_NO", repk=nomgd)
                zk24(jrefn-1+2) = nomgd(1:5)//"C"
                nume_equa_prev = nume_equa
            end if
            call jeveuo(chamno//".REFE", "E", jrefn)
            zk24(jrefn-1+2) = nume_equa_new
!
            call jeveuo(chamno//'.VALE', 'L', jcn)
            call jelira(chamno//'.VALE', 'LONMAX', nbvale)
            if (nbvdim .ne. nbvale) then
                ASSERT(.false.)
            end if
!
! ------- SAUVEGARDE DES VALEURS
!
            b_n = to_blas_int(nbvale)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jcn), b_incx, copie_travail, b_incy)
!
! ------- DESTRUCTION CHAMNO A VALEURS REELLES
!
            call jedetr(chamno//'.VALE')
!
! ------- CREATION CHAMNO A VALEURS COMPLEXES
!
            call wkvect(chamno//'.VALE', 'V V C', nbvale, jcn)
            do ivale = 1, nbvale
                zc(jcn+ivale-1) = dcmplx(copie_travail(ivale), 0.d0)
            end do
!
! ------- CHANGEMENT DE LA REFERENCE A LA GRANDEUR
!
            call sdchgd(chamno, 'C')
        end if
    end do
!
    AS_DEALLOCATE(vr=copie_travail)
    call jedema()
end subroutine
