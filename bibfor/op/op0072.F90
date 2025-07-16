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
subroutine op0072()
!
!  CALCUL PROJECTION VECTEUR SUR BASE DE RITZ
!
!-----------------------------------------------------------------------
!
    implicit none
!
!
!
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/idensd.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rrlds.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/trlds.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/zdotc.h"
!
    character(len=1) :: typvec
    character(len=8) :: nomres, basemo, vectas, nomtyp, maill1, maill2, k8bid
    character(len=14) :: numgen, numdd1, numdd2
    character(len=16) :: typres, nomcom, typbas, matri2
    character(len=19) :: proch1, proch2, nume2, nomcha
    character(len=24) :: matric, deeq, typeba, valk(24)
    complex(kind=8) :: cbid, dcmplx
    real(kind=8) :: zero
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadesc, iadvec, iamatr, iarefe
    integer(kind=8) :: iavale, ibid, icod, idbase, iddeeq, idvec1, idvec2
    integer(kind=8) :: idvec3, idvec4, idvect, iret, j
    integer(kind=8) :: n0, n1, n2, n3, n4, nbid, nbmode, tmod(1)
    integer(kind=8) :: neq
    real(kind=8) :: bid, ebid, pij
    integer(kind=8), pointer :: smde(:) => null()
    character(len=24), pointer :: refa(:) => null()
    character(len=24), pointer :: refe(:) => null()
    integer(kind=8), pointer :: nequ(:) => null()
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    call jemarq()
    call infmaj()
    zero = 0.d0
!
! --- RECUPERATION DES ARGUMENTS DE LA COMMANDE
!
    call getres(nomres, typres, nomcom)
    call getvid(' ', 'NUME_DDL_GENE', scal=numgen, nbret=n0)
    call getvid(' ', 'VECT_ASSE', scal=vectas, nbret=n1)
    call getvid(' ', 'VECT_ASSE_GENE', scal=vectas, nbret=n3)
    call getvid(' ', 'BASE', scal=basemo, nbret=n4)
    call getvtx(' ', 'TYPE_VECT', scal=nomtyp, nbret=n2)
    call gettco(basemo, typbas, .true._1)
!
! --- RECUPERATION DU NB DE MODES
!
    call rsorac(basemo, 'LONUTI', ibid, bid, k8bid, &
                cbid, ebid, 'ABSOLU', tmod, 1, &
                nbid)
    nbmode = tmod(1)
!
!
    call jeveuo(numgen//'.SMOS.SMDE', 'L', vi=smde)
!
! --- VERIFICATION DE LA CONFORMITE DES NUMEROTATIONS
!     DES MODES ET DU VECTEUR ASSEMBLE
!
    call jeveuo(vectas//'           .VALE', 'L', iadvec)
    call jeveuo(vectas//'           .REFE', 'L', vk24=refe)
    call jelira(vectas//'           .VALE', 'TYPE', cval=typvec)
    call dismoi('TYPE_BASE', basemo, 'RESU_DYNA', repk=typeba, arret='C')
!
    if (typbas(1:9) .eq. 'MODE_MECA') then
        proch1 = refe(2) (1:19)
        call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric, arret='C')
        if (typeba(1:1) .eq. ' ') then
            call exisd('MATR_ASSE', matric, iret)
            if (iret .ne. 0) then
                call dismoi('NUME_EQUA', matric, 'MATR_ASSE', repk=proch2)
            else
                proch2 = proch1
            end if
        else
            call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=nume2)
            proch2 = nume2(1:14)//'.NUME'
        end if
!
        if (.not. idensd('NUME_EQUA', proch1, proch2)) then
            call utmess('I', 'ALGORITH9_41')
        end if
!
    else if (typbas(1:9) .eq. 'MODE_GENE') then
        numdd1 = refe(2) (1:14)
        proch1 = numdd1//'.NUME'
        call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric)
        matri2 = matric(1:16)
        call jeveuo(matri2//'   .REFA', 'L', vk24=refa)
        numdd2 = refa(2) (1:14)
        if (numdd1 .ne. numdd2) then
            call utmess('I', 'ALGORITH9_41')
        end if
    end if
!
! --- RECUPERATION DU NOMBRE D'EQUATIONS DU SYSTEME PHYSIQUE
!
!
    if ((typbas(1:9) .eq. 'MODE_MECA')) then
        call dismoi('NB_EQUA', vectas, 'CHAM_NO', repi=neq)
        deeq = proch1//'.DEEQ'
    else if (typbas(1:9) .eq. 'MODE_GENE') then
        call jeveuo(numdd1//'.NUME.NEQU', 'L', vi=nequ)
        neq = nequ(1)
        deeq = numdd1//'.NUME.DEEQ'
    end if
!
    call jeveuo(deeq, 'L', iddeeq)
!
! --- CREATION DE L OBJET VECT_GENE RESULTAT
!
    if (typvec .eq. 'R') then
        call wkvect(nomres//'           .VALE', 'G V R', nbmode, iavale)
    else
        call wkvect(nomres//'           .VALE', 'G V C', nbmode, iavale)
        call wkvect('&&OP0072.VECTASC1', 'V V C', neq, idvec3)
    end if
    call wkvect(nomres//'           .REFE', 'G V K24', 2, iarefe)
    call wkvect(nomres//'           .DESC', 'G V I', 3, iadesc)
    call jeecra(nomres//'           .DESC', 'DOCU', cval='VGEN')
!
! --- REMPLISSAGE DU .REFE ET .VALE
!
    zk24(iarefe) = basemo
    zk24(iarefe+1) = numgen//'.NUME     '
    zi(iadesc) = 1
    zi(iadesc+1) = nbmode
!
!   LE STOCKAGE EST-IL DIAGONAL ?
    if (smde(4) .eq. smde(1)) then
        zi(iadesc+2) = 1
    else
        zi(iadesc+2) = 2
    end if
    call wkvect('&&OP0072.BASEMO', 'V V R', nbmode*neq, idbase)
!
    if ((typbas(1:9) .eq. 'MODE_MECA')) then
!       --- VERIFIER QUE LES MAILLAGES DU CHAMP A PROJETER
!         - LES DEFORMEES MODALES SONT IDENTIQUES
!         - 1. MAILLAGE DE REFERENCE POUR LA BASE
        call rsexch('F', basemo, 'DEPL', 1, nomcha, &
                    iret)
        call dismoi('NOM_MAILLA', nomcha, 'CHAM_NO', repk=maill1)
!       - 2. MAILLAGE DE REFERENCE POUR LE CHAM_NO
        call dismoi('NOM_MAILLA', vectas, 'CHAM_NO', repk=maill2)
        if (maill1 .ne. maill2) then
            valk(1) = proch2
            valk(2) = maill2
            valk(3) = proch1
            valk(4) = maill1
            call utmess('F', 'ALGORITH12_62', nk=4, valk=valk)
        end if
    end if
!
! --- CONVERSION DE BASEMO A LA NUMEROTATION NU
    call copmod(basemo, numer=proch1, bmodr=zr(idbase))
!
    if (nomtyp(1:4) .eq. 'FORC') then
!
! --- PROJECTION D UN VECTEUR DE TYPE FORCE
!
        call wkvect('&&OP0072.VECTASSE', 'V V R', neq, idvect)
        do i = 1, nbmode
!
! --------- RECOPIE DU IEME MODE
!
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(idbase+(i-1)*neq), b_incx, zr(idvect), b_incy)
!
! ------- PRODUIT SCALAIRE VECTASS * MODE
!
            if (typvec .eq. 'R') then
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(iavale+i-1) = ddot(b_n, zr(idvect), b_incx, zr(iadvec), b_incy)
            else
                do j = 1, neq
                    zc(idvec3+j-1) = dcmplx(zr(idvect+j-1), zero)
                end do
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zc(iavale+i-1) = zdotc(b_n, zc(idvec3), b_incx, zc(iadvec), b_incy)
            end if
        end do
    else
!
! --- PROJECTION D UN VECTEUR DE TYPE DEPL OU VITE
!
        call wkvect('&&OP0072.VECTASS1', 'V V R', neq, idvec1)
        call wkvect('&&OP0072.VECTASS2', 'V V R', neq, idvec2)
        if (typvec .eq. 'C') then
            call wkvect('&&OP0072.VECTASC2', 'V V R', neq, idvec4)
        end if
        call wkvect('&&OP0072.MATRNORM', 'V V R', nbmode*nbmode, iamatr)
!
! ----- CALCUL DE TMODE*MODE
!
        do i = 1, nbmode
!
! ----- RECOPIE DU IEME MODE
!
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(idbase+(i-1)*neq), b_incx, zr(idvec1), b_incy)
!
!-------- PRODUIT SCALAIRE MODE(I)*MODE(J)
!
            do j = i, nbmode
!
! ------- RECOPIE DU JEME MODE
!
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(idbase+(j-1)*neq), b_incx, zr(idvec2), b_incy)
!
! --------- PRODUIT SCALAIRE MODE(I)*MODE(J)
!
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                pij = ddot(b_n, zr(idvec1), b_incx, zr(idvec2), b_incy)
                zr(iamatr+i+(j-1)*nbmode-1) = pij
                zr(iamatr+j+(i-1)*nbmode-1) = pij
            end do
        end do
!
! ----- CALCUL DE LA PROJECTION
!
        do i = 1, nbmode
!
! ------- RECOPIE DU IEME MODE
!
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(idbase+(i-1)*neq), b_incx, zr(idvec1), b_incy)
!
! ------- PRODUIT SCALAIRE VECTASS * MODE
!
            if (typvec .eq. 'R') then
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(idvec2+i-1) = ddot(b_n, zr(idvec1), b_incx, zr(iadvec), b_incy)
            else
                do j = 1, neq
                    zc(idvec3+j-1) = dcmplx(zr(idvec1+j-1), zero)
                end do
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zc(idvec4+i-1) = zdotc(b_n, zc(idvec3), b_incx, zc(iadvec), b_incy)
            end if
        end do
!
! ----- FACTORISATION ET RESOLUTION SYSTEME
!
        call trlds(zr(iamatr), nbmode, nbmode, icod)
        if (icod .ne. 0) then
            call utmess('F', 'ALGORITH9_42')
        end if
        if (typvec .eq. 'R') then
            call rrlds(zr(iamatr), nbmode, nbmode, zr(idvec2), 1)
            b_n = to_blas_int(nbmode)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(idvec2), b_incx, zr(iavale), b_incy)
        else
            call utmess('F', 'ALGORITH9_1')
        end if
    end if
!
    call jedema()
end subroutine
