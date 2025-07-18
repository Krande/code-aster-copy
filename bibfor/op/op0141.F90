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
subroutine op0141()
    implicit none
!
!     OPERATEUR DE CALCUL DU MAC DE DEUX BASES MODALES
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/copmod.h"
#include "asterfort/dcapno.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/idensd.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mcmult.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/rsorac.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/tbimex.h"
#include "asterfort/tbimpr.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/zcopy.h"
!
    integer(kind=8) :: n1, n2, n3, ibid, nbmod1, nbmod2, neq, idbas1
    integer(kind=8) :: idbas2, idbas3, idvec3, i, j, nbpara, inom, ityp, ind, imatra
    integer(kind=8) :: idvec1, idvec2, ifm, niv, neq1, llneq2, iret
    integer(kind=8) :: iddl, indv, tmod(1), ieq
    real(kind=8) :: rbid, pij, pii, pjj
    complex(kind=8) :: cbid, dcmplx, ztemp, dconjg
    character(len=1) :: typsca
    character(len=8) :: table, base1, base2, k8b, matras, rep
    character(len=14) :: nu, numdd1, numdd2, numdda
    character(len=16) :: nomcmd, typcon, typba1, typba2, matri1, matri2, depl
    character(len=19) :: matr, pronu1, pronu2, pronua
    character(len=24) :: chamol
    aster_logical :: c1, c2, zcmplx, ieri
    integer(kind=8), pointer :: deeq(:) => null()
    integer(kind=8), pointer :: nllneq1(:) => null()
    blas_int :: b_incx, b_incy, b_n
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
    data depl/'DEPL            '/
!
!     --- RECUPERATION DU RESULTAT ET DU MODE A TRAITER ---
    call jemarq()
!
    call getres(table, typcon, nomcmd)
! CREATION DE LA TABLE CONTENANT LE MAC
    call tbcrsd(table, 'G')
    call titre()
!     ---RECUPERATION DU NIVEAU D'IMPRESSION---
!
    call infmaj()
    call infniv(ifm, niv)
!
! RECUPERATION DE LA MATRICE ASSEMBLEE SI ELLE EXISTE
    call getvid(' ', 'MATR_ASSE', scal=matras, nbret=n1)
    if (n1 .ne. 0) then
! COOL ELLE EXISTE
        call mtdscr(matras)
        matr = matras
        call jeveuo(matr//'.&INT', 'E', imatra)
        call dismoi('NOM_NUME_DDL', matras, 'MATR_ASSE', repk=numdda)
        call dismoi('NB_EQUA', matras, 'MATR_ASSE', repi=neq)
    else
! PAS COOL ELLE EXISTE PAS
        matr = ' '
    end if
!
    ieri = .false.
    call getvtx(' ', 'IERI', scal=rep, nbret=n1)
    if (n1 .eq. 1) then
        if (rep .eq. 'OUI') ieri = .true.
    end if
!
! RECUPERATION DES BASES DE MODES
    call getvid(' ', 'BASE_1', scal=base1, nbret=n2)
    call getvid(' ', 'BASE_2', scal=base2, nbret=n3)
!
    c1 = .false.
    c2 = .false.
!
    call dcapno(base1, depl, 1, chamol)
    call jelira(chamol, 'TYPE', cval=typsca)
    if (typsca .eq. 'C') c1 = .true.
!
    call dcapno(base2, depl, 1, chamol)
    call jelira(chamol, 'TYPE', cval=typsca)
    if (typsca .eq. 'C') c2 = .true.
!
    if (c1 .and. c2) then
        zcmplx = .true.
    else
        zcmplx = .false.
        if (c1 .or. c2) then
            call utmess('F', 'ALGELINE5_71')
        end if
    end if
!
! RECUPERATION DU TYPE ET DU NBRE DE MODES DES BASES
    call gettco(base1, typba1)
    call rsorac(base1, 'LONUTI', 0, rbid, k8b, &
                cbid, rbid, 'ABSOLU', tmod, 1, &
                ibid)
    nbmod1 = tmod(1)
    call gettco(base2, typba2)
    call rsorac(base2, 'LONUTI', 0, rbid, k8b, &
                cbid, rbid, 'ABSOLU', tmod, 1, &
                ibid)
    nbmod2 = tmod(1)
!
! RECUPERATION DE LA NUMEROTATION DES BASES
!
    if ((typba1(1:9) .eq. 'MODE_MECA') .or. (typba1(1:9) .eq. 'MODE_GENE')) then
        call dismoi('REF_RIGI_PREM', base1, 'RESU_DYNA', repk=matri1, arret='C')
        call exisd('MATR_ASSE', matri1, iret)
        if (iret .ne. 0) then
            call dismoi('NOM_NUME_DDL', matri1, 'MATR_ASSE', repk=numdd1)
        else
            call dismoi('NUME_DDL', base1, 'RESU_DYNA', repk=numdd1)
        end if
    else
!       On passe par la numerotation du REFD
        call dismoi('NUME_DDL', base1, 'RESU_DYNA', repk=numdd1)
    end if
    call exisd('NUME_DDL', numdd1, iret)
    if (iret .ne. 1) then
        call utmess('F', 'CALCESSAI0_14', sk=base1)
    end if
!
    call jeveuo(numdd1//'.NUME.NEQU', 'L', vi=nllneq1)
    neq1 = nllneq1(1)
!
!
    if ((typba2(1:9) .eq. 'MODE_MECA') .or. (typba2(1:9) .eq. 'MODE_GENE')) then
        call dismoi('REF_RIGI_PREM', base2, 'RESU_DYNA', repk=matri2, arret='C')
        call exisd('MATR_ASSE', matri2, iret)
        if (iret .ne. 0) then
            call dismoi('NOM_NUME_DDL', matri2, 'MATR_ASSE', repk=numdd2)
        else
            call dismoi('NUME_DDL', base2, 'RESU_DYNA', repk=numdd2)
        end if
    else
        call dismoi('NUME_DDL', base2, 'RESU_DYNA', repk=numdd2)
    end if
    call exisd('NUME_DDL', numdd2, iret)
    if (iret .ne. 1) then
        call utmess('F', 'CALCESSAI0_14', sk=base2)
    end if
    call jeveuo(numdd2//'.NUME.NEQU', 'L', llneq2)
!
! ---- Verification : les deux nume_ddl doivent etre identiques
    pronu1 = (numdd1//'.NUME')
    pronu2 = (numdd2//'.NUME')
    if (.not. idensd('NUME_EQUA', pronu1, pronu2)) then
        call utmess('F', 'ALGELINE2_80')
    end if
!
! --- Verification : le nume_ddl doit etre celui de la MATR_ASSE
    if (matr .ne. ' ') then
        pronua = (numdda//'.NUME')
        if (.not. idensd('NUME_EQUA', pronu1, pronua)) then
            call utmess('F', 'ALGELINE2_81')
        end if
        nu = numdda(1:14)
        call jeveuo(nu//'.NUME.DEEQ', 'L', vi=deeq)
    else
        nu = numdd1(1:14)
        neq = neq1
    end if
!
! INITIALISATION DE LA TABLE DES MACS
    if (zcmplx) then
        nbpara = 3
    else
        nbpara = 4
    end if
    call wkvect('&&OP0141.TYPE_PARA', 'V V K8 ', nbpara, ityp)
    call wkvect('&&OP0141.NOM_PARA', 'V V K16', nbpara, inom)
    do i = 1, 2
        zk8(ityp+i-1) = 'I'
    end do
    if (zcmplx) then
        call wkvect('&&OP0141.MAC', 'V V R', 1, indv)
    else
        call wkvect('&&OP0141.MAC', 'V V R', 2, indv)
! MATRICE GENERALISEE EN PLUS POUR LES MODES REELS
        zk16(inom+3) = 'Y1_W_Y2'
        zk8(ityp+3) = 'R'
    end if
    zk8(ityp+2) = 'R'
    zk16(inom) = 'NUME_MODE_1'
    zk16(inom+1) = 'NUME_MODE_2'
    if (ieri) then
        zk16(inom+2) = 'IERI'
    else
        zk16(inom+2) = 'MAC'
    end if
    call tbajpa(table, nbpara, zk16(inom), zk8(ityp))
!
    call wkvect('&&OP0141.IJ', 'V V I', 2, ind)
!
    if (zcmplx) then
!
        call wkvect('&&OP0141.BASE1', 'V V C', nbmod1*neq, idbas1)
        call wkvect('&&OP0141.BASE2', 'V V C', nbmod2*neq, idbas2)
        call wkvect('&&OP0141.BASE3', 'V V C', neq, idbas3)
        call wkvect('&&OP0141.TEMP1', 'V V C', neq, idvec1)
        call wkvect('&&OP0141.TEMP2', 'V V C', neq, idvec2)
        call wkvect('&&OP0141.TEMP3', 'V V C', neq, idvec3)
!
        call copmod(base1, numer=nu, bmodz=zc(idbas1))
        call copmod(base2, numer=nu, bmodz=zc(idbas2))
!
! BOUCLE DE CALCUL DES MACS
        do i = 1, nbmod1
            pii = 0.d0
            if (matr .ne. ' ') then
                call mcmult('ZERO', imatra, zc(idbas1+(i-1)*neq), zc(idvec1), 1, &
                            .true._1)
!
                do iddl = 1, neq
                    if (deeq(2*iddl) .le. 0) zc(idvec1-1+iddl) = dcmplx(0.d0, 0.d0)
                end do
!
            else
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call zcopy(b_n, zc(idbas1+(i-1)*neq), b_incx, zc(idvec1), b_incy)
            end if
!
! PB AVEC ZDOTC DE BLAS POUR CERTAIN COMPILO -> CALCUL DIRECT
            ztemp = dcmplx(0.0d0, 0.0d0)
            do iddl = 1, neq
                ztemp = ztemp+zc(idbas1+(i-1)*neq-1+iddl)*dconjg(zc(idvec1-1+iddl))
            end do
            pii = abs(ztemp)
!
            zi(ind) = i
!
            do j = 1, nbmod2
                pij = 0.d0
                pjj = 0.d0
                if (matr .ne. ' ') then
                    call mcmult('ZERO', imatra, zc(idbas2+(j-1)*neq), zc(idvec2), 1, &
                                .true._1)
!
                    do iddl = 1, neq
                        if (deeq(2*iddl) .le. 0) zc(idvec2-1+iddl) = dcmplx(0.d0, 0.d0)
                    end do
!
                else
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call zcopy(b_n, zc(idbas2+(j-1)*neq), b_incx, zc(idvec2), b_incy)
                end if
!
                ztemp = dcmplx(0.0d0, 0.0d0)
                do iddl = 1, neq
                    ztemp = ztemp+zc(idbas2+(j-1)*neq-1+iddl)*dconjg(zc(idvec2-1+iddl))
                end do
                pjj = abs(ztemp)
!
                if (ieri) then
                    do iddl = 1, neq
                        zc(idbas3-1+iddl) = zc( &
                                            idbas1+(i-1)*neq-1+iddl)-zc(idbas2+(j-1)*neq-1+iddl)
                    end do
                    call mcmult('ZERO', imatra, zc(idbas3), zc(idvec3), 1, &
                                .true._1)
                    do iddl = 1, neq
                        if (deeq(2*iddl) .le. 0) zc(idvec3-1+iddl) = dcmplx(0.d0, 0.d0)
                    end do
!
                    ztemp = dcmplx(0.0d0, 0.0d0)
                    do iddl = 1, neq
                        ztemp = ztemp+zc(idbas3-1+iddl)*dconjg(zc(idvec3-1+iddl))
                    end do
                    pij = abs(ztemp)
!
                    pij = (pij**2)/(pii**2+pjj**2)
                else
                    ztemp = dcmplx(0.0d0, 0.0d0)
                    do iddl = 1, neq
                        ztemp = ztemp+zc(idbas1+(i-1)*neq-1+iddl)*dconjg(zc(idvec2-1+iddl))
                    end do
                    pij = abs(ztemp)
                    pij = (pij**2)/(pii*pjj)
                end if
!
                zi(ind+1) = j
                zr(indv) = pij
                call tbajli(table, nbpara, zk16(inom), zi(ind), zr(indv), &
                            [cbid], k8b, 0)
            end do
        end do
!
    else
!
        call wkvect('&&OP0141.BASE1', 'V V R', nbmod1*neq, idbas1)
        call wkvect('&&OP0141.BASE2', 'V V R', nbmod2*neq, idbas2)
        call wkvect('&&OP0141.BASE3', 'V V R', neq, idbas3)
        call wkvect('&&OP0141.TEMP1', 'V V R', neq, idvec1)
        call wkvect('&&OP0141.TEMP2', 'V V R', neq, idvec2)
        call wkvect('&&OP0141.TEMP3', 'V V R', neq, idvec3)
!
        call copmod(base1, numer=nu, bmodr=zr(idbas1))
        call copmod(base2, numer=nu, bmodr=zr(idbas2))
!
! BOUCLE DE CALCUL DES MACS
        do i = 1, nbmod1
            pii = 0.d0
            if (matr .ne. ' ') then
                call mrmult('ZERO', imatra, zr(idbas1+(i-1)*neq), zr(idvec1), 1, &
                            .true._1)
                call zerlag(neq, deeq, vectr=zr(idvec1))
            else
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(idbas1+(i-1)*neq), b_incx, zr(idvec1), b_incy)
            end if
!
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            pii = abs(ddot(b_n, zr(idbas1+(i-1)*neq), b_incx, zr(idvec1), b_incy))
!
            zi(ind) = i
!
            do j = 1, nbmod2
                pij = 0.d0
                pjj = 0.d0
                if (matr .ne. ' ') then
                    call mrmult('ZERO', imatra, zr(idbas2+(j-1)*neq), zr(idvec2), 1, &
                                .true._1)
                    call zerlag(neq, deeq, vectr=zr(idvec2))
                else
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(idbas2+(j-1)*neq), b_incx, zr(idvec2), b_incy)
                end if
!
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                pjj = abs(ddot(b_n, zr(idbas2+(j-1)*neq), b_incx, zr(idvec2), b_incy))
!
                if (ieri) then
                    do ieq = 1, neq
                        zr(idbas3-1+ieq) = zr(idbas1+neq*(i-1)-1+ieq)-zr(idbas2+neq*(i-1)-1+ieq)
                    end do
                    call mrmult('ZERO', imatra, zr(idbas3), zr(idvec3), 1, &
                                ASTER_TRUE)
                    call zerlag(neq, deeq, vectr=zr(idvec3))
!
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    pij = abs(ddot(b_n, zr(idbas3), b_incx, zr(idvec3), b_incy))
!
                    pij = (pij**2)/(pii**2+pjj**2)
!  POUR LA MATRICE GENERALISEE : Y1_W_Y2
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    rbid = abs(ddot(b_n, zr(idbas1+(i-1)*neq), b_incx, zr(idvec2), b_incy))
!
                    rbid = (rbid**2)/(pii*pjj)
                    zr(indv+1) = sqrt(rbid*pii*pjj)
                else
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    pij = abs(ddot(b_n, zr(idbas1+(i-1)*neq), b_incx, zr(idvec2), b_incy))
!
                    pij = (pij**2)/(pii*pjj)
                    zr(indv+1) = sqrt(pij*pii*pjj)
                end if
!
                zi(ind+1) = j
                zr(indv) = pij
                call tbajli(table, nbpara, zk16(inom), zi(ind), zr(indv), &
                            [cbid], k8b, 0)
            end do
        end do
!
    end if
!  FIN TEST SUR TYPE DE VECTEURS (C/R)
!
    if (niv .ge. 2) then
        call tbimpr(table, 'TABLEAU', ifm, 3, zk16(inom), &
                    0, '1PE12.5')
        if (nbpara .eq. 4) then
            write (ifm, *) ' '
            write (ifm, 1000) zk16(inom+2)
            call tbimex(table, ifm, 4, zk16(inom), 'EXCEL', &
                        '1PE12.5')
            write (ifm, *) ' '
        end if
    end if
1000 format('AFFICHAGE ', a4, ' ET MATRICE GENERALISEE : ')
!
    call jedema()
end subroutine
