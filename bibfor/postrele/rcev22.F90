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
subroutine rcev22(nbinti, kinti, iocc, csigm, cinst, &
                  ccont, lfatig, flexio, lrocht, cnoc, &
                  cresu, cpres, lsymm)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rc32axis.h"
#include "asterfort/rc32my.h"
#include "asterfort/tbexip.h"
#include "asterfort/tbextb.h"
#include "asterfort/tbexv1.h"
#include "asterfort/tbliva.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbinti
    aster_logical :: lfatig, flexio, lrocht, lsymm
    character(len=16) :: kinti
    character(len=24) :: csigm, cinst, ccont, cnoc, cresu, cpres
!     OPERATEUR POST_RCCM, TYPE_RESU_MECA='EVOLUTION'
!     LECTURE DU MOT CLE FACTEUR "TRANSITOIRE"
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ibid, n1, nbinst, kinst, ncmpr, i, j, k, l, l1, l2, l3, ndim
    integer(kind=8) :: nbabsc, jabsc, nbxcoo, jxcoo, nbycoo, jycoo
    integer(kind=8) :: jsigm, jinst, ncmp, iret, jsioe, iocc, nbins0
    integer(kind=8) :: jnocc, jresu, nbcycl, jresp
    parameter(ncmp=6)
    real(kind=8) :: r8b, prec(2), momen0, momen1, vale(2)
    real(kind=8) :: momen0_axis(4), momen1_axis(4), momen2_axis(4), rho
    complex(kind=8) :: cbid
    aster_logical :: exist, cfait
    character(len=8) :: k8b, crit(2), nocmp(ncmp)
    character(len=16) :: motclf, valek(2), table, tabl0, tabfle, tabfl0, tabpre
    character(len=16) :: tabpr0
    character(len=19) :: nomf
    character(len=24) :: instan, abscur, coord_x, coord_y
    character(len=24) :: valk(7)
    real(kind=8), allocatable :: cont_flexio(:, :), cont_pressi(:, :), contraintes(:, :)
!    real(kind=8), pointer :: cont_flexio(:) => null()
!    real(kind=8), pointer :: cont_pressi(:) => null()
!    real(kind=8), pointer :: contraintes(:) => null()
! DEB ------------------------------------------------------------------
    call jemarq()
!
    motclf = 'TRANSITOIRE'
!
    if (lsymm) then
        call getvr8(' ', 'RAYON_MOYEN', scal=rho, nbret=n1)
    end if
!
    nocmp(1) = 'SIXX'
    nocmp(2) = 'SIYY'
    nocmp(3) = 'SIZZ'
    nocmp(4) = 'SIXY'
    nocmp(5) = 'SIXZ'
    nocmp(6) = 'SIYZ'
!
    valek(1) = 'INST            '
    valek(2) = 'ABSC_CURV       '
!
    prec(1) = 1.0d-06
    prec(2) = 1.0d-06
    crit(1) = 'RELATIF'
    crit(2) = 'RELATIF'
!
    instan = '&&RCEV22.INSTANT'
    abscur = '&&RCEV22.ABSC_CURV'
    coord_x = '&&RCEVO2.COOR_X'
    coord_y = '&&RCEVO2.COOR_Y'
!
! --- RECHERCHE DU NOMBRE D'INSTANTS A COMBINER
!
    nbinst = 0
    cfait = .false.
!
    call getvid(motclf, 'TABL_RESU_MECA', iocc=iocc, scal=tabl0, nbret=n1)
    call getvid(motclf, 'TABL_SIGM_THER', iocc=iocc, scal=tabfl0, nbret=n1)
    if (n1 .ne. 0) flexio = .true.
    call getvid(motclf, 'TABL_RESU_PRES', iocc=iocc, scal=tabpr0, nbret=n1)
    if (n1 .ne. 0) then
        lrocht = .true.
        call wkvect(cpres, 'V V K8', 1, jresp)
        zk8(jresp-1+1) = tabpr0
    end if
!
    call getvr8(motclf, 'INST', iocc=iocc, scal=r8b, nbret=n1)
    if (n1 .ne. 0) then
        nbins0 = -n1
    else
        call getvid(motclf, 'LIST_INST', iocc=iocc, scal=nomf, nbret=n1)
        if (n1 .ne. 0) then
            call jelira(nomf//'.VALE', 'LONMAX', nbins0)
        else
            if (nbinti .eq. 1) then
                table = tabl0
                tabfle = tabfl0
                tabpre = tabpr0
            else
                cfait = .true.
                table = '&&RCEV22.RESU_ME'
                tabfle = '&&RCEV22.SIGM_TH'
                tabpre = '&&RCEV22.RESU_PR'
                call tbextb(tabl0, 'V', table, 1, 'INTITULE', &
                            'EQ', [ibid], [r8b], [cbid], kinti, &
                            [r8b], k8b, iret)
                if (iret .eq. 10) then
                    valk(1) = 'INTITULE'
                    valk(2) = tabl0
                    call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
                else if (iret .eq. 20) then
                    valk(1) = tabl0
                    valk(2) = 'INTITULE'
                    call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
                end if
                if (flexio) then
                    call tbextb(tabfl0, 'V', tabfle, 1, 'INTITULE', &
                                'EQ', [ibid], [r8b], [cbid], kinti, &
                                [r8b], k8b, iret)
                    if (iret .eq. 10) then
                        valk(1) = 'INTITULE'
                        valk(2) = tabfl0
                        call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
                    else if (iret .eq. 20) then
                        valk(1) = tabfl0
                        valk(2) = 'INTITULE'
                        call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
                    end if
                end if
                if (lrocht) then
                    call tbextb(tabpr0, 'V', tabpre, 1, 'INTITULE', &
                                'EQ', [ibid], [r8b], [cbid], kinti, &
                                [r8b], k8b, iret)
                    if (iret .eq. 10) then
                        valk(1) = 'INTITULE'
                        valk(2) = tabpr0
                        call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
                    else if (iret .eq. 20) then
                        valk(1) = tabpr0
                        valk(2) = 'INTITULE'
                        call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
                    end if
                end if
            end if
            call tbexip(table, valek(1), exist, k8b)
            if (.not. exist) then
                valk(1) = table
                valk(2) = 'INTITULE'
                valk(3) = kinti
                valk(4) = valek(1)
                call utmess('F', 'POSTRCCM_17', nk=4, valk=valk)
            end if
            call tbexv1(table, valek(1), instan, 'V', nbins0, &
                        k8b)
            call jedetr(instan)
        end if
    end if
    nbinst = nbinst+nbins0
!
! --- PRESENCE DES COMPOSANTES DANS LA TABLE
!
    if (.not. cfait) then
        if (nbinti .eq. 1) then
            table = tabl0
            tabfle = tabfl0
            tabpre = tabpr0
        else
            cfait = .true.
            table = '&&RCEV22.RESU_ME'
            tabfle = '&&RCEV22.SIGM_TH'
            tabpre = '&&RCEV22.RESU_PR'
            call tbextb(tabl0, 'V', table, 1, 'INTITULE', &
                        'EQ', [ibid], [r8b], [cbid], kinti, &
                        [r8b], k8b, iret)
            if (iret .eq. 10) then
                valk(1) = 'INTITULE'
                valk(2) = tabl0
                call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
            else if (iret .eq. 20) then
                valk(1) = tabl0
                valk(2) = 'INTITULE'
                call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
            end if
            if (flexio) then
                call tbextb(tabfl0, 'V', tabfle, 1, 'INTITULE', &
                            'EQ', [ibid], [r8b], [cbid], kinti, &
                            [r8b], k8b, iret)
                if (iret .eq. 10) then
                    valk(1) = 'INTITULE'
                    valk(2) = tabfl0
                    call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
                else if (iret .eq. 20) then
                    valk(1) = tabfl0
                    valk(2) = 'INTITULE'
                    call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
                end if
            end if
            if (lrocht) then
                call tbextb(tabpr0, 'V', tabpre, 1, 'INTITULE', &
                            'EQ', [ibid], [r8b], [cbid], kinti, &
                            [r8b], k8b, iret)
                if (iret .eq. 10) then
                    valk(1) = 'INTITULE'
                    valk(2) = tabpr0
                    call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
                else if (iret .eq. 20) then
                    valk(1) = tabpr0
                    valk(2) = 'INTITULE'
                    call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
                end if
            end if
        end if
    end if
!
    ncmpr = 6
    do i = 1, 4
        call tbexip(table, nocmp(i), exist, k8b)
        if (.not. exist) then
            valk(1) = table
            valk(2) = 'INTITULE'
            valk(3) = kinti
            valk(4) = nocmp(i)
            call utmess('F', 'POSTRCCM_17', nk=4, valk=valk)
        end if
        if (flexio) then
            call tbexip(tabfle, nocmp(i), exist, k8b)
            if (.not. exist) then
                valk(1) = tabfle
                valk(2) = 'INTITULE'
                valk(3) = kinti
                valk(4) = nocmp(i)
                call utmess('F', 'POSTRCCM_17', nk=4, valk=valk)
            end if
        end if
        if (lrocht) then
            call tbexip(tabpre, nocmp(i), exist, k8b)
            if (.not. exist) then
                valk(1) = tabpre
                valk(2) = 'INTITULE'
                valk(3) = kinti
                valk(4) = nocmp(i)
                call utmess('F', 'POSTRCCM_17', nk=4, valk=valk)
            end if
        end if
    end do
    call tbexip(table, nocmp(5), exist, k8b)
    if (.not. exist) ncmpr = 4
!
! --- ON RECUPERE L'ABSC_CURV DANS LA TABLE
!
    call tbexip(table, valek(2), exist, k8b)
    if (.not. exist) then
        valk(1) = table
        valk(2) = 'INTITULE'
        valk(3) = kinti
        valk(4) = valek(2)
        call utmess('F', 'POSTRCCM_17', nk=4, valk=valk)
    end if
    call tbexv1(table, valek(2), abscur, 'V', nbabsc, &
                k8b)
    call tbexv1(table, 'COOR_X          ', coord_x, 'V', nbxcoo, &
                k8b)
    call tbexv1(table, 'COOR_Y          ', coord_y, 'V', nbycoo, &
                k8b)
!
    call jeveuo(abscur, 'L', jabsc)
    call jeveuo(coord_x, 'L', jxcoo)
    call jeveuo(coord_y, 'L', jycoo)
! --- AU CAS OU LA LIGNE COUPE EST HORIZONTALE OU VERTICALE
    if (nbxcoo .eq. 1) then
        do i = 1, nbabsc
            zr(jxcoo-1+i) = zr(jxcoo)
        end do
    end if
    if (nbycoo .eq. 1) then
        do i = 1, nbabsc
            zr(jycoo-1+i) = zr(jycoo)
        end do
    end if
    allocate (contraintes(ncmpr, nbabsc))
    allocate (cont_flexio(ncmpr, nbabsc))
    allocate (cont_pressi(ncmpr, nbabsc))
!    AS_ALLOCATE(vr=contraintes, size=nbabsc)
!    AS_ALLOCATE(vr=cont_flexio, size=nbabsc)
!    AS_ALLOCATE(vr=cont_pressi, size=nbabsc)
!
! --- CREATION DES OBJETS DE TRAVAIL
!
    if (lsymm) then
        ndim = 9*nbinst*ncmp
    else
        ndim = 6*nbinst*ncmp
    end if
    call wkvect(csigm, 'V V R', ndim, jsigm)
    call wkvect(cinst, 'V V R', nbinst, jinst)
    call wkvect(cnoc, 'V V I', nbinst, jnocc)
    call wkvect(cresu, 'V V K8', nbinst, jresu)
    if (lfatig) then
        ndim = 2*nbinst*ncmp
        call wkvect(ccont, 'V V R', ndim, jsioe)
    end if
!
! --- RECUPERATION DES INFORMATIONS
!
    call getvis(motclf, 'NB_OCCUR', iocc=iocc, scal=nbcycl, nbret=n1)
!
! --- ON RECUPERE LES INSTANTS DANS LA TABLE
!
    call getvr8(motclf, 'INST', iocc=iocc, scal=r8b, nbret=n1)
    if (n1 .ne. 0) then
        nbins0 = -n1
        call wkvect(instan, 'V V R', nbins0, kinst)
        call getvr8(motclf, 'INST', iocc=iocc, nbval=nbins0, vect=zr(kinst), &
                    nbret=n1)
        call getvr8(motclf, 'PRECISION', iocc=iocc, scal=prec(1), nbret=n1)
        call getvtx(motclf, 'CRITERE', iocc=iocc, scal=crit(1), nbret=n1)
    else
        call getvid(motclf, 'LIST_INST', iocc=iocc, scal=nomf, nbret=n1)
        if (n1 .ne. 0) then
            call jelira(nomf//'.VALE', 'LONMAX', nbins0)
            call jeveuo(nomf//'.VALE', 'L', kinst)
            call getvr8(motclf, 'PRECISION', iocc=iocc, scal=prec(1), nbret=n1)
            call getvtx(motclf, 'CRITERE', iocc=iocc, scal=crit(1), nbret=n1)
        else
            prec(1) = 1.0d-06
            crit(1) = 'RELATIF'
            call tbexip(table, valek(1), exist, k8b)
            if (.not. exist) then
                valk(1) = table
                valk(2) = 'INTITULE'
                valk(3) = kinti
                valk(4) = valek(1)
                call utmess('F', 'POSTRCCM_17', nk=4, valk=valk)
            end if
            call tbexv1(table, valek(1), instan, 'V', nbins0, &
                        k8b)
            call jeveuo(instan, 'L', kinst)
        end if
    end if
!
! ---
    do i = 1, nbins0
!
        zr(jinst+i-1) = zr(kinst+i-1)
        zi(jnocc-1+i) = nbcycl
        zk8(jresu-1+i) = tabl0
!
        vale(1) = zr(kinst+i-1)
!
        do j = 1, ncmpr
!
            do k = 1, nbabsc
                vale(2) = zr(jabsc+k-1)
!
                call tbliva(table, 2, valek, [ibid], vale, &
                            [cbid], k8b, crit, prec, nocmp(j), &
                            k8b, ibid, contraintes(j, k), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = table
                    valk(2) = nocmp(j)
                    valk(3) = valek(1)
                    valk(4) = valek(2)
                    call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                valr=vale)
                end if
!
                if (flexio) then
                    call tbliva(tabfle, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, cont_flexio(j, k), cbid, k8b, &
                                iret)
                    if (iret .ne. 0) then
                        valk(1) = tabfle
                        valk(2) = nocmp(j)
                        valk(3) = valek(1)
                        valk(4) = valek(2)
                        call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                    valr=vale)
                    end if
                end if
!
                if (lrocht) then
                    call tbliva(tabpre, 2, valek, [ibid], vale, &
                                [cbid], k8b, crit, prec, nocmp(j), &
                                k8b, ibid, cont_pressi(j, k), cbid, k8b, &
                                iret)
                    if (iret .ne. 0) then
                        valk(1) = tabpre
                        valk(2) = nocmp(j)
                        valk(3) = valek(1)
                        valk(4) = valek(2)
                        call utmess('F', 'POSTRCCM_2', nk=4, valk=valk, nr=2, &
                                    valr=vale)
                    end if
                end if
!
            end do
!
            if (lfatig) then
                l = ncmp*(i-1)+j
                zr(jsioe-1+l) = contraintes(j, 1)
                l = ncmp*nbinst+ncmp*(i-1)+j
                zr(jsioe-1+l) = contraintes(j, nbabsc)
            end if
!
            if (.not. lsymm) then
                call rc32my(nbabsc, zr(jabsc), contraintes(j, :), momen0, momen1)
                momen1 = 0.5d0*momen1
                !
                l = ncmp*(i-1)+j
                zr(jsigm-1+l) = momen0
                l = ncmp*nbinst+ncmp*(i-1)+j
                zr(jsigm-1+l) = momen1
                !
                if (flexio) then
                    call rc32my(nbabsc, zr(jabsc), cont_flexio(j, :), momen0, momen1)
                    momen1 = 0.5d0*momen1
                else
                    momen0 = 0.d0
                    momen1 = 0.d0
                end if
                l = 2*ncmp*nbinst+ncmp*(i-1)+j
                zr(jsigm-1+l) = momen0
                l = 3*ncmp*nbinst+ncmp*(i-1)+j
                zr(jsigm-1+l) = momen1
                !
                if (lrocht) then
                    call rc32my(nbabsc, zr(jabsc), cont_pressi(j, :), momen0, momen1)
                    momen1 = 0.5d0*momen1
                else
                    momen0 = r8vide()
                    momen1 = r8vide()
                end if
                l = 4*ncmp*nbinst+ncmp*(i-1)+j
                zr(jsigm-1+l) = momen0
                l = 5*ncmp*nbinst+ncmp*(i-1)+j
                zr(jsigm-1+l) = momen1
            end if
!
        end do
        if (lsymm) then
            call rc32axis(nbabsc, zr(jabsc), zr(jxcoo), zr(jycoo), contraintes, &
                          momen0_axis, momen1_axis, momen2_axis, rho)
            l1 = ncmp*(i-1)
            l2 = ncmp*nbinst+ncmp*(i-1)
            l3 = 2*ncmp*nbinst+ncmp*(i-1)
            do j = 1, 4
                zr(jsigm-1+l1+j) = momen0_axis(j)
                zr(jsigm-1+l2+j) = momen1_axis(j)
                zr(jsigm-1+l3+j) = momen2_axis(j)
            end do
!
            if (flexio) then
                call rc32axis(nbabsc, zr(jabsc), zr(jxcoo), zr(jycoo), cont_flexio, &
                              momen0_axis, momen1_axis, momen2_axis, rho)
            else
                momen0_axis = 0.d0
                momen1_axis = 0.d0
                momen2_axis = 0.d0
            end if
            l1 = 3*ncmp*nbinst+ncmp*(i-1)
            l2 = 4*ncmp*nbinst+ncmp*(i-1)
            l3 = 5*ncmp*nbinst+ncmp*(i-1)
            do j = 1, 4
                zr(jsigm-1+l1+j) = momen0_axis(j)
                zr(jsigm-1+l2+j) = momen1_axis(j)
                zr(jsigm-1+l3+j) = momen2_axis(j)
            end do
!
            if (lrocht) then
                call rc32axis(nbabsc, zr(jabsc), zr(jxcoo), zr(jycoo), cont_pressi, &
                              momen0_axis, momen1_axis, momen2_axis, rho)
            else
                momen0_axis = 0.d0
                momen1_axis = 0.d0
                momen2_axis = 0.d0
            end if
            l1 = 6*ncmp*nbinst+ncmp*(i-1)
            l2 = 7*ncmp*nbinst+ncmp*(i-1)
            l3 = 8*ncmp*nbinst+ncmp*(i-1)
            do j = 1, 4
                zr(jsigm-1+l1+j) = momen0_axis(j)
                zr(jsigm-1+l2+j) = momen1_axis(j)
                zr(jsigm-1+l3+j) = momen2_axis(j)
            end do
        end if
!
    end do
!
    call jedetr(instan)
    call jedetr(abscur)
    call jedetr(coord_x)
    call jedetr(coord_y)
    deallocate (contraintes)
    deallocate (cont_flexio)
    deallocate (cont_pressi)
!    AS_DEALLOCATE(vr=contraintes)
!    AS_DEALLOCATE(vr=cont_flexio)
!    AS_DEALLOCATE(vr=cont_pressi)
    if (nbinti .ne. 1) then
        call detrsd('TABLE', table)
        if (flexio) call detrsd('TABLE', tabfle)
        if (lrocht) call detrsd('TABLE', tabpre)
    end if
!
    call jedema()
end subroutine
