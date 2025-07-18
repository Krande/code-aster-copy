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
subroutine rc32mu()
    implicit none
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3200
!     LECTURE DU MOT CLE FACTEUR "RESU_MECA_UNIT"
!
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rc32my.h"
#include "asterfort/rcver1.h"
#include "asterfort/rcveri.h"
#include "asterfort/tbexip.h"
#include "asterfort/tbexv1.h"
#include "asterfort/tbliva.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ibid, ns(13), nbabsc, jabsc, iret, jmune, jmuno, i, j, k, l, ndim
    integer(kind=8) :: ncmp, nb, kk
    parameter(ncmp=6)
    real(kind=8) :: prec, momen0, momen1
    complex(kind=8) :: cbid
    aster_logical :: exist
    character(len=8) :: k8b, crit, nocmp(ncmp), tbsig(13)
    character(len=16) :: motclf, valek
    character(len=24) :: abscur
    character(len=24) :: valk(7)
    real(kind=8), pointer :: contraintes(:) => null()
! DEB ------------------------------------------------------------------
    call jemarq()
!
    motclf = 'RESU_MECA_UNIT'
    call getfac(motclf, nb)
!
    if (nb .eq. 0) goto 999
!
    prec = 1.0d-06
    crit = 'RELATIF'
!
    call getvid(motclf, 'TABL_MX', iocc=1, scal=tbsig(4), nbret=ns(4))
    if (ns(4) .eq. 0) then
        call getvid(motclf, 'TABL_MX_TUBU', iocc=1, scal=tbsig(4), nbret=ns(4))
    end if
    call rcveri(tbsig(4))
!
    call getvid(motclf, 'TABL_FX', iocc=1, scal=tbsig(1), nbret=ns(1))
    if (ns(1) .eq. 0) then
        call getvid(motclf, 'TABL_FX_TUBU', iocc=1, scal=tbsig(1), nbret=ns(1))
    end if
    if (ns(1) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(1))
!
    call getvid(motclf, 'TABL_FY', iocc=1, scal=tbsig(2), nbret=ns(2))
    if (ns(2) .eq. 0) then
        call getvid(motclf, 'TABL_FY_TUBU', iocc=1, scal=tbsig(2), nbret=ns(2))
    end if
    if (ns(2) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(2))
!
    call getvid(motclf, 'TABL_FZ', iocc=1, scal=tbsig(3), nbret=ns(3))
    if (ns(3) .eq. 0) then
        call getvid(motclf, 'TABL_FZ_TUBU', iocc=1, scal=tbsig(3), nbret=ns(3))
    end if
    if (ns(3) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(3))
!
    call getvid(motclf, 'TABL_MY', iocc=1, scal=tbsig(5), nbret=ns(5))
    if (ns(5) .eq. 0) then
        call getvid(motclf, 'TABL_MY_TUBU', iocc=1, scal=tbsig(5), nbret=ns(5))
    end if
    if (ns(5) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(5))
!
    call getvid(motclf, 'TABL_MZ', iocc=1, scal=tbsig(6), nbret=ns(6))
    if (ns(6) .eq. 0) then
        call getvid(motclf, 'TABL_MZ_TUBU', iocc=1, scal=tbsig(6), nbret=ns(6))
    end if
    if (ns(6) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(6))
!
    call getvid(motclf, 'TABL_FX_CORP', iocc=1, scal=tbsig(7), nbret=ns(7))
    if (ns(7) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(7))
!
    call getvid(motclf, 'TABL_FY_CORP', iocc=1, scal=tbsig(8), nbret=ns(8))
    if (ns(8) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(8))
!
    call getvid(motclf, 'TABL_FZ_CORP', iocc=1, scal=tbsig(9), nbret=ns(9))
    if (ns(9) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(9))
!
    call getvid(motclf, 'TABL_MX_CORP', iocc=1, scal=tbsig(10), nbret=ns(10))
    if (ns(10) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(10))
!
    call getvid(motclf, 'TABL_MY_CORP', iocc=1, scal=tbsig(11), nbret=ns(11))
    if (ns(11) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(11))
!
    call getvid(motclf, 'TABL_MZ_CORP', iocc=1, scal=tbsig(12), nbret=ns(12))
    if (ns(12) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(12))
!
    call getvid(motclf, 'TABL_PRES', iocc=1, scal=tbsig(13), nbret=ns(13))
    if (ns(13) .ne. 0) call rcver1('MECANIQUE', tbsig(4), tbsig(13))
!
! --- ON RECUPERE L'ABSC_CURV DANS LA TABLE 'TABL_MX'
!
    valek = 'ABSC_CURV       '
    call tbexip(tbsig(4), valek, exist, k8b)
    if (.not. exist) then
        valk(1) = tbsig(4)
        valk(2) = valek
        call utmess('F', 'POSTRCCM_1', nk=2, valk=valk)
    end if
    abscur = '&&RC32MU.ABSC_CURV'
    call tbexv1(tbsig(4), valek, abscur, 'V', nbabsc, &
                k8b)
    call jeveuo(abscur, 'L', jabsc)
!
    AS_ALLOCATE(vr=contraintes, size=nbabsc)
!
    nocmp(1) = 'SIXX'
    nocmp(2) = 'SIYY'
    nocmp(3) = 'SIZZ'
    nocmp(4) = 'SIXY'
    nocmp(5) = 'SIXZ'
    nocmp(6) = 'SIYZ'
!
! --- 13 TABLES A  ( 6 COMPOSANTES + 6 LINEARISEES + 6 M_0 + 6 M_1 )
    ndim = 13*(6+6+6+6)
    call wkvect('&&RC3200.MECA_UNIT .ORIG', 'V V R', ndim, jmuno)
    call wkvect('&&RC3200.MECA_UNIT .EXTR', 'V V R', ndim, jmune)
!
! --- LES PROFILS DE CONTRAINTES ISSUS DES CALCULS MECANIQUES UNITAIRES
!
    do kk = 1, ndim
        zr(jmuno+kk-1) = 0.d0
        zr(jmune+kk-1) = 0.d0
    end do
!
    do i = 1, 13
!
        if (ns(i) .eq. 0) goto 10
!
        call tbexip(tbsig(i), valek, exist, k8b)
        if (.not. exist) then
            valk(1) = tbsig(i)
            valk(2) = valek
            call utmess('F', 'POSTRCCM_1', nk=2, valk=valk)
        end if
        do j = 1, ncmp
!
            do k = 1, nbabsc
                call tbliva(tbsig(i), 1, valek, [ibid], zr(jabsc+k-1), &
                            [cbid], k8b, crit, [prec], nocmp(j), &
                            k8b, ibid, contraintes(k), cbid, k8b, &
                            iret)
                if (iret .ne. 0) then
                    valk(1) = tbsig(i)
                    valk(2) = nocmp(j)
                    valk(3) = valek
                    call utmess('F', 'POSTRCCM_44', nk=3, valk=valk, sr=zr(jabsc+k-1))
                end if
            end do
!
            l = ncmp*(i-1)+j
            zr(jmuno-1+l) = contraintes(1)
            zr(jmune-1+l) = contraintes(nbabsc)
!
            call rc32my(nbabsc, zr(jabsc), contraintes, momen0, momen1)
!
            l = 13*ncmp+ncmp*(i-1)+j
            zr(jmuno-1+l) = momen0-0.5d0*momen1
            zr(jmune-1+l) = momen0+0.5d0*momen1
!
            l = 2*13*ncmp+ncmp*(i-1)+j
            zr(jmuno-1+l) = momen0
            zr(jmune-1+l) = momen0
!
            l = 3*13*ncmp+ncmp*(i-1)+j
            zr(jmuno-1+l) = 0.5d0*momen1
            zr(jmune-1+l) = 0.5d0*momen1
!
        end do
!
10      continue
    end do
!
    call jedetr(abscur)
    AS_DEALLOCATE(vr=contraintes)
!
999 continue
    call jedema()
end subroutine
