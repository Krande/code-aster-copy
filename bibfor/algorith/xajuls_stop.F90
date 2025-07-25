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
subroutine xajuls_stop(noma, cnslt, jconx1, jconx2, ima)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/conare.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "blas/ddot.h"
!
    character(len=8), intent(in) :: noma
    character(len=19), intent(in) :: cnslt
    integer(kind=8), intent(in) :: jconx1
    integer(kind=8), intent(in) :: jconx2
    integer(kind=8), intent(in) :: ima
!
! person_in_charge: sam.cuvilliez at edf.fr
!
! ---------------------------------------------------------------------
!
!   XFEM : sous-routine de xajuls (reajustement des level set)
!
! ---------------------------------------------------------------------
!
!   But : Dans xajuls, les criteres de reajustement quadratiques
!         peuvent conduire a une division par zero. Dans ce cas, on
!         autorise le court-circuit d'iteration de la boucle sur les
!         aretes de la maille courante ima sous certaines conditions
!         tres particulieres :
!         1. tous les noeuds sommets de ima doivent verifier abs(lst)
!            < r8prem
!         2. les level-sets doivent avoir ete calculees depuis le
!            catalogue de formes geometriques de DEFI_FISS_XFEM dans
!            le cas 2D FORM_FISS='SEGMENT'
!         3. pour chaque arete de la maille courante ima, on compare
!            la longueur de son projete orthogonal sur le segment
!            (qui constitue la fissure) a la longueur de ce segment.
!            Au moins un de ces projetes doit avoir une longueur
!            comparable (critere relatif) a celle du segment. De
!            cette maniere on s'assure que la maille se trouve
!            "loin" de la fissure (a moins que le maillage ne soit
!            extrement grossier)
!
! ---------------------------------------------------------------------
!
!   in / noma   : nom du maillage
!   in / cnslt  : nom du cham_no_s de lst
!   in / jconx1 : adresse connectivite du maillage noma
!   in / jconx2 : adresse LONCUM(connectivite) du maillage noma
!   in / ima    : numero de la maille courante dans la boucle de xajuls
!
! ---------------------------------------------------------------------
!
    integer(kind=8) :: mxval, nbret, ibid, nbar, ia, nbnos, nunoa, nunob, ino
    integer(kind=8) :: i, nuno, na, nb
    integer(kind=8) :: ar(12, 3)
    real(kind=8) :: norm_s, norm_p, pscal, r8pre, diffe, crit
    real(kind=8) :: vect1(3), vect2(3), v_seg(3), v_are(3), v_pro(3)
    character(len=8) :: k8typm
    character(len=16) :: geofi
    aster_logical :: l_crit
    integer(kind=8), pointer :: vi_tym(:) => null()
    real(kind=8), pointer :: vr_lts(:) => null()
    real(kind=8), pointer :: vr_coo(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
!   rq : crit est le critere relatif en deca duquel on considere que
!   norm_p == norm_s (cf "verif 3" plus bas). Le choix de la valeur
!   de crit n'est base sur aucun argument geometrique et est discutable
    parameter(crit=1.d-1)
!
! ---------------------------------------------------------------------
!
    call jemarq()
!
! --
!   recuperation des objets necessaires
! --
!
    call jeveuo(cnslt//'.CNSV', 'L', vr=vr_lts)
    call jeveuo(noma//'.TYPMAIL', 'L', vi=vi_tym)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vr_coo)
!
! --
!   verif 1 : tous les noeuds sommets doivent avoir abs(lst) < r8prem
! --
!
    call jenuno(jexnum('&CATA.TM.NOMTM', vi_tym(ima)), k8typm)
    call conare(k8typm, ar, nbar)
!   nombre de noeuds sommets = nombre d'aretes
    nbnos = nbar
    r8pre = r8prem()
    do ino = 1, nbnos
        nuno = zi(jconx1-1+zi(jconx2+ima-1)+ino-1)
        ASSERT(abs(vr_lts(nuno)) .lt. r8pre)
    end do
!
! --
!   verif 2 : les level-sets on ete calculees depuis FORM_FISS == SEGMENT
!             (seul cas tolere)
! --
!
    mxval = 0
    nbret = 0
    geofi = ''
    call getvtx('DEFI_FISS', 'FORM_FISS', iocc=1, nbval=mxval, vect=geofi, &
                nbret=nbret)
    if (nbret .eq. -1) then
        call getvtx('DEFI_FISS', 'FORM_FISS', iocc=1, scal=geofi, nbret=nbret)
        if (geofi .eq. 'SEGMENT') then
            call getvr8('DEFI_FISS', 'PFON_ORIG', iocc=1, nbval=3, vect=vect1, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'PFON_EXTR', iocc=1, nbval=3, vect=vect2, &
                        nbret=ibid)
        else
            ASSERT(.false.)
        end if
    else
        ASSERT(.false.)
    end if
!
! --
!   verif 3 : pour chaque arete de la maille ima, on compare la longueur de son
!             projete orthogonal sur le segment qui constitue la fissure a la
!             longueur de ce segment. Au moins un de ces projetes doit avoir
!             une longueur comparable (critere relatif) a celle du segment.
! --
!
!   v_seg : vecteur associe au segment fissure, de norme norm_s
    v_seg(:) = vect2(:)-vect1(:)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    norm_s = ddot(b_n, v_seg, b_incx, v_seg, b_incy)
    norm_s = sqrt(norm_s)
    v_seg(:) = v_seg(:)/norm_s
!
!   boucle sur les aretes de ima
    l_crit = .false.
    do ia = 1, nbar
        na = ar(ia, 1)
        nb = ar(ia, 2)
        nunoa = zi(jconx1-1+zi(jconx2+ima-1)+na-1)
        nunob = zi(jconx1-1+zi(jconx2+ima-1)+nb-1)
!       v_are : vecteur associe a l'arete courante
        do i = 1, 3
            v_are(i) = vr_coo(3*(nunob-1)+i)-vr_coo(3*(nunoa-1)+i)
        end do
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pscal = ddot(b_n, v_are, b_incx, v_seg, b_incy)
!       v_pro : vecteur projete orthogonal de v_are sur v_seg, de norme norm_p
        v_pro(:) = pscal*v_seg(:)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        norm_p = ddot(b_n, v_pro, b_incx, v_pro, b_incy)
        norm_p = sqrt(norm_p)
!       ecart relatif entre norm_p et norm_s
        diffe = abs(norm_p-norm_s)/norm_s
        if (diffe .lt. crit) then
            l_crit = .true.
            exit
        end if
    end do
!
    if (.not. l_crit) then
        ASSERT(.false.)
    end if
!
    call jedema()
!
end subroutine
