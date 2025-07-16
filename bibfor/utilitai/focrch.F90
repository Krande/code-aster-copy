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
subroutine focrch(nomfon, resu, noeud, parax, paray, &
                  base, i0, intitu, ind, listr, &
                  sst, nsst, ier)
!
    use DynaGene_module
!
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
    integer(kind=8) :: nsst, i0, ind, ier
    character(len=1) :: base
    character(len=16) :: parax, paray
    character(len=8) :: sst, noeud
    character(len=19) :: nomfon, resu, listr
    character(len=24) :: intitu
!
!     RECUPERATION D'UNE FONCTION DANS UNE STRUCTURE "TRAN_GENE"
!     POUR UN NOEUD DE CHOC
!     ------------------------------------------------------------------
! IN  : NOMFON : NOM DE LA FONCTION
! IN  : RESU   : NOM DE LA STRUCTURE RESULTAT
! IN  : NOEUD  : NOEUD DE CHOC
! IN  : PARAX  : PARAMETRE DE LA FONCTION EN X
! IN  : PARAY  : PARAMETRE DE LA FONCTION EN Y
! IN  : BASE   : 'GLOBALE'
! IN  : INT    : PRISE EN COMPTE D'UN NOM DE LIAISON
! IN  : INTITU : NOM D'UNE LIAISON
! IN  : IND    : PRISE EN COMPTE D'UNE LISTE DE PARAMETRES
! IN  : LISTR  : LISTE DE PARAMETRES
! IN  : SST    : NOM DE LA SOUS-STRUCTURE
! IN  : NSST   : PRISE EN COMPTE DU NOM D'UNE SOUS-STRUCTURE
! OUT : IER    : CODE RETOUR, = 0 : OK
!     ------------------------------------------------------------------
    character(len=8) :: k8b
    character(len=16) :: nomcmd
    character(len=19) :: fonct1, fonct2
!     ----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ic, inl, ie, ival, jinst, iparax
    integer(kind=8) :: jparx, jpary, jval, jvalx, jvaly, start, nbvint, iparay
    integer(kind=8) :: lfon, lg, lpro, lval, nbnoli, nbinst, nbpara
    integer(kind=8) :: nbval, jvint, i_bloc, shift, length
    character(len=24), pointer :: nlname(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    integer(kind=8), pointer :: vindx(:) => null()
    real(kind=8), pointer :: vr(:) => null()
    type(DynaGene) :: dyna_gene
    blas_int :: b_incx, b_incy, b_n
!
!-----------------------------------------------------------------------
    call jemarq()
    ier = 999
    call getres(k8b, k8b, nomcmd)
!
    call jeveuo(resu//'.DESC', 'L', vi=desc)
    nbnoli = desc(3)
!
    call jeveuo(resu(1:16)//'.NL.INTI', 'L', vk24=nlname)
!
    ic = 1
    if (i0 .ne. 0) then
        do inl = 1, nbnoli
            if (nlname((inl-1)*5+1) .eq. intitu) goto 4
        end do
        call utmess('A', 'UTILITAI_86', sk=intitu)
        goto 999
4       continue
        if (nsst .eq. 0) then
            if (nlname((inl-1)*5+2) (1:8) .eq. noeud) goto 16
            ic = 2
            if (nlname((inl-1)*5+3) (1:8) .eq. noeud) goto 16
            lg = max(1, lxlgut(noeud))
            call utmess('A', 'UTILITAI_87', sk=noeud(1:lg))
            goto 999
        else
            if (nlname((inl-1)*5+4) (1:8) .eq. sst) goto 116
            if (nlname((inl-1)*5+5) (1:8) .eq. sst) goto 116
            call utmess('A', 'UTILITAI_88')
            goto 999
116         continue
            if (nlname((inl-1)*5+2) .ne. noeud .and. nlname((inl-1)*5+3) .ne. noeud) then
                lg = max(1, lxlgut(noeud))
                call utmess('A', 'UTILITAI_89', sk=noeud(1:lg))
                goto 999
            end if
            if (nlname((inl-1)*5+2) (1:8) .eq. noeud .and. nlname((inl-1)*5+4) (1:8) .eq. sst) &
                goto 16
            ic = 2
            if (nlname((inl-1)*5+3) (1:8) .eq. noeud .and. nlname((inl-1)*5+5) (1:8) .eq. sst) &
                goto 16
            lg = max(1, lxlgut(noeud))
            call utmess('A', 'UTILITAI_90', sk=noeud(1:lg))
            goto 999
        end if
    end if
!     --- RECHERCHE DU NOEUD_1 DE CHOC ---
    do inl = 1, nbnoli
        if (nlname((inl-1)*5+2) (1:8) .eq. noeud) goto 16
    end do
!     --- RECHERCHE DU NOEUD_2 DE CHOC ---
    ic = 2
    do inl = 1, nbnoli
        if (nlname((inl-1)*5+3) (1:8) .eq. noeud) goto 16
    end do
    lg = max(1, lxlgut(noeud))
    call utmess('A', 'UTILITAI_87', sk=noeud(1:lg))
    goto 999
16  continue
!
    call dyna_gene%init(resu(1:8))
    nbinst = dyna_gene%length
!
    if (dyna_gene%n_bloc .eq. 0) then
        call jeveuo(resu(1:16)//'.NL.VINT', 'L', jvint)
    else
        jvint = 0
    end if
    call jeveuo(resu(1:16)//'.NL.VIND', 'L', vi=vindx)
    start = vindx(inl)-1
    nbvint = vindx(nbnoli+1)-1
!
    if (dyna_gene%n_bloc .eq. 0) then
        call jeveuo(resu//'.DISC', 'L', jinst)
    else
        if (parax(1:4) .eq. 'INST' .or. paray(1:4) .eq. 'INST' .or. ind .ne. 0) then
            call wkvect('&&FOCRCH.INST', 'V V R', nbinst, jinst)
            do i_bloc = 1, dyna_gene%n_bloc
                call dyna_gene%get_values(dyna_gene%disc, i_bloc, shift, length, vr=vr)
                b_n = to_blas_int(length)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vr, b_incx, zr(jinst+shift), b_incy)
            end do
        end if
    end if
!
    if (parax(1:4) .eq. 'INST') then
        jvalx = jinst
        goto 20
    else if (parax(1:2) .eq. 'FN') then
        jparx = jvint+start
    else if (parax(1:3) .eq. 'FT1') then
        jparx = jvint+start+1
    else if (parax(1:3) .eq. 'FT2') then
        jparx = jvint+start+2
    else if (parax(1:2) .eq. 'VN') then
        jparx = jvint+start+3
    else if (parax(1:3) .eq. 'VT1') then
        jparx = jvint+start+4
    else if (parax(1:3) .eq. 'VT2') then
        jparx = jvint+start+5
    else if (parax(1:5) .eq. 'DXLOC') then
        if (ic .eq. 1) then
            jparx = jvint+start+6
        else
            jparx = jvint+start+9
        end if
    else if (parax(1:5) .eq. 'DYLOC') then
        if (ic .eq. 1) then
            jparx = jvint+start+7
        else
            jparx = jvint+start+10
        end if
    else if (parax(1:5) .eq. 'DZLOC') then
        if (ic .eq. 1) then
            jparx = jvint+start+8
        else
            jparx = jvint+start+11
        end if
    else if (parax(1:4) .eq. 'VINT') then
        read (parax(5:7), '(I2)') iparax
        jparx = jvint+start+iparax-1
    else
        lg = max(1, lxlgut(parax(1:8)))
        call utmess('A', 'UTILITAI_91', sk=parax(1:lg))
        goto 999
    end if
    call wkvect('&&FOCRCH.PARAX', 'V V R', nbinst, jvalx)
    if (dyna_gene%n_bloc .eq. 0) then
        b_n = to_blas_int(nbinst)
        b_incx = to_blas_int(nbvint)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jparx), b_incx, zr(jvalx), b_incy)
    else
        do i_bloc = 1, dyna_gene%n_bloc
            call dyna_gene%get_values(dyna_gene%vint, i_bloc, shift, length, vr=vr)
            vr => vr(1+jparx:)
            b_n = to_blas_int(length)
            b_incx = to_blas_int(nbvint)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vr, b_incx, zr(jvalx+shift), b_incy)
        end do
    end if
20  continue
!
    if (paray(1:4) .eq. 'INST') then
        jvaly = jinst
        goto 22
    else if (paray(1:2) .eq. 'FN') then
        jpary = jvint+start
    else if (paray(1:3) .eq. 'FT1') then
        jpary = jvint+start+1
    else if (paray(1:3) .eq. 'FT2') then
        jpary = jvint+start+2
    else if (paray(1:2) .eq. 'VN') then
        jpary = jvint+start+3
    else if (paray(1:3) .eq. 'VT1') then
        jpary = jvint+start+4
    else if (paray(1:3) .eq. 'VT2') then
        jpary = jvint+start+5
    else if (paray(1:5) .eq. 'DXLOC') then
        if (ic .eq. 1) then
            jpary = jvint+start+6
        else
            jpary = jvint+start+9
        end if
    else if (paray(1:5) .eq. 'DYLOC') then
        if (ic .eq. 1) then
            jpary = jvint+start+7
        else
            jpary = jvint+start+10
        end if
    else if (paray(1:5) .eq. 'DZLOC') then
        if (ic .eq. 1) then
            jpary = jvint+start+8
        else
            jpary = jvint+start+11
        end if
    else if (paray(1:4) .eq. 'VINT') then
        read (paray(5:7), '(I3)') iparay
        jpary = jvint+start+iparay-1
    else
        lg = max(1, lxlgut(paray(1:8)))
        call utmess('A', 'UTILITAI_91', sk=paray(1:lg))
        goto 999
    end if
    call wkvect('&&FOCRCH.PARAY', 'V V R', nbinst, jvaly)
    if (dyna_gene%n_bloc .eq. 0) then
        b_n = to_blas_int(nbinst)
        b_incx = to_blas_int(nbvint)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jpary), b_incx, zr(jvaly), b_incy)
    else
        do i_bloc = 1, dyna_gene%n_bloc
            call dyna_gene%get_values(dyna_gene%vint, i_bloc, shift, length, vr=vr)
            vr => vr(1+jpary:)
            b_n = to_blas_int(length)
            b_incx = to_blas_int(nbvint)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vr, b_incx, zr(jvaly+shift), b_incy)
        end do
    end if
22  continue
!
    if (ind .eq. 0) then
        ASSERT(lxlgut(nomfon) .le. 24)
        call wkvect(nomfon//'.PROL', base//' V K24', 6, lpro)
        zk24(lpro) = 'FONCTION'
        zk24(lpro+1) = 'LIN LIN '
        zk24(lpro+2) = parax
        zk24(lpro+3) = paray
        zk24(lpro+4) = 'EE'
        zk24(lpro+5) = nomfon
!
        nbval = nbinst*2
        call wkvect(nomfon//'.VALE', base//' V R', nbval, lval)
        lfon = lval+nbinst
        do ival = 0, nbinst-1
            zr(lval+ival) = zr(jvalx+ival)
            zr(lfon+ival) = zr(jvaly+ival)
        end do
        ier = 0
!
    else
        fonct1 = '&&FOCRCH.FONCT1'
        ASSERT(lxlgut(fonct1) .le. 24)
        call wkvect(fonct1//'.PROL', 'V V K24', 6, lpro)
        zk24(lpro) = 'FONCTION'
        zk24(lpro+1) = 'LIN LIN '
        zk24(lpro+2) = 'INST'
        zk24(lpro+3) = parax
        zk24(lpro+4) = 'EE'
        zk24(lpro+5) = fonct1
        nbval = nbinst*2
        call wkvect(fonct1//'.VALE', 'V V R', nbval, lval)
        lfon = lval+nbinst
        do ival = 0, nbinst-1
            zr(lval+ival) = zr(jinst+ival)
            zr(lfon+ival) = zr(jvalx+ival)
        end do
!
        fonct2 = '&&FOCRCH.FONCT2'
        ASSERT(lxlgut(fonct2) .le. 24)
        call wkvect(fonct2//'.PROL', 'V V K24', 6, lpro)
        zk24(lpro) = 'FONCTION'
        zk24(lpro+1) = 'LIN LIN '
        zk24(lpro+2) = 'INST'
        zk24(lpro+3) = paray
        zk24(lpro+4) = 'EE'
        zk24(lpro+5) = fonct2
        nbval = nbinst*2
        call wkvect(fonct2//'.VALE', 'V V R', nbval, lval)
        lfon = lval+nbinst
        do ival = 0, nbinst-1
            zr(lval+ival) = zr(jinst+ival)
            zr(lfon+ival) = zr(jvaly+ival)
        end do
!
        call jeveuo(listr//'.VALE', 'L', jval)
        call jelira(listr//'.VALE', 'LONUTI', nbpara)
!
        ASSERT(lxlgut(nomfon) .le. 24)
        call wkvect(nomfon//'.PROL', base//' V K24', 6, lpro)
        zk24(lpro) = 'FONCTION'
        zk24(lpro+1) = 'LIN LIN '
        zk24(lpro+2) = parax
        zk24(lpro+3) = paray
        zk24(lpro+4) = 'EE'
        zk24(lpro+5) = nomfon
!
        nbval = nbpara*2
        call wkvect(nomfon//'.VALE', base//' V R', nbval, lval)
        lfon = lval+nbpara
        do ival = 0, nbpara-1
            call fointe('F ', fonct1, 1, 'INST', zr(jval+ival), &
                        zr(lval+ival), ie)
            call fointe('F ', fonct2, 1, 'INST', zr(jval+ival), &
                        zr(lfon+ival), ie)
        end do
!
        call jedetr(fonct1//'.PROL')
        call jedetr(fonct1//'.VALE')
        call jedetr(fonct2//'.PROL')
        call jedetr(fonct2//'.VALE')
        ier = 0
    end if
    if (parax(1:4) .ne. 'INST') call jedetr('&&FOCRCH.PARAX')
    if (paray(1:4) .ne. 'INST') call jedetr('&&FOCRCH.PARAY')
    if ((parax(1:4) .eq. 'INST' .or. paray(1:4) .eq. 'INST' .or. ind .ne. 0) .and. &
        dyna_gene%n_bloc .ne. 0) call jedetr('&&FOCRCH.INST')
!
    call dyna_gene%free()
999 continue
    call jedema()
end subroutine
