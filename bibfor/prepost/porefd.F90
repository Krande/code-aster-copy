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

subroutine porefd(trange, noeu, cmp, nomrez)

    use DynaGene_module

    implicit none
#include "jeveux.h"
#include "nldef.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    character(len=*) :: trange, noeu, cmp, nomrez
!
!     POST-TRAITEMENT DE "RELA_EFFO_DEPL"
!
! ----------------------------------------------------------------------
    integer(kind=8) ::      nbpt, nbred, inume, nbnoli, nbvint, inl, start
    integer(kind=8) ::    i, ii, ic, nbpara, bloc_ini, i_bloc
    parameter(nbpara=8)
    real(kind=8) :: para(nbpara), xmax, temd, temf, temm
    complex(kind=8) :: c16b
    character(len=8) :: nomres, typara(nbpara), valek(3)
    character(len=16) :: nopara(nbpara), nomk16
    character(len=19) :: nomk19
    character(len=24) :: identifier
    integer(kind=8) :: nlin
    real(kind=8), pointer :: disc(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    integer(kind=8), pointer :: rdindx(:) => null()
!
    integer(kind=8), pointer :: nltype(:) => null()
    integer(kind=8), pointer :: vindx(:) => null()
    character(len=24), pointer :: nlname(:) => null()
    real(kind=8), pointer :: vint(:) => null()
    type(DynaGene) :: dyna_gene
!
    data nopara/'RELATION', 'NOEUD', 'CMP',&
     &              'PHASE', 'INST_INIT', 'INST_FIN',&
     &              'MAXI', 'INST_MAXI'/
    data typara/'K8', 'K8', 'K8', 'I', 'R', 'R', 'R', 'R'/
!     ------------------------------------------------------------------
!
    call jemarq()
    c16b = (0.d0, 0.d0)
    nomk16 = '                '
    nomk19 = '                   '
    nomk19(1:8) = trange
    nomk16(1:8) = trange
    nomres = nomrez

    call dyna_gene%init(trange(1:8))

!
    call tbcrsd(nomres, 'G')
    call tbajpa(nomres, nbpara, nopara, typara)
!
    call jeveuo(nomk19//'.DESC', 'L', vi=desc)
    nbnoli = desc(3)

    call jeveuo(nomk16//'.NL.TYPE', 'L', vi=nltype)
    call jeveuo(nomk16//'.NL.VIND', 'L', vi=vindx)
    call jeveuo(nomk16//'.NL.INTI', 'L', vk24=nlname)
    nbvint = vindx(nbnoli+1)-1

    AS_ALLOCATE(vi=rdindx, size=nbnoli)
    nbred = 0
    do i = 1, nbnoli
        if (nltype(i) .eq. NL_FX_RELATIONSHIP) then
            nbred = nbred+1
            rdindx(nbred) = i
        end if
    end do

    do inume = 1, nbred
        i = rdindx(inume)
        identifier = nlname((i-1)*5+1)
        if ((identifier(1:8) .eq. noeu) .and. (identifier(9:16) .eq. cmp)) goto 12
    end do
    call utmess('F', 'PREPOST4_57')

12  continue
    AS_DEALLOCATE(vi=rdindx)
    inl = i
!
    valek(1) = identifier(17:24)
    valek(2) = noeu
    valek(3) = cmp
!
!     --- RECHERCHE DU MAXIMUM DE LA FONCTION ---
    start = vindx(inl)

    call jelibe(nomk16//'.NL.TYPE')
    call jelibe(nomk16//'.NL.VIND')
    call jelibe(nomk16//'.NL.INTI')
!
!     --- RECHERCHE DES PHASES NON-LINEAIRE ---

    if (dyna_gene%n_bloc .eq. 0) then
        bloc_ini = 0
    else
        bloc_ini = 1
    end if

    do i_bloc = bloc_ini, dyna_gene%n_bloc
        call dyna_gene%get_values(dyna_gene%vint, i_bloc, length=nbpt, vr=vint)
        do i = 1, nbpt
            nlin = nint(vint((i-1)*nbvint+start-1+3))
            if (nlin .eq. 1) goto 20
        end do
    end do
    goto 500
!
20  continue
!
    ii = 0
    ic = 0
    do i_bloc = bloc_ini, dyna_gene%n_bloc
        call dyna_gene%get_values(dyna_gene%vint, i_bloc, length=nbpt, vr=vint)
        call dyna_gene%get_values(dyna_gene%disc, i_bloc, vr=disc)

        do i = 1, nbpt
            nlin = nint(vint((i-1)*nbvint+start-1+3))
            if (nlin .eq. 1 .and. ic .eq. 0) then
                xmax = vint((i-1)*nbvint+start-1+1)
                ic = 1
                ii = ii+1
                temd = disc(i)
                temm = disc(i)
            else if (nlin .eq. 1) then
                if (abs(vint((i-1)*nbvint+start-1+1)) .gt. abs(xmax)) then
                    xmax = vint((i-1)*nbvint+start-1+1)
                    temm = disc(i)
                end if
            else if (nlin .eq. 0 .and. ic .eq. 1) then
                ic = 0
                temf = disc(i)
                para(1) = temd
                para(2) = temf
                para(3) = xmax
                para(4) = temm
                call tbajli(nomres, nbpara, nopara, [ii], para, &
                            [c16b], valek, 0)
            end if
        end do
    end do
    if (ic .eq. 1) then
        temf = disc(nbpt)
        para(1) = temd
        para(2) = temf
        para(3) = xmax
        para(4) = temm
        call tbajli(nomres, nbpara, nopara, [ii], para, &
                    [c16b], valek, 0)
    end if
!
500 continue

    call dyna_gene%free()
!
    call jedema()
end subroutine
