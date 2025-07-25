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

subroutine mdinit(basemo, nbmode, nbnoli, depgen, vitgen, &
                  vint, ier, tinit, reprise, accgen, &
                  index)
!
    use DynaGene_module
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/extrac.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nlget.h"
#include "asterfort/nltype.h"
#include "asterfort/utmess.h"
    character(len=8) :: basemo
    integer(kind=8) :: nbmode, nbnoli
    real(kind=8) :: depgen(*), vitgen(*)
    real(kind=8), pointer :: vint(:)
    integer(kind=8) :: ier
    real(kind=8) :: tinit
    aster_logical, optional, intent(out) :: reprise
    real(kind=8), optional, intent(out) :: accgen(*)
    integer(kind=8), optional, intent(out) :: index
!
!
! DONNEES INITIALES
!
! IN  : BASEMO : NOM DU CONCEPT BASE MODALE
! IN  : NBMODE : NOMBRE DE MODES
! IN  : NBNOLI : NOMBRE DE NON LINEARITES
! OUT : DEPGEN : DEPLACEMENTS GENERALISES
! OUT : VITGEN : VITESSES GENERALISEES
! OUT : VINT   : VARIABLES INTERNES
!                (ON RETOURNE UNE VALEUR UNIQUEMENT SI nbnoli>0 ET QU'ON
!                 EST DANS UN CAS DE REPRISE)
! OUT : IER    : CODE RETOUR
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: im, ic
    character(len=19) :: nomdep, nomvit
    character(len=8) :: tran, crit, inter
! --------------------------------------------------------------------------------------------------
    integer(kind=8)               :: jdesc, jrefe, n1, jvind, vmessi(2)
    integer(kind=8)               :: nbinst, nc, ni, np, nt, nbvint, nbnoli0, nl_type, i_bloc, shift
    real(kind=8)          :: prec
    character(len=8)      :: sd_nl
    character(len=24)     :: no1_name, no2_name, vmessk(6), nltype_k0, nltype_k1
    integer(kind=8), pointer :: types(:) => null()
    real(kind=8), pointer :: acce(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    real(kind=8), pointer :: depl(:) => null()
    real(kind=8), pointer :: vite(:) => null()
    real(kind=8), pointer :: v_vint(:) => null()
    real(kind=8), pointer :: depi(:) => null()
    real(kind=8), pointer :: viti(:) => null()
    character(len=24), pointer :: inti(:) => null()
    type(DynaGene) :: dyna_gene
! --------------------------------------------------------------------------------------------------
    call jemarq()
    sd_nl = '&&OP29NL'
    ier = 0
    if (present(reprise)) reprise = .false.
!
    nbvint = 0
    if (nbnoli .gt. 0) nbvint = size(vint)
!
!     --- DEPLACEMENT ---
    call getvid('ETAT_INIT', 'DEPL', iocc=1, scal=nomdep, nbret=n1)
    if (n1 .ne. 0) then
        call jeveuo(nomdep//'.VALE', 'L', vr=depi)
!
!        --- VERIF COMPATIBILITE DES BASES DE PROJECTION
        call jeveuo(nomdep//'.REFE', 'L', jrefe)
        if (zk24(jrefe) (1:8) .ne. basemo) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_42')
        end if
        call jeveuo(nomdep//'.DESC', 'L', jdesc)
        if (zi(jdesc+1) .ne. nbmode) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_43')
        end if
        do im = 1, nbmode
            depgen(im) = depi(im)
        end do
    end if
!
!     --- VITESSE ---
    call getvid('ETAT_INIT', 'VITE', iocc=1, scal=nomvit, nbret=n1)
    if (n1 .ne. 0) then
        call jeveuo(nomvit//'.VALE', 'L', vr=viti)
!
!        --- VERIF COMPATIBILITE DES BASES DE PROJECTION
        call jeveuo(nomvit//'.REFE', 'L', jrefe)
        if (zk24(jrefe) (1:8) .ne. basemo) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_42')
        end if
        call jeveuo(nomvit//'.DESC', 'L', jdesc)
        if (zi(jdesc+1) .ne. nbmode) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_43')
        end if
        do im = 1, nbmode
            vitgen(im) = viti(im)
        end do
    end if
!
!     --- CAS D UNE REPRISE ---
    call getvid('ETAT_INIT', 'RESULTAT', iocc=1, scal=tran, nbret=nt)
    if (nt .ne. 0) then
!       recuperation du descripteur
        call jeveuo(tran//'           .DESC', 'L', jdesc)
        nbnoli0 = zi(jdesc+2)
        if (nbnoli0 .ne. nbnoli) then
            vmessi(1) = nbnoli0
            vmessi(2) = nbnoli
            call utmess('F', 'ALGORITH5_82', ni=2, vali=vmessi)
        end if

        if (nbnoli .ne. 0) then
!           recuperation des donnees sur les non linearites
            call jeveuo(tran//'        .NL.INTI', 'L', vk24=inti)
            call jeveuo(tran//'        .NL.TYPE', 'L', vi=types)
            do ic = 1, nbnoli
                nltype_k0 = nltype(types(ic))
                call nlget(sd_nl, _NL_TYPE, iocc=ic, iscal=nl_type)
                nltype_k1 = nltype(nl_type)
                call nlget(sd_nl, _NO1_NAME, iocc=ic, kscal=no1_name)
                call nlget(sd_nl, _NO2_NAME, iocc=ic, kscal=no2_name)
                if ((nltype_k0 .ne. nltype_k1) &
                    .or. (inti((ic-1)*5+2) .ne. no1_name) &
                    .or. (inti((ic-1)*5+3) .ne. no2_name)) then
                    vmessk(1) = nltype_k0
                    vmessk(2) = nltype_k1
                    vmessk(3) = inti((ic-1)*5+2)
                    vmessk(4) = no1_name
                    vmessk(5) = inti((ic-1)*5+3)
                    vmessk(6) = no2_name
                    call utmess('F', 'ALGORITH5_83', nk=6, valk=vmessk)
                end if
            end do
        end if
!       recuperation des champs depl vite vint
!           les calculs sont donc deja fait au pas de recuperation
        call getvtx('ETAT_INIT', 'CRITERE', iocc=1, scal=crit, nbret=nc)
        call getvr8('ETAT_INIT', 'PRECISION', iocc=1, scal=prec, nbret=np)
        if (nc .eq. 0) crit = 'RELATIF'
        if (np .eq. 0) prec = 1.d-6

        call getvr8('ETAT_INIT', 'INST_INIT', iocc=1, scal=tinit, nbret=ni)

        call dyna_gene%init(tran)

        if (ni .eq. 0) then
            call dyna_gene%get_values(dyna_gene%disc, dyna_gene%n_bloc, shift, nbinst, vr=disc)
            tinit = disc(nbinst)
        else
            call dyna_gene%get_values_by_disc(dyna_gene%disc, tinit, shift, nbinst, vr=disc)
        end if
!       Deplacement
        inter = 'NON'

        call dyna_gene%get_current_bloc(dyna_gene%disc, i_bloc)
        call dyna_gene%get_values(dyna_gene%depl, i_bloc, vr=depl)
        call extrac(inter, prec, crit, nbinst, disc, &
                    tinit, depl, nbmode, depgen, ier, index)
        index = shift+index
        if (ier .ne. 0) then
            call utmess('F', 'ALGORITH5_46')
        end if
!       Vitesse
        call dyna_gene%get_values(dyna_gene%vite, i_bloc, vr=vite)
        inter = 'NON'
        call extrac(inter, prec, crit, nbinst, disc, &
                    tinit, vite, nbmode, vitgen, ier)
        if (ier .ne. 0) then
            call utmess('F', 'ALGORITH5_47')
        end if
!       Acceleration
        if (present(accgen)) then
            call dyna_gene%get_values(dyna_gene%acce, i_bloc, vr=acce)
            inter = 'NON'
            call extrac(inter, prec, crit, nbinst, disc, &
                        tinit, acce, nbmode, accgen, ier)
            if (ier .ne. 0) then
                call utmess('F', 'ALGORITH5_47')
            end if
        end if
!       Variables internes
        if (nbnoli .gt. 0) then
            call jeveuo(tran//'        .NL.VIND', 'L', jvind)
            nbvint = zi(jvind+nbnoli)-1
            call dyna_gene%get_values(dyna_gene%vint, i_bloc, vr=v_vint)
            inter = 'NON'
            call extrac(inter, prec, crit, nbinst, disc, &
                        tinit, v_vint, nbvint, vint, ier)
            if (ier .ne. 0) then
                call utmess('F', 'ALGORITH5_48')
            end if
        end if

        call dyna_gene%free

        if (present(reprise)) reprise = .true.
    end if
!
    call jedema()
!
end subroutine
