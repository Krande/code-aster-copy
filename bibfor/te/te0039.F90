! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine te0039(option, nomte)
!
use te0047_type
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dis_choc_frot_syme.h"
#include "asterfort/discret_sief.h"
#include "asterfort/infdis.h"
#include "asterfort/infted.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/terefe.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/lteatt.h"
#include "blas/dcopy.h"
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! IN OPTION    : K16 :  OPTION DE CALCUL
!      'SIEQ_ELNO'       'SIEQ_ELGA'
!     'FORC_NODA'  'REFE_FORC_NODA'
!     'FONL_NOEU'
! IN NOMTE     : K16 : NOM DU TYPE ELEMENT
!     DISCRETS :
!        'MECA_DIS_T_N'      'MECA_DIS_T_L'     'MECA_DIS_TR_N'
!        'MECA_DIS_TR_L'     'MECA_2D_DIS_T_N'  'MECA_2D_DIS_T_L'
!        'MECA_2D_DIS_TR_N'  'MECA_2D_DIS_TR_L'
! --------------------------------------------------------------------------------------------------
!
    type(te0047_dscr) :: for_discret
!
!   Les variables internes
    integer, parameter  :: nbvari=9
    real(kind=8)        :: varmo(nbvari), varpl(nbvari)
!
    integer :: lorien, lmater, ii, jj
    integer :: ivectu, icontg, neq
    integer :: iplouf, infodi, itype, ibid
    integer :: igeom, ideplm, ideplp, icompo, jdc, irep, ifono, ilogic
!
    real(kind=8) :: pgl(3, 3), force(3)
    real(kind=8) :: fs(12), ugp(12), dug(12), ulp(12), dul(12), dvl(12), dpe(12), dve(12)
    real(kind=8) :: sim(12), sip(12), fono(12)
    real(kind=8) :: klv(78)
    real(kind=8) :: forref, momref
    real(kind=8) :: r8bid
!
    character(len=8)  :: k8bid
    character(len=16) :: kmess(5)
!
    aster_logical, parameter    :: Predic=ASTER_FALSE
!
! --------------------------------------------------------------------------------------------------
    infodi = 1
!   On verifie que les caracteristiques ont ete affectees
!   Le code du discret
    call infdis('CODE', ibid, r8bid, nomte)
!   Le code stoké dans la carte
    call infdis('TYDI', infodi, r8bid, k8bid)
    if (infodi .ne. ibid) then
        call utmess('F+', 'DISCRETS_25', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'F+')
    endif
!   Discret de type raideur
    call infdis('DISK', infodi, r8bid, k8bid)
    if (infodi .eq. 0) then
        call utmess('A+', 'DISCRETS_27', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'A+')
    endif
!   Matrice de raideur symetrique ou pas, pour les discrets
    call infdis('SYMK', infodi, r8bid, k8bid)
!   Récupere les informations sur les elements
    for_discret%option = option
    for_discret%nomte  = nomte
    call infted(for_discret%nomte, infodi, &
                for_discret%nbt, for_discret%nno, for_discret%nc, for_discret%ndim, itype)
    neq = for_discret%nno*for_discret%nc
!
    if (option(1:14) .eq. 'REFE_FORC_NODA') then
        call jevech('PVECTUR', 'E', ivectu)
        if (lteatt('MODELI','DTR')) then
            call terefe('EFFORT_REFE', 'MECA_DISCRET', forref)
            call terefe('MOMENT_REFE', 'MECA_DISCRET', momref)
            do  ii = 1, for_discret%nno
                do jj = 1, 3
                    zr(ivectu+(ii-1)*for_discret%nc+jj-1)=forref
                enddo
                do jj = 4, for_discret%nc
                    zr(ivectu+(ii-1)*for_discret%nc+jj-1)=momref
                enddo
            enddo
        else if (lteatt('MODELI','2DT')) then
            call terefe('EFFORT_REFE', 'MECA_DISCRET', forref)
            do  ii = 1, for_discret%nno
                zr(ivectu+(ii-1)*for_discret%nc)=forref
                zr(ivectu+(ii-1)*for_discret%nc+1)=forref
            enddo
        else if (lteatt('MODELI','2TR')) then
            call terefe('EFFORT_REFE', 'MECA_DISCRET', forref)
            call terefe('MOMENT_REFE', 'MECA_DISCRET', momref)
            do ii = 1, for_discret%nno
                zr(ivectu+(ii-1)*for_discret%nc)=forref
                zr(ivectu+(ii-1)*for_discret%nc+1)=forref
                zr(ivectu+(ii-1)*for_discret%nc+2)=momref
            enddo
        else if (lteatt('MODELI','DIT')) then
            call terefe('EFFORT_REFE', 'MECA_DISCRET', forref)
            do ii = 1, for_discret%nno
                zr(ivectu+(ii-1)*for_discret%nc)=forref
                zr(ivectu+(ii-1)*for_discret%nc+1)=forref
                zr(ivectu+(ii-1)*for_discret%nc+2)=forref
            enddo
        else
            kmess(1) = option
            kmess(2) = nomte
            kmess(3) = 'TE0039'
            call utmess('F', 'DISCRETS_15', nk=2, valk=kmess)
        endif
    else if (option.eq.'FONL_NOEU') then
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PCOMPOR', 'L', icompo)
        call jevech('PMATERC', 'L', lmater)
        if (lteatt('MODELI','DTR') .or. lteatt('MODELI','DIT')) then
            !   PARAMETRES EN ENTREE
            call jevech('PCAORIE', 'L', lorien)
            call matrot(zr(lorien), for_discret%pgl)
            ! DEPLACEMENTS DANS LE REPERE GLOBAL
            !   UGM = DEPLACEMENT PRECEDENT
            !   DUG = INCREMENT DE DEPLACEMENT
            !   UGP = DEPLACEMENT COURANT
            do ii = 1, neq
                dug(ii) = zr(ideplp+ii-1)
                ugp(ii) = zr(ideplm+ii-1) + dug(ii)
            enddo
            ! Deplacements dans le repere local
            !   ULM = DEPLACEMENT PRECEDENT    = PLG * UGM
            !   DUL = INCREMENT DE DEPLACEMENT = PLG * DUG
            !   ULP = DEPLACEMENT COURANT      = PLG * UGP
            if (for_discret%ndim .eq. 3) then
                call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, dug, dul)
                call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, ugp, ulp)
            else if (for_discret%ndim.eq.2) then
                call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, dug, dul)
                call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, ugp, ulp)
            endif
            ! Seul le cas symetrique est traite
            call infdis('SYMK', iplouf, r8bid, k8bid)
            if (iplouf .ne. 1) then
                kmess(1) = option
                kmess(2) = nomte
                kmess(3) = 'TE0039'
                kmess(4) = ' '
                call utmess('F', 'DISCRETS_12', nk=4, valk=kmess)
            endif
            !
            call jevech('PCADISK', 'L', jdc)
            call infdis('REPK', irep, r8bid, k8bid)
            call dcopy(for_discret%nbt, zr(jdc), 1, klv, 1)
            if (irep .eq. 1) then
                call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
            endif
            if (zk16(icompo) .eq. 'DIS_CHOC') then
                varmo(:) = 0.0; dvl(:)   = 0.0; dpe(:)   = 0.0; dve(:)   = 0.0
                ! Relation de comportement de choc : forces nodales
                call jevech('PVECTUR', 'E', ifono)
                do  ii = 1, neq
                    zr(ifono+ii-1) = 0.0
                    sim(ii)        = 0.0
                enddo
                !
                ilogic = 0; force(1:3) = 0.0
                call discret_sief(for_discret, klv, ulp, sim, ilogic, sip, fono, force)
                call dis_choc_frot_syme(for_discret, zi(lmater), ulp, zr(igeom), klv, &
                                        dvl, dpe, dve, Predic, force, varmo, varpl)
                ilogic = 2
                call discret_sief(for_discret, klv, ulp, sim, ilogic, sip, zr(ifono), force)
                do ii = 1, neq
                    zr(ifono+ii-1) = zr(ifono+ii-1)-fono(ii)
                enddo
                if (for_discret%nno .eq. 2) then
                    do ii = 1, for_discret%nc
                        zr(ifono+ii-1) = 0.0
                    enddo
                endif
            endif
        endif
    else if (option .eq. 'FORC_NODA') then
        call jevech('PCONTMR', 'L', icontg)
        call jevech('PVECTUR', 'E', ivectu)
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                fs(ii) = zr(icontg+ii-1)
            enddo
        else
            do  ii = 1, for_discret%nc
                fs(ii)                = -zr(icontg+ii-1)
                fs(ii+for_discret%nc) =  zr(icontg+ii+for_discret%nc-1)
            enddo
        endif
        call jevech('PCAORIE', 'L', lorien)
        !
        call matrot(zr(lorien), pgl)
        if (for_discret%ndim .eq. 3) then
            call utpvlg(for_discret%nno, for_discret%nc, pgl, fs, zr(ivectu))
        else if (for_discret%ndim.eq.2) then
            call ut2vlg(for_discret%nno, for_discret%nc, pgl, fs, zr(ivectu))
        endif
    else
        ASSERT(.false.)
    endif
!
end subroutine
