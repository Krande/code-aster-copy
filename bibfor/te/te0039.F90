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
subroutine te0039(option, nomte)
!
    use te0047_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
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
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/lteatt.h"
#include "blas/dcopy.h"
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! IN OPTION    : K16 :  OPTION DE CALCUL
!     'SIEQ_ELNO'  'SIEQ_ELGA'
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
    integer(kind=8), parameter :: nbvari = 9
    real(kind=8) :: varmo(nbvari), varpl(nbvari)
!
    integer(kind=8) :: lorien, lmater, ii, jj
    integer(kind=8) :: ivectu, jvSief, neq
    integer(kind=8) :: iplouf, infodi, itype, ibid
    integer(kind=8) :: igeom, jdc, irep, ifono, ilogic, jvDisp
!
    real(kind=8) :: pgl(3, 3), force(3)
    real(kind=8) :: fs(12), ugp(12), ulp(12), dvl(12), dpe(12), dve(12)
    real(kind=8) :: sim(12), sip(12), fono(12)
    real(kind=8) :: klv(78), kgv(78)
    real(kind=8) :: forref, momref
    real(kind=8) :: r8bid
!
    character(len=8) :: k8bid
    character(len=16) :: kmess(5)
    character(len=16), pointer :: compor(:) => null()
!
    aster_logical, parameter :: Predic = ASTER_FALSE
    blas_int :: b_incx, b_incy, b_n
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
    end if
!   Discret de type raideur
    call infdis('DISK', infodi, r8bid, k8bid)
    if (infodi .eq. 0) then
        call utmess('A+', 'DISCRETS_27', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'A+')
    end if
!   Matrice de raideur symetrique ou pas, pour les discrets
    call infdis('SYMK', infodi, r8bid, k8bid)
!   Récupere les informations sur les elements
    for_discret%option = option
    for_discret%nomte = nomte
    call infted(for_discret%nomte, infodi, for_discret%nbt, for_discret%nno, for_discret%nc, &
                for_discret%ndim, itype)
    neq = for_discret%nno*for_discret%nc
!
    if (option(1:14) .eq. 'REFE_FORC_NODA') then
        call jevech('PVECTUR', 'E', ivectu)
        if (lteatt('MODELI', 'DTR')) then
            call terefe('EFFORT_REFE', 'MECA_DISCRET', forref)
            call terefe('MOMENT_REFE', 'MECA_DISCRET', momref)
            do ii = 1, for_discret%nno
                do jj = 1, 3
                    zr(ivectu+(ii-1)*for_discret%nc+jj-1) = forref
                end do
                do jj = 4, for_discret%nc
                    zr(ivectu+(ii-1)*for_discret%nc+jj-1) = momref
                end do
            end do
        else if (lteatt('MODELI', '2DT')) then
            call terefe('EFFORT_REFE', 'MECA_DISCRET', forref)
            do ii = 1, for_discret%nno
                zr(ivectu+(ii-1)*for_discret%nc) = forref
                zr(ivectu+(ii-1)*for_discret%nc+1) = forref
            end do
        else if (lteatt('MODELI', '2TR')) then
            call terefe('EFFORT_REFE', 'MECA_DISCRET', forref)
            call terefe('MOMENT_REFE', 'MECA_DISCRET', momref)
            do ii = 1, for_discret%nno
                zr(ivectu+(ii-1)*for_discret%nc) = forref
                zr(ivectu+(ii-1)*for_discret%nc+1) = forref
                zr(ivectu+(ii-1)*for_discret%nc+2) = momref
            end do
        else if (lteatt('MODELI', 'DIT')) then
            call terefe('EFFORT_REFE', 'MECA_DISCRET', forref)
            do ii = 1, for_discret%nno
                zr(ivectu+(ii-1)*for_discret%nc) = forref
                zr(ivectu+(ii-1)*for_discret%nc+1) = forref
                zr(ivectu+(ii-1)*for_discret%nc+2) = forref
            end do
        else
            kmess(1) = option
            kmess(2) = nomte
            kmess(3) = 'TE0039'
            call utmess('F', 'DISCRETS_15', nk=2, valk=kmess)
        end if
    else if (option .eq. 'FONL_NOEU') then
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PDEPLAR', 'L', jvDisp)
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PMATERC', 'L', lmater)
        if (lteatt('MODELI', 'DTR') .or. lteatt('MODELI', 'DIT')) then
!   PARAMETRES EN ENTREE
            call jevech('PCAORIE', 'L', lorien)
            call matrot(zr(lorien), for_discret%pgl)
! DÉPLACEMENTS DANS LE REPÈRE GLOBAL
!   UGM = DEPLACEMENT PRECEDENT
!   UGP = DEPLACEMENT COURANT
            do ii = 1, neq
                ugp(ii) = zr(jvDisp+ii-1)
            end do
! DÉPLACEMENTS DANS LE REPÈRE LOCAL
!   ULM = DEPLACEMENT PRECEDENT    = PLG * UGM
!   DUL = INCREMENT DE DEPLACEMENT = PLG * DUG
!   ULP = DEPLACEMENT COURANT      = PLG * UGP
            if (for_discret%ndim .eq. 3) then
                call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, ugp, ulp)
            else if (for_discret%ndim .eq. 2) then
                call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, ugp, ulp)
            end if
! Seul le cas symetrique est traite
            call infdis('SYMK', iplouf, r8bid, k8bid)
            if (iplouf .ne. 1) then
                kmess(1) = option
                kmess(2) = nomte
                kmess(3) = 'TE0039'
                kmess(4) = ' '
                call utmess('F', 'DISCRETS_12', nk=4, valk=kmess)
            end if
            !
            call jevech('PCADISK', 'L', jdc)
            call infdis('REPK', irep, r8bid, k8bid)
! irep = 1 = matrice en repère global
            if (irep .eq. 1) then
! Matrice dans le repère global
                b_n = to_blas_int(for_discret%nbt)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(jdc), b_incx, kgv, b_incy)
! Matrice dans le repère local
                call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, kgv, klv)
            else
! Matrice dans le repère local
                b_n = to_blas_int(for_discret%nbt)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
! Matrice dans le repère global
                call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, kgv)
            end if
            !
            if (compor(RELA_NAME) .eq. 'DIS_CHOC') then
                varmo(:) = 0.0; dvl(:) = 0.0; dpe(:) = 0.0; dve(:) = 0.0
! Relation de comportement de choc : forces nodales
                call jevech('PVECTUR', 'E', ifono)
                do ii = 1, neq
                    zr(ifono+ii-1) = 0.0
                    sim(ii) = 0.0
                end do
                !
                ilogic = 0; force(1:3) = 0.0
                call discret_sief(for_discret, klv, ulp, sim, ilogic, &
                                  sip, fono, force)
                call dis_choc_frot_syme(for_discret, zi(lmater), ulp, zr(igeom), klv, &
                                        kgv, dvl, dpe, dve, Predic, &
                                        force, varmo, varpl)
                ilogic = 2
                call discret_sief(for_discret, klv, ulp, sim, ilogic, &
                                  sip, zr(ifono), force)
                do ii = 1, neq
                    zr(ifono+ii-1) = zr(ifono+ii-1)-fono(ii)
                end do
                if (for_discret%nno .eq. 2) then
                    do ii = 1, for_discret%nc
                        zr(ifono+ii-1) = 0.0
                    end do
                end if
            end if
        end if
    else if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PVECTUR', 'E', ivectu)
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                fs(ii) = zr(jvSief+ii-1)
            end do
        else
            do ii = 1, for_discret%nc
                fs(ii) = -zr(jvSief+ii-1)
                fs(ii+for_discret%nc) = zr(jvSief+ii+for_discret%nc-1)
            end do
        end if
        call jevech('PCAORIE', 'L', lorien)
        !
        call matrot(zr(lorien), pgl)
        if (for_discret%ndim .eq. 3) then
            call utpvlg(for_discret%nno, for_discret%nc, pgl, fs, zr(ivectu))
        else if (for_discret%ndim .eq. 2) then
            call ut2vlg(for_discret%nno, for_discret%nc, pgl, fs, zr(ivectu))
        end if
    else
        ASSERT(.false.)
    end if
!
end subroutine
