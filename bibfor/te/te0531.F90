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
subroutine te0531(option, nomte)
!
    use pipeElem_module
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/jevech.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lteatt.h"
#include "asterfort/pmfinfo.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "jeveux.h"
#include "MultiFiber_type.h"
!
    character(len=16), intent(in) :: nomte, option
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: BARRE, DKT, DST, GRILLE_*, POU_D_*, Q4G, TUYAU
!
! Options: EPVC_ELGA, EPME_ELGA, EPSP_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: npg
    integer(kind=8) :: nbcmp, jnbspi, jvEpsiMeca, nbLayer, icomp, jvMater, iret
    integer(kind=8) :: nbsp, ipg, ksp, i, jvEpsi, icomp2, jvSigm, icomp3
    integer(kind=8) :: nbgf, icp, isdcom, ig, nbfig, ifib, icaba
    integer(kind=8) :: nbsec, tygrfi, nbcarm, nug(10)
    real(kind=8) :: epsth, sigma(6), trsig, e, nu, c1, c2, a
    character(len=4)  :: fami
    character(len=8)  :: materi
    character(len=16) :: elasKeyword
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: elasID
    aster_logical :: lmeca, pmf, grille, lPipe, barre, coque, lplas
!
! --------------------------------------------------------------------------------------------------
!
    nbcmp = 1
    materi = ' '
    lmeca = .false.
    lplas = .false.
    fami = 'RIGI'

! - Get modelisation
    call elrefe_info(fami=fami, npg=npg)
    grille = lteatt('GRILLE', 'OUI')
    pmf = lteatt('TYPMOD2', 'PMF')
    lPipe = lteatt('TUYAU', 'OUI')
    barre = (nomte .eq. 'MECA_BARRE')
    if (grille) then
        coque = .false.
    else
        coque = lteatt('COQUE', 'OUI')
    end if

! - Get material properties
    call tecach('NNO', 'PMATERC', 'L', iret, iad=jvMater)

! - Get (mechanical) strains
    call jevech('PDEFOPG', 'E', jvEpsiMeca)

! - Get (total) strains
    if (option .eq. 'EPME_ELGA' .or. option .eq. 'EPSP_ELGA') then
        lmeca = .true.
        call jevech('PDEFORR', 'L', jvEpsi)
    end if

! - Get stress and access to elastic material parameters
    if (option .eq. 'EPSP_ELGA') then
        lplas = .true.
        call jevech('PCONTRR', 'L', jvSigm)
        call get_elas_id(zi(jvMater), elasID, elasKeyword)
        if (elasKeyword .ne. 'ELAS') then
            call utmess('F', 'ELEMENTS_49', valk=[elasKeyword, option])
        end if
    end if

! - Get number of sub-points
    nbsp = 0
    if (coque) then
        call jevech('PNBSP_I', 'L', jnbspi)
        nbLayer = zi(jnbspi-1+1)
        nbsp = 3*nbLayer
        nbsec = 3
        nbsp = nbLayer*nbsec

    elseif (lPipe) then
        call pipeGetSubPoints(nspg_=nbsp)

    elseif (grille .or. barre) then
        nbsp = 1

    elseif (pmf) then

    else
        ASSERT(ASTER_FALSE)
    end if

! - Compute option
    if (coque .or. lPipe) then
        if (lmeca) then
            nbcmp = 6
        end if
        do ipg = 1, npg
            do ksp = 1, nbsp
                call verift(fami, ipg, ksp, '+', zi(jvMater), &
                            epsth_=epsth)
                icomp = jvEpsiMeca+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                if (lmeca) then
                    icomp2 = jvEpsi+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                    do i = 1, 2
                        zr(icomp+i) = zr(icomp2+i)-epsth
                    end do
                    do i = 3, 6
                        zr(icomp+i) = zr(icomp2+i)
                    end do
                    if (lplas) then
                        call get_elas_para(fami, zi(jvMater), '+', ipg, ksp, &
                                           elasID, elasKeyword, &
                                           e_=e, nu_=nu)
                        c1 = (1.d0+nu)/e
                        c2 = nu/e
                        icomp3 = jvSigm+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                        do i = 1, nbcmp
                            sigma(i) = zr(icomp3+i)
                        end do
                        trsig = sigma(1)+sigma(2)+sigma(3)
!                       soustraction des deformations elastiques
                        do i = 1, 2
                            zr(icomp+i) = zr(icomp+i)-(c1*sigma(i)-c2*trsig)
                        end do
!                       on laisse la composante 3 nulle
                        zr(icomp+4) = zr(icomp+4)-c1*sigma(4)
!                       les composantes EPXZ et EPYZ mises à zero,
!                       le calcul produirait des resultats faux car
!                       SIXZ et SIYZ sont nulles dans le champ de contrainte
                        do i = 5, 6
                            zr(icomp+i) = 0.D0
                        end do
                    end if
                else
                    zr(icomp+1) = epsth
                end if
            end do
        end do
!
    else if (grille .or. barre) then
        if (barre .and. lplas) then
            call jevech('PCAGNBA', 'L', icaba)
            a = zr(icaba)
        else
            a = 1.d0
        end if
        do ipg = 1, npg
            do ksp = 1, nbsp
                call verift(fami, ipg, ksp, '+', zi(jvMater), &
                            epsth_=epsth)
                icomp = jvEpsiMeca+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                if (lmeca) then
                    icomp2 = jvEpsi+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                    zr(icomp+1) = zr(icomp2+1)-epsth
                    if (lplas) then
                        call get_elas_para(fami, zi(jvMater), '+', ipg, ksp, &
                                           elasID, elasKeyword, &
                                           e_=e)
                        c1 = 1.d0/e
                        icomp3 = jvSigm+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
!                       soustraction des deformations elastiques
                        zr(icomp+1) = zr(icomp+1)-c1*zr(icomp3+1)/a
                    end if
                else
                    zr(icomp+1) = epsth
                end if
            end do
        end do
!
    else if (pmf) then
!       Récupération des caractéristiques des fibres
        call pmfinfo(nbsp, nbgf, tygrfi, nbcarm, nug)
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jeveuo(compor(MULTCOMP), 'L', isdcom)
        do ipg = 1, npg
            ksp = 0
            do ig = 1, nbgf
                icp = isdcom-1+(nug(ig)-1)*MULTI_FIBER_SIZEK
!               nombre de fibres de ce groupe
                read (zk24(icp+MULTI_FIBER_NBFI), '(I24)') nbfig
                materi = zk24(icp+MULTI_FIBER_MATER) (1:8)
                do ifib = 1, nbfig
                    ksp = ksp+1
                    ASSERT(ksp .le. nbsp)

                    call verift(fami, ipg, ksp, '+', zi(jvMater), &
                                materi_=materi, epsth_=epsth)
                    icomp = jvEpsiMeca+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                    if (lmeca) then
                        icomp2 = jvEpsi+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                        zr(icomp+1) = zr(icomp2+1)-epsth
                        if (lplas) then
                            call get_elas_para(fami, zi(jvMater), '+', ipg, ksp, &
                                               elasID, elasKeyword, &
                                               e_=e)
                            c1 = 1.d0/e
                            icomp3 = jvSigm+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
!                           soustraction des deformations elastiques
                            zr(icomp+1) = zr(icomp+1)-c1*zr(icomp3+1)
                        end if
                    else
                        zr(icomp+1) = epsth
                    end if
                end do
            end do
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
