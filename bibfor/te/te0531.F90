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

subroutine te0531(option, nomte)
!
    implicit none
!
#include "jeveux.h"
!
#include "MultiFiber_type.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/verift.h"
#include "asterfort/tecach.h"
#include "asterfort/jevech.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lteatt.h"
#include "asterfort/pmfinfo.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvarc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    character(len=16)  :: nomte, option
! ======================================================================
!
!    OPTIONS  : 'EPVC_ELGA' 'EPME_ELGA' 'EPSP_ELGA'
!    ELEMENTS : DKT, GRILLE, PMF, TUYAU, BARRE
!
!    POUR ELVC_ELGA :
!    1 COMPOSANTES :
!    EPTHER= DILATATION THERMIQUE (LONGI)   : ALPHA*(T-TREF)
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................

    integer(kind=8) :: npg
    integer(kind=8) :: nbcmp, jnbspi, idefo, nbcou, icomp, imate, iret
    integer(kind=8) :: nbsp, ipg, ksp, i, idefto, icomp2, idsig, icomp3
    integer(kind=8) :: nbgf, icp, isdcom, ig, nbfig, ifib, icaba
    integer(kind=8) :: nbsec, iret1, nbpar, nbv, icodre(2), tygrfi, nbcarm, nug(10)
    real(kind=8) :: epsth, sigma(6), trsig, temp, valres(2), e, nu, c1, c2, a
    character(len=4)  :: fami
    character(len=8)  :: materi, nompar, nomres(2)
    character(len=32) :: phenom
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: lmeca, pmf, grille, tuyau, barre, coque, lplas

    nbcmp = 1
    materi = ' '
    lmeca = .false.
    lplas = .false.
    fami = 'RIGI'
!
    call elrefe_info(fami=fami, npg=npg)

    grille = lteatt('GRILLE', 'OUI')
    pmf = lteatt('TYPMOD2', 'PMF')
    tuyau = lteatt('TUYAU', 'OUI')
    barre = (nomte .eq. 'MECA_BARRE')
    if (grille) then
        coque = .false.
    else
        coque = lteatt('COQUE', 'OUI')
    end if
!
    call tecach('NNO', 'PMATERC', 'L', iret, iad=imate)
    call jevech('PDEFOPG', 'E', idefo)

    if (option .eq. 'EPME_ELGA' .or. option .eq. 'EPSP_ELGA') then
        lmeca = .true.
        call jevech('PDEFORR', 'L', idefto)
    end if
!
    if (option .eq. 'EPSP_ELGA') then
        lplas = .true.
        call jevech('PCONTRR', 'L', idsig)
!
        nbv = 2
        nomres(1) = 'E'
        nomres(2) = 'NU'
        nbpar = 1
        nompar = 'TEMP'
!
        call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
        if (phenom .ne. 'ELAS') call utmess('F', 'ELEMENTS_49', valk=[phenom, option])
    end if
!
    if (.not. grille .and. .not. barre) then
        call jevech('PNBSP_I', 'L', jnbspi)
    else
        nbsp = 1
    end if
!
    if (coque .or. tuyau) then
!
        nbcou = zi(jnbspi-1+1)
        if (lplas) then
            call utmess('A', 'ELEMENTS3_13')
        end if
        if (tuyau) then
            nbsec = zi(jnbspi-1+2)
        else
            nbsec = 3
        end if
        nbsp = nbcou*nbsec
        if (lmeca) nbcmp = 6
!
        do ipg = 1, npg
            do ksp = 1, nbsp
!
                call verift(fami, ipg, ksp, '+', zi(imate), &
                            epsth_=epsth)
!
                icomp = idefo+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                if (lmeca) then
                    icomp2 = idefto+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                    do i = 1, 2
                        zr(icomp+i) = zr(icomp2+i)-epsth
                    end do
                    do i = 3, 6
                        zr(icomp+i) = zr(icomp2+i)
                    end do
                    if (lplas) then
                        call rcvarc(' ', 'TEMP', '+', fami, ipg, &
                                    ksp, temp, iret1)
                        call rcvalb(fami, ipg, ksp, '+', zi(imate), &
                                    ' ', phenom, nbpar, nompar, [temp], &
                                    nbv, nomres, valres, icodre, 1)
!
                        e = valres(1)
                        nu = valres(2)
                        c1 = (1.d0+nu)/e
                        c2 = nu/e
!
                        icomp3 = idsig+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
!
                        do i = 1, nbcmp
                            sigma(i) = zr(icomp3+i)
                        end do
!
                        trsig = sigma(1)+sigma(2)+sigma(3)
!
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
!
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
!
        do ipg = 1, npg
            do ksp = 1, nbsp
!
                call verift(fami, ipg, ksp, '+', zi(imate), &
                            epsth_=epsth)
!
                icomp = idefo+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                if (lmeca) then
                    icomp2 = idefto+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                    zr(icomp+1) = zr(icomp2+1)-epsth
                    if (lplas) then
                        call rcvarc(' ', 'TEMP', '+', fami, ipg, &
                                    ksp, temp, iret1)
                        call rcvalb(fami, ipg, ksp, '+', zi(imate), &
                                    ' ', phenom, nbpar, nompar, [temp], &
                                    nbv, nomres, valres, icodre, 1)
!
                        e = valres(1)
                        c1 = 1.d0/e
!
                        icomp3 = idsig+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
!                       soustraction des deformations elastiques
                        zr(icomp+1) = zr(icomp+1)-c1*zr(icomp3+1)/a
                    end if
                else
                    zr(icomp+1) = epsth
                end if
!
            end do
        end do
!
    else if (pmf) then
!       Récupération des caractéristiques des fibres
        call pmfinfo(nbsp, nbgf, tygrfi, nbcarm, nug)
!
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
!
                    ksp = ksp+1
                    ASSERT(ksp .le. nbsp)

                    call verift(fami, ipg, ksp, '+', zi(imate), &
                                materi_=materi, epsth_=epsth)

                    icomp = idefo+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                    if (lmeca) then
                        icomp2 = idefto+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1
                        zr(icomp+1) = zr(icomp2+1)-epsth
                        if (lplas) then
                            call rcvarc(' ', 'TEMP', '+', fami, ipg, &
                                        ksp, temp, iret1)
                            call rcvalb(fami, ipg, ksp, '+', zi(imate), &
                                        ' ', phenom, nbpar, nompar, [temp], &
                                        nbv, nomres, valres, icodre, 1)
!
                            e = valres(1)
                            c1 = 1.d0/e
!
                            icomp3 = idsig+nbcmp*nbsp*(ipg-1)+nbcmp*(ksp-1)-1

!                           soustraction des deformations elastiques
                            zr(icomp+1) = zr(icomp+1)-c1*zr(icomp3+1)
                        end if
                    else
                        zr(icomp+1) = epsth
                    end if
!
                end do
            end do
        end do
    else
        ASSERT(.false.)
    end if
!
end subroutine
