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
subroutine te0332(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=*) :: option, nomte
!     FONCTION REALISEE :
!
!         CALCUL DU TAUX DE CROISSANCE DE CAVITES SELON UNE LOI DE
!         RICE ET TRACEY EN COMPORTEMENT NON-LINEAIRE.
!         ELEMENTS ISOPARAMETRIQUES 2D.
!
!         OPTION : 'RICE_TRACEY'
!
! ENTREE  --->  OPTION : NOM DE L'OPTION DE CALCUL
!         --->  NOMTE  : NOM DU TYPE D'ELEMENT
!
!     ------------------------------------------------------------------
!
!
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: optcal(12), rela_name
    real(kind=8) :: sig(6), triax, volu, rsr0, numema, depseq
    real(kind=8) :: poids, r, volume, dvol, sigm, sigeq
    real(kind=8) :: dfdx(9), dfdy(9)
    real(kind=8) :: cong(4), varigp, varigm, sdrsrp, sdrsrm, crois
    integer(kind=8) :: nno, kp, npg, k, iritra, jtab(7)
    integer(kind=8) :: issopt, ima, nbvari, ipopp, ndim, nnos, jgano
    integer(kind=8) :: ipoids, ivf, idfde, ii, iret, iadzi, ivarmg, iazk24
    integer(kind=8) :: igeom, icong, ivarpg, isdrmr, isdrpr, kq
    aster_logical :: laxi
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!     RECUPERATION DU NUMERO DE LA MAILLE :
!     -------------------------------------
    call tecael(iadzi, iazk24, noms=0)
    ima = zi(iadzi)
    numema = dble(ima)
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    poids = 0.d0
    triax = 0.d0
    rsr0 = 0.d0
    volu = 0.d0
    volume = 0.d0
    dvol = 0.d0
    depseq = 0.d0
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTPR', 'L', icong)
    call jevech('PVARIMR', 'L', ivarmg)
    call jevech('PVARIPR', 'L', ivarpg)
    call jevech('PSDRMR', 'L', isdrmr)
    call jevech('PSOUSOP', 'L', issopt)
    call tecach('OOO', 'PVARIPR', 'L', iret, nval=7, &
                itab=jtab)
    nbvari = max(jtab(6), 1)*jtab(7)
    call jevech('PCOMPOR', 'L', vk16=compor)
    rela_name = compor(RELA_NAME)
!
    if ((rela_name .eq. 'VMIS_ISOT_TRAC') .or. (rela_name .eq. 'VMIS_ISOT_LINE') .or. &
        (rela_name .eq. 'LEMAITRE') .or. (rela_name .eq. 'VMIS_ECMI_TRAC') .or. &
        (rela_name .eq. 'VMIS_ECMI_LINE') .or. (rela_name .eq. 'VISC_CIN1_CHAB') .or. &
        (rela_name .eq. 'VISC_CIN2_CHAB')) then
        ipopp = 1
    else
        call utmess('F', 'ELEMENTS3_74', sk=rela_name)
    end if
!
    call jevech('PRICTRA', 'E', iritra)
    call jevech('PSDRPR', 'E', isdrpr)
!
    optcal(1) = zk24(issopt) (1:16)
    optcal(2) = zk24(issopt) (17:19)
!
    do ii = 1, 4
        cong(ii) = 0.d0
    end do
    varigm = 0.d0
    varigp = 0.d0
!
!     --- BOUCLE SUR POINTS DE GAUSS SUIVANT OPTIONS DE CALCUL ---
!
    if ((optcal(1) .eq. 'SIGM_ELMOY') .and. (optcal(2) .eq. 'NON')) then
        do kp = 1, npg
            k = (kp-1)*nno
            r = 0.d0
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            if (laxi) then
                do ii = 1, nno
                    r = r+zr(igeom+2*ii-2)*zr(ivf+k+ii-1)
                end do
                poids = poids*r
            end if
            dvol = poids
            volume = volume+dvol
            do ii = 1, 4
                cong(ii) = cong(ii)+dvol*zr(icong+4*kp+ii-5)
            end do
            varigm = varigm+dvol*zr(ivarmg+nbvari*(kp-1)+ipopp-1)
            varigp = varigp+dvol*zr(ivarpg+nbvari*(kp-1)+ipopp-1)
        end do
!        ------- SIGXX MOYENNEE SUR L'ELEMENT -----------
        sig(1) = cong(1)/volume
!        ------- SIGYY MOYENNEE SUR L'ELEMENT -----------
        sig(2) = cong(2)/volume
!        ------- SIGZZ MOYENNEE SUR L'ELEMENT -----------
        sig(3) = cong(3)/volume
!        ------- SIGXY MOYENNEE SUR L'ELEMENT -----------
        sig(4) = cong(4)/volume
!        ------- EPSPEQ MOYENNEE SUR L'ELEMENT ----------
        varigm = varigm/volume
        varigp = varigp/volume
!
        sigm = (sig(1)+sig(2)+sig(3))/3.d0
        sigeq = ( &
                sig(1)-sigm)*(sig(1)-sigm)+(sig(2)-sigm)*(sig(2)-sigm)+(sig(3)-sigm)*(sig(3)-si&
                &gm)+2.d0*sig(4)*sig(4 &
                )
        sigeq = sqrt(1.5d0*sigeq)
        triax = sigm/sigeq
        volu = volume
        depseq = varigp-varigm
        do kq = 1, npg
            zr(isdrpr+kq-1) = zr(isdrmr+kq-1)
        end do
!
!
    elseif ((optcal(1) .eq. 'SIGM_ELGA') .and. (optcal(2) .eq. 'OUI')) &
        then
        do kp = 1, npg
            r = 0.d0
            sdrsrm = zr(isdrmr+kp-1)
            k = (kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            if (laxi) then
                do ii = 1, nno
                    r = r+zr(igeom+2*ii-2)*zr(ivf+k+ii-1)
                end do
                poids = poids*r
            end if
            volume = poids
            do ii = 1, 4
                cong(ii) = zr(icong+4*kp+ii-5)
            end do
            varigm = zr(ivarmg+nbvari*(kp-1)+ipopp-1)
            varigp = zr(ivarpg+nbvari*(kp-1)+ipopp-1)
!
            sigm = (cong(1)+cong(2)+cong(3))/3.d0
            sigeq = ( &
                    cong(1)-sigm)*(cong(1)-sigm)+(cong(2)-sigm)*(cong(2)-sigm)+(cong(3)-sigm)*(c&
                    &ong(3)-sigm)+2.d0*cong(4)*cong(4 &
                    )
            sigeq = sqrt(1.5d0*sigeq)
            triax = sigm/sigeq
            volu = volume
            depseq = varigp-varigm
            sdrsrp = sdrsrm+0.283d0*sign(1.d0, triax)*exp(1.5d0*abs( &
                                                          triax))*depseq
            zr(isdrpr+kp-1) = sdrsrp
            crois = exp(sdrsrp)
            if (crois .gt. rsr0) then
                rsr0 = crois
                volu = volume
            end if
!
        end do
!
!
    elseif ((optcal(1) .eq. 'SIGM_ELMOY') .and. (optcal(2) .eq. 'OUI')) &
        then
        do kp = 1, npg
            r = 0.d0
            k = (kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            if (laxi) then
                do ii = 1, nno
                    r = r+zr(igeom+2*ii-2)*zr(ivf+k+ii-1)
                end do
                poids = poids*r
            end if
            dvol = poids
            volume = volume+dvol
            do ii = 1, 4
                cong(ii) = cong(ii)+zr(icong+4*kp+ii-5)*dvol
            end do
            varigm = varigm+dvol*zr(ivarmg+nbvari*(kp-1)+ipopp-1)
            varigp = varigp+dvol*zr(ivarpg+nbvari*(kp-1)+ipopp-1)
        end do
!        ------- SIGXX MOYENNEE SUR L'ELEMENT -----------
        sig(1) = cong(1)/volume
!        ------- SIGYY MOYENNEE SUR L'ELEMENT -----------
        sig(2) = cong(2)/volume
!        ------- SIGZZ MOYENNEE SUR L'ELEMENT -----------
        sig(3) = cong(3)/volume
!        ------- SIGXY MOYENNEE SUR L'ELEMENT -----------
        sig(4) = cong(4)/volume
!        ------- EPSPEQ MOYENNEE SUR L'ELEMENT ----------
        varigm = varigm/volume
        varigp = varigp/volume
!
        sigm = (sig(1)+sig(2)+sig(3))/3.d0
        sigeq = ( &
                sig(1)-sigm)*(sig(1)-sigm)+(sig(2)-sigm)*(sig(2)-sigm)+(sig(3)-sigm)*(sig(3)-si&
                &gm)+2.d0*sig(4)*sig(4 &
                )
        sigeq = sqrt(1.5d0*sigeq)
        triax = sigm/sigeq
        volu = volume
        depseq = varigp-varigm
        sdrsrm = zr(isdrmr)
        sdrsrp = sdrsrm+0.283d0*sign(1.d0, triax)*exp(1.5d0*abs(triax))*depseq
        do kq = 1, npg
            zr(isdrpr+kq-1) = sdrsrp
        end do
        rsr0 = exp(sdrsrp)
!
!
    elseif ((optcal(1) .eq. 'SIGM_ELGA') .and. (optcal(2) .eq. 'NON')) &
        then
        do kp = 1, npg
            k = (kp-1)*nno
            r = 0.d0
            do ii = 1, 4
                cong(ii) = zr(icong+(4*kp)-5+ii)
            end do
            varigm = zr(ivarmg+nbvari*(kp-1)+ipopp-1)
            varigp = zr(ivarpg+nbvari*(kp-1)+ipopp-1)
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            if (laxi) then
                do ii = 1, nno
                    r = r+zr(igeom+2*ii-2)*zr(ivf+k+ii-1)
                end do
                poids = poids*r
            end if
            dvol = poids
            volume = volume+dvol
            sigm = (cong(1)+cong(2)+cong(3))/3.d0
            sigeq = ( &
                    cong(1)-sigm)*(cong(1)-sigm)+(cong(2)-sigm)*(cong(2)-sigm)+(cong(3)-sigm)*(c&
                    &ong(3)-sigm)+2.d0*cong(4)*cong(4 &
                    )
            sigeq = sqrt(1.5d0*sigeq)
            triax = triax+dvol*sigm/sigeq
            depseq = depseq+(varigp-varigm)*dvol
!
        end do
!
        triax = triax/volume
        depseq = depseq/volume
        volu = volume
        do kq = 1, npg
            zr(isdrpr+kq-1) = zr(isdrmr+kq-1)
        end do
!
!
    else
!       OPTION DE CALCUL NON VALIDE
        ASSERT(.false.)
    end if
!
!
    zr(iritra) = triax
    zr(iritra+1) = rsr0
    zr(iritra+2) = volu
    zr(iritra+3) = numema
    zr(iritra+4) = depseq
!
!     DESTRUCTION DES OBJETS CREES DANS LA BASE
!
end subroutine
