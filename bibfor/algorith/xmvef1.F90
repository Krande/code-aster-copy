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

subroutine xmvef1(ndim, jnne, jnnm, ndeple, nnc,&
                  hpg, ffc, ffe,&
                  ffm, jacobi, dlagrc, dlagrf,&
                  coeffr, lpenaf, coefff, tau1,&
                  tau2, rese, mproj, coefcr,&
                  jeu, nsinge, nsingm, fk_escl, fk_mait,&
                  nvit, nconta, jddle, jddlm,&
                  nfhe, nfhm, lmulti, heavn, heavfa,&
                  vtmp)
!
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "asterfort/indent.h"
#include "asterfort/xplma2.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalc_code.h"
    integer :: ndim, nnc, jnne(3), jnnm(3)
    integer :: nsinge, nsingm, nvit, jddle(2), jddlm(2), nfhe
    integer :: nfhm, heavn(*), heavfa(*)
    real(kind=8) :: hpg, ffc(9), ffe(20), ffm(20), jacobi
    real(kind=8) :: dlagrc, dlagrf(2), jeu
    real(kind=8) :: coefff, coeffr, coefcr
    real(kind=8) :: fk_escl(27,3,3), fk_mait(27,3,3)
    real(kind=8) :: tau1(3), tau2(3), rese(3), mproj(3, 3), vtmp(336)
    integer :: nconta, ndeple
    aster_logical :: lpenaf, lmulti
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE XFEMGG - CALCUL ELEM.)
!
! CALCUL DU SECOND MEMBRE POUR LE FROTTEMENT
! CAS AVEC CONTACT
!
! ----------------------------------------------------------------------
! ROUTINE SPECIFIQUE A L'APPROCHE <<GRANDS GLISSEMENTS AVEC XFEM>>,
! TRAVAIL EFFECTUE EN COLLABORATION AVEC I.F.P.
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DU PROBLEME
! IN  NNE    : NOMBRE DE NOEUDS DE LA MAILLE ESCLAVE
! IN  NNES   : NOMBRE DE NOEUDS SOMMETS DE LA MAILLE ESCLAVE
! IN  NNC    : NOMBRE DE NOEUDS DE CONTACT
! IN  NNM    : NOMBRE DE NOEUDS DE LA MAILLE MAITRE
! IN  NFAES  : NUMERO DE LA FACETTE DE CONTACT ESCLAVE
! IN  CFACE  : MATRICE DE CONECTIVITE DES FACETTES DE CONTACT
! IN  HPG    : POIDS DU POINT INTEGRATION DU POINT DE CONTACT
! IN  FFC    : FONCTIONS DE FORME DU POINT DE CONTACT DANS ELC
! IN  FFE    : FONCTIONS DE FORME DU POINT DE CONTACT DANS ESC
! IN  FFM    : FONCTIONS DE FORME DE LA PROJECTION DU PTC DANS MAIT
! IN  JACOBI : JACOBIEN DE LA MAILLE AU POINT DE CONTACT
! IN  JPCAI  : POINTEUR VERS LE VECT DES ARRETES ESCLAVES INTERSECTEES
! IN  COEFFA : COEF_REGU_FROT
! IN  COEFFF : COEFFICIENT DE FROTTEMENT DE COULOMB
! IN  TAU1   : PREMIERE TANGENTE
! IN  TAU2   : SECONDE TANGENTE
! IN  RESE   : PROJECTION DE LA BOULE UNITE POUR LE FROTTEMENT
! IN  MPROJ  : MATRICE DE L'OPERATEUR DE PROJECTION
! IN  DLAGRF : LAGRANGES DE FROTTEMENT AU POINT D'INT??GRATION
! IN  TYPMAI : NOM DE LA MAILLE ESCLAVE D'ORIGINE (QUADRATIQUE)
! IN  NSINGE : NOMBRE DE FONCTION SINGULIERE ESCLAVE
! IN  NSINGM : NOMBRE DE FONCTION SINGULIERE MAITRE
! IN  NVIT   : POINT VITAL OU PAS
! IN  INADH  : POINT ADHERENT OU PAS
! I/O VTMP   : VECTEUR SECOND MEMBRE ELEMENTAIRE DE CONTACT/FROTTEMENT
! ----------------------------------------------------------------------
    integer :: i, j, k, ii, pli, iin, nddle, hea_fa(2), alp
    integer :: nne, nnes, nnem, nnm, nnms, ddles, ddlem, ddlms, ddlmm
    real(kind=8) :: vectt(3), tt(2), vv, t, iescl(3), imait(3)
! ----------------------------------------------------------------------
!
! --- INITIALISATIONS
!
!  CETTE ROUTINE N AUTORISE QU UNE FACETTE MONOFISSUREE
!  ON DIMENSIONNE LES CHAMPS DE SIGNES SELON CETTE HYPOTHESE
    iescl(1) = 1.d0
    iescl(2) = -1.d0
    imait(1) = 1.d0
    imait(2) = 1.d0
    if (.not.lmulti) then
      hea_fa(1)=xcalc_code(1,he_inte=[-1])
      hea_fa(2)=xcalc_code(1,he_inte=[+1])
    endif
!
    nne=jnne(1)
    nnes=jnne(2)
    nnem=jnne(3)
    nnm=jnnm(1)
    nnms=jnnm(2)
    ddles=jddle(1)
    ddlem=jddle(2)
    ddlms=jddlm(1)
    ddlmm=jddlm(2)
    nddle = ddles*nnes+ddlem*nnem
!
    vectt(:) = 0.d0
    tt(:) = 0.d0
!
! --- CALCUL DE RESE.C(*,I)
!
    do i = 1, ndim
        do k = 1, ndim
            vectt(i) = rese(k)*mproj(k,i) + vectt(i)
        end do
    end do
!
! --- CALCUL DE T.(T-P)
!
    do i = 1, ndim
        t = dlagrf(1)*tau1(i)+dlagrf(2)*tau2(i)-rese(i)
        tt(1)= t*tau1(i)+tt(1)
        if (ndim .eq. 3) tt(2)= t*tau2(i)+tt(2)
    end do
!
! --------------------- CALCUL DE [L1_FROT]-----------------------------
!
    if (nnm .ne. 0) then
!
        do j = 1, ndim
            do i = 1, ndeple
                if (lmulti) then
                    iescl(2)=xcalc_heav(heavn(nfhe*(i-1)+1),&
                                            heavfa(1),&
                                            heavn(nfhe*nne+nfhm*nnm+i))
                else
                    iescl(2)=xcalc_heav(heavn(i),&
                                            hea_fa(1),&
                                            heavn(nfhe*nne+nfhm*nnm+i))
                endif
! --- BLOCS ES,CL ; ES,EN ; (ES,SI)
                if (nconta .eq. 3 .and. ndim .eq. 3) then
                    vv = jacobi*hpg*coefff*(dlagrc-coefcr*jeu)*vectt( j)
                else
                    vv = jacobi*hpg*coefff*dlagrc*vectt(j)
                endif
                call indent(i, ddles, ddlem, nnes, iin)
                ii = iin + j
                vtmp(ii) = vtmp(ii)-vv*ffe(i)
                ii = ii + ndim
                vtmp(ii) = vtmp(ii)-vv*iescl(2)*ffe(i)
                do alp = 1, nsinge*ndim
                    ii = iin + 2*ndim + alp
                    vtmp(ii) = vtmp(ii)+fk_escl(i,alp,j)* vv
                end do
            end do
            do i = 1, nnm
                if (lmulti) then
                    imait(2)=xcalc_heav(heavn(nfhe*nne+nfhm*(i-1)+1),&
                                            heavfa(2),&
                                            heavn((1+nfhe)*nne+nfhm*nnm+i))
                else
                    imait(2)=xcalc_heav(heavn(nne+i),&
                                        hea_fa(2),&
                                        heavn((1+nfhe)*nne+nfhm*nnm+i))
                endif
                if (nconta .eq. 3 .and. ndim .eq. 3) then
                    vv = jacobi*hpg*coefff* (dlagrc-coefcr*jeu)*vectt( j)
                else
                    vv = jacobi*hpg*coefff* dlagrc*vectt(j)
                endif
                call indent(i, ddlms, ddlmm, nnms, iin)
                ii = nddle + iin + j
                vtmp(ii) = vtmp(ii)+vv*ffm(i)
                ii = ii + ndim
                vtmp(ii) = vtmp(ii)+vv*imait(2)*ffm(i)
                do alp = 1, nsingm*ndim
                    ii = ii + 2*ndim + alp
                    vtmp(ii) = vtmp(ii)+vv*fk_mait(i,alp,j)
                end do
            end do
        end do
    else
!
        do j = 1, ndim
            do i = 1, ndeple
! --- BLOCS ES,SI
                if (nconta .eq. 3 .and. ndim .eq. 3) then
                    vv = jacobi*hpg*coefff* (dlagrc-coefcr*jeu)*vectt( j)
                else
                    vv = jacobi*hpg*coefff* dlagrc*vectt(j)
                endif
                call indent(i, ddles, ddlem, nnes, iin)
                do alp=1,ndim*nsinge
                  ii = iin + alp
                  vtmp(ii) = vtmp(ii)+fk_escl(i,alp,j)* vv
                enddo
            end do
        end do
    endif
!
! --------------------- CALCUL DE [L3]----------------------------------
!
    if (nvit .eq. 1) then
        do i = 1, nnc
            call xplma2(ndim, nne, nnes, ddles, i,&
                        nfhe, pli)
            do j = 1, ndim-1
                ii = pli+j
                if (lpenaf) then
                    vtmp(ii) = vtmp(ii)+jacobi*hpg*tt(j)*ffc(i)
                else
                    if (nconta .eq. 3 .and. ndim .eq. 3) then
                        vtmp(ii) = vtmp(ii)+jacobi*hpg*tt(j)*ffc(i)/coeffr
                    else
                        vtmp(ii) = vtmp(ii)+jacobi*hpg*tt(j)*ffc(i)*coefff* dlagrc/coeffr
                    endif
                endif
            end do
        end do
    endif
!
end subroutine
