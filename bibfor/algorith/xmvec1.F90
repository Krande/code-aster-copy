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
!
subroutine xmvec1(ndim, jnne, ndeple, nnc, jnnm,&
                  hpg, ffc, ffe, ffm, jacobi,&
                  dlagrc, coefcr, coefcp, lpenac, jeu,&
                  norm, nsinge, nsingm, fk_escl, fk_mait,&
                  jddle, jddlm, nfhe, nfhm, lmulti,&
                  heavno, heavn, heavfa, vtmp)
!
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "asterfort/indent.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xplma2.h"
    integer :: ndim, jnne(3), jnnm(3), nnc
    integer :: nsinge, nsingm
    real(kind=8) :: hpg, ffc(9), jacobi, ffe(20), ffm(20)
    real(kind=8) :: dlagrc, jeu, norm(3), coefcr, coefcp
    real(kind=8) :: fk_escl(27, 3, 3), fk_mait(27, 3, 3)
    real(kind=8) :: vtmp(336)
    integer :: ndeple, jddle(2), jddlm(2)
    integer :: nfhe, nfhm, heavno(8), heavfa(*), heavn(*)
    aster_logical :: lpenac, lmulti
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE XFEMGG - CALCUL ELEM.)
!
! VECTEUR SECOND MEMBRE SI CONTACT AVEC COMPLIANCE (XFEM)
!
! ----------------------------------------------------------------------
! ROUTINE SPECIFIQUE A L'APPROCHE <<GRANDS GLISSEMENTS AVEC XFEM>>,
! TRAVAIL EFFECTUE EN COLLABORATION AVEC I.F.P.
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DU PROBLEME
! IN  NNE    : NOMBRE DE NOEUDS DE LA MAILLE ESCLAVE
! IN  NNES   : NOMBRE DE NOEUDS SOMMETS DE LA MAILLE ESCLAVE
! IN  NNM    : NOMBRE DE NOEUDS DE LA MAILLE MAITRE
! IN  HPG    : POIDS DU POINT INTEGRATION DU POINT DE CONTACT
! IN  FFC    : FONCTIONS DE FORME DU POINT DE CONTACT DANS ELC
! IN  FFE    : FONCTIONS DE FORME DU POINT DE CONTACT DANS ESC
! IN  FFM    : FONCTIONS DE FORME DE LA PROJECTION DU PTC DANS MAIT
! IN  JACOBI : JACOBIEN DE LA MAILLE AU POINT DE CONTACT
! IN  DLAGRC : LAGRANGE DE CONTACT AU POINT D'INT??GRATION
! IN  COEFCA : COEF_REGU_CONT
! IN  JEU    : VALEUR DU JEU
! IN  NORM   : VALEUR DE LA NORMALE AU POINT DE CONTACT
! IN  NSINGE  : NOMBRE DE FONCTION SINGULIERE ESCLAVE
! IN  NSINGM  : NOMBRE DE FONCTION SINGULIERE MAITRE
! I/O VTMP   : VECTEUR SECOND MEMBRE ELEMENTAIRE DE CONTACT/FROTTEMENT
! ----------------------------------------------------------------------
    integer :: i, j, ii, pl, iin, nddle
    integer :: nne, nnes, nnem, nnm, nnms, ddles, ddlem, ddlms, ddlmm
    integer :: ifh, iddl, hea_fa(2), alp
    real(kind=8) :: vv, iescl(6), imait(6)
! ----------------------------------------------------------------------
!
!
! --- INITIALISATION
!
    iescl(1) = 1
    iescl(2) =-1
    imait(1) = 1
    imait(2) = 1
!    DEFINITION A LA MAIN DE LA TOPOLOGIE DE SOUS-DOMAINE PAR FACETTE (SI NFISS=1)
    if (.not.lmulti) then
        hea_fa(1)=xcalc_code(1,he_inte=[-1])
        hea_fa(2)=xcalc_code(1,he_inte=[+1])
    endif
!
! --------------------- CALCUL DE [L1_CONT]-----------------------------
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
    if (nnm .ne. 0) then
!
        do j = 1, ndim
            do i = 1, ndeple
                call indent(i, ddles, ddlem, nnes, iin)
                if (lpenac) then
                    vv = hpg*jacobi*dlagrc* norm(j)
                else
                    vv = hpg*jacobi*(dlagrc-coefcr*jeu) *norm( j)
                endif
                if (lmulti) then
                    do ifh = 1, nfhe
                        iescl(1+ifh)=xcalc_heav(heavn(nfhe*(i-1)+ifh),&
                                                heavfa(1),&
                                                heavn(nfhe*nne+nfhm*nnm+i))
                    end do
                else
                    iescl(2)=xcalc_heav(heavn(i),&
                                            hea_fa(1),&
                                            heavn(nfhe*nne+nfhm*nnm+i))
                endif
                do iddl = 1, 1+nfhe
                    ii = iin + (iddl-1)*ndim + j
                    vtmp(ii) = -iescl(iddl)*vv* ffe(i)
                end do
                do alp = 1, ndim*nsinge
                    ii = iin + (1+nfhe+nsinge-1)*ndim + alp
                    vtmp(ii) = vtmp(ii)+fk_escl(i,alp,j)*vv
                enddo
            end do
            do i = 1, nnm
                call indent(i, ddlms, ddlmm, nnms, iin)
                iin = iin + nddle
                if (lpenac) then
                    vv = hpg*jacobi*dlagrc* norm(j)
                else
                    vv = hpg*jacobi*(dlagrc-coefcr*jeu) *norm( j)
                endif
                if (lmulti) then
                    do ifh = 1, nfhm
                        imait(1+ifh)=xcalc_heav(heavn(nfhe*nne+nfhm*(i-1)+ifh),&
                                                heavfa(2),&
                                                heavn((1+nfhe)*nne+nfhm*nnm+i))
                    end do
                else
                    imait(2)=xcalc_heav(heavn(nne+i),&
                                            hea_fa(2),&
                                            heavn((1+nfhe)*nne+nfhm*nnm+i))
                endif
                do iddl = 1, 1+nfhm
                    ii = iin + (iddl-1)*ndim + j
                    vtmp(ii) = imait(iddl)*vv*ffm(i)
                end do
                do alp = 1, ndim*nsingm
                    ii = iin + (1+nfhm+nsingm-1)*ndim + alp
                    vtmp(ii) = vtmp(ii)+fk_mait(i,alp,j)*vv
                enddo
            end do
        end do
    else
        do j = 1, ndim
            do i = 1, ndeple
! --- BLOCS ES,SI
                if (lpenac) then
                    vv = hpg*jacobi*dlagrc* ffe(i)*norm(j)
                else
                    vv = hpg*jacobi*(dlagrc-coefcr*jeu) * ffe(i)*norm( j)
                endif
                call indent(i, ddles, ddlem, nnes, iin)
                do alp = 1, ndim*nsinge
                    ii = iin + alp
                    vtmp(ii) = vtmp(ii)+fk_escl(i,alp,j)*vv
                enddo
            end do
        end do
    endif
!
! --------------------- CALCUL DE [L2]----------------------------------
!
    do i = 1, nnc
        call xplma2(ndim, nne, nnes, ddles, i,&
                    nfhe, pl)
        if (lmulti) pl = pl + (heavno(i)-1)*ndim
        if (lpenac) then
            vtmp(pl) = -hpg*jacobi*(dlagrc/coefcp+jeu) *ffc(i)
        else
            vtmp(pl) = -hpg*jacobi*jeu*ffc(i)
        endif
    end do
!
end subroutine
