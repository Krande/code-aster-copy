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
subroutine xmmjec(ndim, jnnm, jnne, ndeple, nsinge,&
                  nsingm, ffe, ffm, norm, jgeom,&
                  jdepde, fk_escl, fk_mait, jddle, jddlm,&
                  nfhe, nfhm, lmulti, heavn, heavfa,&
                  jeuca)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/indent.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
!
    integer :: ndim
    real(kind=8) :: norm(3)
    real(kind=8) :: ffe(20), ffm(20)
    real(kind=8) :: jeuca
    integer :: jgeom, jdepde
    integer :: jnnm(3), jnne(3), jddle(2), jddlm(2)
    integer :: nsinge, nsingm, nfhe, nfhm, heavfa(*), heavn(*)
    real(kind=8) :: fk_escl(27, 3, 3), fk_mait(27, 3, 3)
    aster_logical :: lmulti
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (METHODE XFEM-GG - TE)
!
! CALCUL DU JEU
!
! ----------------------------------------------------------------------
!
!
! IN  NDIM   : DIMENSION DU PROBLEME
! IN  NDDL   : NOMBRE TOTAL DE DEGRES DE LIBERTE DE LA MAILLE DE CONTACT
! IN  NNE    : NOMBRE DE NOEUDS DE LA MAILLE ESCLAVE
! IN  NNM    : NOMBRE DE NOEUDS DE LA MAILLE MAITRE
! IN  NNC    : NOMBRE DE NOEUDS DE LA MAILLE DE CONTACT
! IN  NNES   : NOMBRE DE NOEUDS SOMMETS DE LA MAILLE ESCLAVE
! IN  NSINGE : NOMBRE DE FONCTIONS SINGULIERE ESCLAVES
! IN  NSINGM : NOMBRE DE FONCTIONS SINGULIERE MAIT RES
! IN  DDLES : NOMBRE DE DDLS D'UN NOEUD SOMMET ESCLAVE
! IN  NORM   : VALEUR DE LA NORMALE
! IN  JGEOM  : POINTEUR JEVEUX SUR GEOMETRIE INITIALE
! IN  JDEPDE : POINTEUR JEVEUX POUR DEPDEL
! OUT JEUCA  : VALEUR DU JEU POUR LES SECONDS MEMBRES DE CONTACT
!
!
!
!
    integer :: idim, inom, inoes, in, nddle
    integer :: ndeple, nne, nnes, nnem, nnm, nnms, ddles, ddlem, ddlms, ddlmm
    integer :: ifh, iddl, hea_fa(2), alp
    real(kind=8) :: pose(3), posm(3), iescl(6), imait(6), pos
    integer :: pl
!
! ----------------------------------------------------------------------
!
!
!
! --- INNITIALISATION
!
    iescl(1) = 1
    iescl(2) = -1
    imait(1) = 1
    imait(2) = 1
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
    jeuca = 0.d0
    pose(:) = 0.d0
    posm(:) = 0.d0
!    DEFINITION A LA MAIN DE LA TOPOLOGIE DE SOUS-DOMAINE PAR FACETTE (SI NFISS=1)
    if (.not.lmulti) then
        hea_fa(1)=xcalc_code(1,he_inte=[-1])
        hea_fa(2)=xcalc_code(1,he_inte=[+1])
    endif
!
! --- CALCUL DE LA POSITION COURANTE DU POINT ESCLAVE
!
    do idim = 1, ndim
        do inoes = 1, ndeple
            call indent(inoes, ddles, ddlem, nnes, in)
            if (nnm .ne. 0) then
                if (lmulti) then
                    do ifh = 1, nfhe
                        iescl(1+ifh)=xcalc_heav(heavn(nfhe*(inoes-1)+ifh),&
                                                heavfa(1),&
                                                heavn(nfhe*nne+nfhm*nnm+inoes))
                    end do
                else
                    iescl(2)=xcalc_heav(heavn(inoes),&
                                            hea_fa(1),&
                                            heavn(nfhe*nne+nfhm*nnm+inoes))
                endif
                pos = zr(jgeom-1+ndim*(inoes-1)+idim)
                do iddl = 1, 1+nfhe
                    pl = in + (iddl-1)*ndim + idim
                    pos = pos + iescl(iddl)*zr(jdepde-1+pl)
                end do
                pose(idim) = pose(idim) + pos*ffe(inoes)
                do alp = 1, ndim*nsinge
                    pl = in + (1+nfhe+nsinge-1)*ndim + alp
                    pose(idim) = pose(idim) - fk_escl(inoes,alp,idim)* zr(jdepde-1+pl)
                end do
!
            else
                do alp = 1, ndim*nsinge
                    pl = in + alp
                    pose(idim) = pose(idim) - fk_escl(inoes,alp,idim)* zr( jdepde-1+pl)
                enddo
            endif
        end do
    end do
!
! --- CALCUL DE LA POSITION COURANTE DU POINT MAITRE
!
    do idim = 1, ndim
        do inom = 1, nnm
            call indent(inom, ddlms, ddlmm, nnms, in)
            in = in + nddle
            if (lmulti) then
                do ifh = 1, nfhm
                    imait(1+ifh)=xcalc_heav(heavn(nfhe*nne+nfhm*(inom-1)+ifh),&
                                            heavfa(2),&
                                            heavn((1+nfhe)*nne+nfhm*nnm+inom))
                end do
            else
                imait(2)=xcalc_heav(heavn(nne+inom),&
                                        hea_fa(2),&
                                        heavn((1+nfhe)*nne+nfhm*nnm+inom))
            endif
            pos = zr(jgeom-1+nne*ndim+(inom-1)*ndim+idim)
            do iddl = 1, 1+nfhm
                pl = in + (iddl-1)*ndim + idim
                pos = pos + imait(iddl)*zr(jdepde-1+pl)
            end do
            posm(idim) = posm(idim) + pos*ffm(inom)
            do alp = 1, ndim*nsingm
                pl = in + (1+nfhm+nsingm-1)*ndim + alp
                posm(idim) = posm(idim) + fk_mait(inom,alp,idim)* zr( jdepde-1+pl)
            enddo
        end do
    end do
!
! --- CALCUL DU JEU
!
    do idim = 1, ndim
        jeuca = jeuca + norm(idim)*(pose(idim)-posm(idim))
    end do
!
end subroutine
