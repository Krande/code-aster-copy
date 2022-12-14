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
subroutine xtdepm(ndim, jnnm, jnne, ndeple, nsinge,&
                  nsingm, ffe, ffm, jdepde, fk_escl,&
                  fk_mait, jddle, jddlm, nfhe, nfhm,&
                  lmulti, heavn, heavfa, ddeple, ddeplm)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/indent.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
    integer :: ndim, jnnm(3), jnne(3), nfhe, nfhm
    integer :: nsinge, nsingm, heavn(*), heavfa(*)
    integer :: jdepde, ndeple, jddle(2), jddlm(2)
    real(kind=8) :: ffm(20), ffe(20)
    real(kind=8) :: ddeple(3), ddeplm(3)
    real(kind=8) :: fk_escl(27, 3, 3), fk_mait(27, 3, 3)
    aster_logical :: lmulti
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE XFEMGG - UTILITAIRE)
!
! CALCUL DES INCREMENTS - DEPLACEMENTS
!
! ----------------------------------------------------------------------
!
!
! IN  NDIM   : DIMENSION DU PROBLEME
! IN  NNE    : NOMBRE DE NOEUDS DE LA MAILLE ESCLAVE
! IN  NNM    : NOMBRE DE NOEUDS DE LA MAILLE MAITRE
! IN  NNES   : NOMBRE DE NOEUDS SOMMETS DE LA MAILLE ESCLAVE
! IN  NSINGE : NOMBRE DE FONCTIONS SINGULIERE ESCLAVES
! IN  NSINGM : NOMBRE DE FONCTIONS SINGULIERE MAIT RES
! IN  DDLES : NOMBRE DE DDLS D'UN NOEUD SOMMET ESCLAVE
! IN  JDEPDE : POINTEUR JEVEUX POUR DEPDEL
! IN  FFE    : FONCTIONS DE FORMES ESCLAVE
! IN  FFM    : FONCTIONS DE FORMES MAITRE
! OUT DDEPLE : INCREMENT DEPDEL DU DEPL. DU POINT DE CONTACT
! OUT DDEPLM : INCREMENT DEPDEL DU DEPL. DU PROJETE DU POINT DE CONTACT
!
!
!
!
    integer :: idim, inoe, inom, pl, in, iddl, hea_fa(2), nddle, nnem
    integer :: nne, nnes, nnm, nnms, ddles, ddlem, ddlms, ddlmm, ifh, alp
    real(kind=8) :: iescl(6), imait(6)
!
! ----------------------------------------------------------------------
!
!
! --- INITIALISATIONS
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
    nddle=ddles*nnes+ddlem*nnem
!
    ddeplm(:) = 0.d0
    ddeple(:) = 0.d0
    iescl(:) = 0.d0
    imait(:) = 0.d0
!
    iescl(1) = 1
    iescl(2) = -1
    imait(1) = 1
    imait(2) = 1
    if (.not.lmulti) then
        hea_fa(1)=xcalc_code(1,he_inte=[-1])
        hea_fa(2)=xcalc_code(1,he_inte=[+1])
    endif
!
!
    do idim = 1, ndim
        do inoe = 1, ndeple
            call indent(inoe, ddles, ddlem, nnes, in)
            if (nnm .ne. 0) then
                if (lmulti) then
                    do ifh = 1, nfhe
                        iescl(1+ifh)=xcalc_heav(heavn(nfhe*(inoe-1)+ifh),&
                                                heavfa(1),&
                                                heavn(nfhe*nne+nfhm*nnm+inoe))
                    end do
                else
                    iescl(2)=xcalc_heav(heavn(inoe),&
                                            hea_fa(1),&
                                            heavn(nfhe*nne+nfhm*nnm+inoe))
                endif
                do iddl = 1, 1+nfhe
                    pl = in + (iddl-1)*ndim + idim
                    ddeple(idim) = ddeple(idim)+ ffe(inoe)*iescl(iddl)*zr(jdepde-1+ pl)
                end do
            endif
            do alp = 1, ndim*nsinge
                pl = in + (1+nfhe+nsinge-1)*ndim + alp
                ddeple(idim) = ddeple(idim) -fk_escl(inoe,alp,idim)*zr(jdepde-1+pl)
            enddo
        end do
    end do
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
            do iddl = 1, 1+nfhm
                pl = in + (iddl-1)*ndim + idim
                ddeplm(idim) = ddeplm(idim) + ffm(inom)*imait(iddl)*zr(jdepde-1+pl)
            end do
            do alp = 1, ndim*nsingm
                pl = in + (1+nfhm+nsingm-1)*ndim + alp
                ddeplm(idim) = ddeplm(idim) +fk_mait(inom,alp,idim)*zr(jdepde-1+pl)
            enddo
        end do
    end do
!
!
end subroutine
