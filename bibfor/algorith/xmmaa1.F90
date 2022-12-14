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

subroutine xmmaa1(ndim, jnne, ndeple, nnc, jnnm,&
                  hpg, ffc, ffe,&
                  ffm, jacobi, coefcr, coefcp,&
                  lpenac, norm, nsinge, nsingm,&
                  fk_escl, fk_mait, jddle, jddlm,&
                  nfhe, nfhm, lmulti, heavno, heavn, heavfa,&
                  mmat)
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "asterfort/indent.h"
#include "asterfort/xplma2.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalc_code.h"
    integer :: ndim, jnne(3), jnnm(3)
    integer :: nsinge, nsingm
    integer :: nfhe, nfhm, heavno(8), heavfa(*), heavn(*)
    real(kind=8) :: mmat(336, 336), norm(3)
    real(kind=8) :: hpg, ffc(8), ffe(20), ffm(20), jacobi
    real(kind=8) :: coefcr, coefcp
    real(kind=8) :: fk_escl(27,3,3), fk_mait(27,3,3)
    integer :: ndeple, nnc, jddle(2), jddlm(2)
    aster_logical :: lpenac, lmulti
!
!
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE XFEMGG - CALCUL ELEM.)
!
! CALCUL DE A ET DE AT
!
!
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! ROUTINE SPECIFIQUE A L'APPROCHE <<GRANDS GLISSEMENTS AVEC XFEM>>,
! TRAVAIL EFFECTUE EN COLLABORATION AVEC I.F.P.
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DU PROBLEME
! IN  JNNE   : MAILLE ESCL : (1) NB NDS
!                            (2) NB NDS SOMMETS
!                            (3) NB NDS MILIEU
! OUT NDEPLE : NOMBRE DE NOEUDS ESCL POSSEDANT DES DDLS DE DEPLACEMENT
! IN  NNC    : NOMBRE DE NOUEDS DE CONTACT
! IN  JNNM   : MAILLE MAIT : (1) NB NDS
!                            (2) NB NDS SOMMETS
!                            (3) NB NDS MILIEU
! IN  HPG    : POIDS DU POINT INTEGRATION DU POINT DE CONTACT
! IN  FFC    : FONCTIONS DE FORME DU POINT DE CONTACT DANS ELC
! IN  FFE    : FONCTIONS DE FORME DU POINT DE CONTACT DANS ESC
! IN  FFM    : FONCTIONS DE FORME DE LA PROJECTION DU PTC DANS MAIT
! IN  JACOBI : JACOBIEN DE LA MAILLE AU POINT DE CONTACT
! IN  COEFCA : COEF_REGU_CONT
! IN  NORM   : VALEUR DE LA NORMALE AU POINT DE CONTACT
! IN  MAILLE : NOM DE LA MAILLE ESCLAVE D'ORIGINE (QUADRATIQUE)
! IN  NSINGE : NOMBRE DE FONCTION SINGULIERE ESCLAVE
! IN  NSINGM : NOMBRE DE FONCTION SINGULIERE MAITRE
! IN  NCONTA : TYPE DE CONTACT (1=P1P1, 2=P1P1A, 3=P2P1)
! IN  JDDLE  : MAILLE ESCL : (1) DDLS D'UN NOEUD SOMMET
!                            (2) DDLS D'UN NOEUD MILIEU
! IN  JDDLM  : MAILLE MAIT : (1) DDLS D'UN NOEUD SOMMET
!                            (2) DDLS D'UN NOEUD MILIEU
! I/O MMAT   : MATRICE ELEMENTAIRE DE CONTACT/FROTTEMENT
!
! ----------------------------------------------------------------------
!
    integer :: i, j, k, l, ii, jj, pl, jjn, iin, nddle
    integer :: nne, nnes, nnm, nnms, ddles, ddlem, ddlms, ddlmm
    integer :: pli, plj, ifh, iddl, jddl, hea_fa(2), alpi, alpj
    real(kind=8) :: mm, iescl(6), jescl(6), imait(6), jmait(6)
!
! ----------------------------------------------------------------------
!
!
! --- INITIALISATION
!
    iescl(1) = 1
    iescl(2) = -1
    jescl(1) = 1
    jescl(2) = -1
    imait(1) = 1
    imait(2) = 1
    jmait(1) = 1
    jmait(2) = 1
!    DEFINITION A LA MAIN DE LA TOPOLOGIE DE SOUS-DOMAINE PAR FACETTE (SI NFISS=1)
    if (.not.lmulti) then
      hea_fa(1)=xcalc_code(1,he_inte=[-1])
      hea_fa(2)=xcalc_code(1,he_inte=[+1])
    endif
!
! --------------------- CALCUL DE [A] ----------------------------------
!
    nne=jnne(1)
    nnes=jnne(2)
    nnm=jnnm(1)
    nnms=jnnm(2)
    ddles=jddle(1)
    ddlem=jddle(2)
    ddlms=jddlm(1)
    ddlmm=jddlm(2)
    nddle = ddles*nnes+ddlem*(nne-nnes)
!
    if (nnm .ne. 0) then
        do k = 1, ndim
            do i = 1, nnc
                call xplma2(ndim, nne, nnes, ddles, i,&
                            nfhe, pl)
                if (lmulti) pl = pl + (heavno(i)-1)*ndim
                do j = 1, ndeple
                    mm = hpg*ffc(i)*jacobi*norm(k)
                    call indent(j, ddles, ddlem, nnes, jjn)
                    if (lmulti) then
                        do ifh = 1, nfhe
                            jescl(1+ifh)=xcalc_heav(heavn(nfhe*(j-1)+ifh),&
                                                    heavfa(1),&
                                                    heavn(nfhe*nne+nfhm*nnm+j))
                        end do
                    else
                        jescl(2)=xcalc_heav(heavn(j),&
                                            hea_fa(1),&
                                            heavn(nfhe*nne+nfhm*nnm+j))
                    endif
                    do jddl = 1, 1+nfhe
                        jj = jjn + (jddl-1)*ndim + k
                        mmat(pl,jj) = mmat(pl,jj)-jescl(jddl)*mm*ffe(j)
                        mmat(jj,pl) = mmat(jj,pl)-jescl(jddl)*mm*ffe(j)
                    end do
                    do alpj = 1, nsinge*ndim
                        jj = jjn + (1+nfhe+1-1)*ndim + alpj
                        mmat(pl,jj) = mmat(pl,jj)+mm*fk_escl(j,alpj,k)
                        mmat(jj,pl) = mmat(jj,pl)+mm*fk_escl(j,alpj,k)
                    end do
                end do
                do j = 1, nnm
                    mm = hpg*ffc(i)*jacobi*norm(k)
                    call indent(j, ddlms, ddlmm, nnms, jjn)
                    jjn = jjn + nddle
                    if (lmulti) then
                        do ifh = 1, nfhm
                            jmait(1+ifh)=xcalc_heav(heavn(nfhe*nne+nfhm*(j-1)+ifh),&
                                                    heavfa(2),&
                                                    heavn((1+nfhe)*nne+nfhm*nnm+j))
                        end do
                    else
                        jmait(2)=xcalc_heav(heavn(nne+j),&
                                            hea_fa(2),&
                                            heavn((1+nfhe)*nne+nfhm*nnm+j))
                    endif
                    do jddl = 1, 1+nfhm
                        jj = jjn + (jddl-1)*ndim + k
                        mmat(pl,jj) = mmat(pl,jj)+jmait(jddl)*mm*ffm(j)
                        mmat(jj,pl) = mmat(jj,pl)+jmait(jddl)*mm*ffm(j)
                    end do
                    do alpj = 1, nsingm*ndim
                        jj = jjn + (1+nfhm+1-1)*ndim + alpj
                        mmat(pl,jj) = mmat(pl,jj)+mm*fk_mait(j,alpj,k)
                        mmat(jj,pl) = mmat(jj,pl)+mm*fk_mait(j,alpj,k)
                    end do
                end do
            end do
        end do
!
! --------------------- CALCUL DE [AU]----------------------------------
!
        do k = 1, ndim
            do l = 1, ndim
                do i = 1, ndeple
                    call indent(i, ddles, ddlem, nnes, iin)
                    do j = 1, ndeple
                        call indent(j, ddles, ddlem, nnes, jjn)
                        if (lpenac) then
                            mm = 0.d0
                        else
                            mm = hpg*coefcr*norm(l)* jacobi*norm(k)
                        endif
                        if (lmulti) then
                            do ifh = 1, nfhe
                                iescl(1+ifh)=xcalc_heav(heavn(nfhe*(i-1)+ifh),&
                                                        heavfa(1),&
                                                        heavn(nfhe*nne+nfhm*nnm+i))
                                jescl(1+ifh)=xcalc_heav(heavn(nfhe*(j-1)+ifh),&
                                                        heavfa(1),&
                                                        heavn(nfhe*nne+nfhm*nnm+j))
                            end do
                        else
                                iescl(2)=xcalc_heav(heavn(i),&
                                                    hea_fa(1),&
                                                    heavn(nfhe*nne+nfhm*nnm+i))
                                jescl(2)=xcalc_heav(heavn(j),&
                                                    hea_fa(1),&
                                                    heavn(nfhe*nne+nfhm*nnm+j))
                        endif
                        do iddl = 1, 1+nfhe
                            ii = iin + (iddl-1)*ndim + l
                            do jddl = 1, 1+nfhe
                                jj = jjn + (jddl-1)*ndim + k
                                mmat(ii,jj) = mmat(ii,jj)+&
                                              iescl(iddl)*jescl(jddl)*mm*ffe(j)*ffe(i)
                            end do
                            do alpj = 1, nsinge*ndim
                                jj = jjn + (1+nfhe+1-1)*ndim + alpj
                                mmat(ii,jj) = mmat(ii,jj)+&
                                              iescl(iddl)*mm*ffe(i)*fk_escl(j,alpj,k)
                            end do
                        end do
                        do alpi = 1, nsinge*ndim
                            ii = iin + (1+nfhe+1-1)*ndim + alpi
                            do jddl = 1, 1+nfhe
                                jj = jjn + (jddl-1)*ndim + k
                                mmat(ii,jj) = mmat(ii,jj)-&
                                              jescl(jddl)*mm*ffe(j)*fk_escl(i,alpi,l)
                            enddo
                            do alpj = 1, nsinge*ndim
                                jj = jjn + (1+nfhe+1-1)*ndim + alpj
                                mmat(ii,jj) = mmat(ii,jj)-&
                                              mm*fk_escl(i,alpi,l)*fk_escl(j,alpj,k)
                            enddo
                        enddo
                    end do
                    do j = 1, nnm
                        call indent(j, ddlms, ddlmm, nnms, jjn)
                        jjn = jjn + nddle
                        if (lpenac) then
                            mm = 0.d0
                        else
                            mm = hpg*coefcr*norm(l)* jacobi*norm(k)
                        endif
                        if (lmulti) then
                            do ifh = 1, nfhe
                                iescl(1+ifh)=xcalc_heav(heavn(nfhe*(i-1)+ifh),&
                                                        heavfa(1),&
                                                        heavn(nfhe*nne+nfhm*nnm+i))
                            end do
                            do ifh = 1, nfhm
                                jmait(1+ifh)=xcalc_heav(heavn(nfhe*nne+nfhm*(j-1)+ifh),&
                                                        heavfa(2),&
                                                        heavn((1+nfhe)*nne+nfhm*nnm+j))
                            end do
                        else
                           iescl(2)=xcalc_heav(heavn(i),&
                                               hea_fa(1),&
                                               heavn(nfhe*nne+nfhm*nnm+i))
                           jmait(2)=xcalc_heav(heavn(nne+j),&
                                               hea_fa(2),&
                                               heavn((1+nfhe)*nne+nfhm*nnm+j))
                        endif
                        do iddl = 1, 1+nfhe
                            ii = iin + (iddl-1)*ndim + l
                            do  jddl = 1, 1+nfhm
                                jj = jjn + (jddl-1)*ndim + k
                                mmat(ii,jj) = mmat(ii,jj)-&
                                        iescl(iddl)*jmait(jddl)*mm*ffm(j)*ffe(i)
                                mmat(jj,ii) = mmat(jj,ii)-&
                                        iescl(iddl)*jmait(jddl)*mm*ffm(j)*ffe(i)
                            end do
                            do alpj = 1, nsingm*ndim
                                jj = jjn + (1+nfhm+1-1)*ndim + alpj
                                mmat(ii,jj) = mmat(ii,jj)-&
                                        iescl(iddl)*mm*ffe(i)*fk_mait(j,alpj,k)
                                mmat(jj,ii) = mmat(jj,ii)-&
                                        iescl(iddl)*mm*ffe(i)*fk_mait(j,alpj,k)
                            enddo
                        end do
                        do alpi = 1, nsinge*ndim
                            ii = iin + (1+nfhe+1-1)*ndim + alpi
                            do jddl = 1, 1+nfhm
                                jj = jjn + (jddl-1)*ndim + k
                                mmat(ii,jj) = mmat(ii,jj)+&
                                        jmait(jddl)*mm*ffm(j)*fk_escl(i,alpi,l)
                                mmat(jj,ii) = mmat(jj,ii)+&
                                        jmait(jddl)*mm*ffm(j)*fk_escl(i,alpi,l)
                            enddo
                            do alpj = 1, nsingm*ndim
                                jj = jjn + (1+nfhm+1-1)*ndim + alpj
                                mmat(ii,jj) = mmat(ii,jj)+&
                                        mm*fk_escl(i,alpi,l)*fk_mait(j,alpj,k)
                                mmat(jj,ii) = mmat(jj,ii)+&
                                        mm*fk_escl(i,alpi,l)*fk_mait(j,alpj,k)
                            enddo
                        enddo
                    end do
                end do
                do i = 1, nnm
                    call indent(i, ddlms, ddlmm, nnms, iin)
                    iin = iin + nddle
                    do j = 1, nnm
                        call indent(j, ddlms, ddlmm, nnms, jjn)
                        jjn = jjn + nddle
                        if (lpenac) then
                            mm = 0.d0
                        else
                            mm = hpg*coefcr*norm(l)* jacobi*norm(k)
                        endif
                        if (lmulti) then
                            do ifh = 1, nfhm
                                imait(1+ifh)=xcalc_heav(heavn(nfhe*nne+nfhm*(i-1)+ifh),&
                                            heavfa(2),&
                                            heavn((1+nfhe)*nne+nfhm*nnm+i))
                                jmait(1+ifh)=xcalc_heav(heavn(nfhe*nne+nfhm*(j-1)+ifh),&
                                            heavfa(2),&
                                            heavn((1+nfhe)*nne+nfhm*nnm+j))
                            end do
                        else
                            imait(2)=xcalc_heav(heavn(nne+i),&
                                                hea_fa(2),&
                                                heavn((1+nfhe)*nne+nfhm*nnm+i))
                            jmait(2)=xcalc_heav(heavn(nne+j),&
                                                hea_fa(2),&
                                                heavn((1+nfhe)*nne+nfhm*nnm+j))
                        endif
                        do iddl = 1, 1+nfhm
                            ii = iin + (iddl-1)*ndim + l
                            do jddl = 1, 1+nfhm
                                jj = jjn + (jddl-1)*ndim + k
                                mmat(ii,jj) = mmat(ii,jj)+&
                                    imait(iddl)*jmait(jddl)*mm*ffm(i)*ffm(j)
                            end do
                            ii = iin + (iddl-1)*ndim + l
                            do alpj = 1, nsingm*ndim
                                jj = jjn + (1+nfhm+1-1)*ndim + alpj
                                mmat(ii,jj) = mmat(ii,jj)+&
                                    imait(iddl)*mm*ffm(i)*fk_mait(j,alpj,k)
                            end do
                        end do
                        do alpi = 1, nsingm*ndim
                            ii = iin + (1+nfhm+1-1)*ndim + alpi
                            do jddl = 1, 1+nfhm
                                jj = jjn + (jddl-1)*ndim + k
                                mmat(ii,jj) = mmat(ii,jj)+&
                                    jmait(jddl)* mm *ffm(j)*fk_mait(i,alpi,l)
                            end do
                            do alpj = 1, nsingm*ndim
                                jj = jjn + (1+nfhm+1-1)*ndim + alpj
                                mmat(ii,jj) = mmat(ii,jj)+&
                                    mm*fk_mait(i,alpi,l)*fk_mait(j,alpj,k)
                            end do
                        enddo
                    end do
                end do
            end do
        end do
    else
!
! --------------------- CALCUL DE [A] ----------------------------------
!
        do k = 1, ndim
            do i = 1, nnc
                call xplma2(ndim, nne, nnes, ddles, i,&
                            nfhe, pl)
                if (lmulti) pl = pl + (heavno(i)-1)*ndim
                do j = 1, ndeple
! --- BLOCS ES:CONT, CONT:ES
                    call indent(j, ddles, ddlem, nnes, jjn)
                    mm = hpg*ffc(i)*jacobi*norm(k)
                    do alpj = 1, ndim*nsinge
                          jj = jjn + alpj
                          mmat(pl,jj) = mmat(pl,jj)+fk_escl(j,alpj,k) * mm
                          mmat(jj,pl) = mmat(jj,pl)+fk_escl(j,alpj,k) * mm
                    enddo
                end do
            end do
        end do
!
! --------------------- CALCUL DE [AU]----------------------------------
!
        do k = 1, ndim
            do l = 1, ndim
                do i = 1, ndeple
                    do j = 1, ndeple
! --- BLOCS ES:ES
                        if (lpenac) then
                            mm = 0.d0
                        else
                            mm = hpg*coefcr*norm(l)* jacobi*norm(k)
                        endif
                        call indent(i, ddles, ddlem, nnes, iin)
                        call indent(j, ddles, ddlem, nnes, jjn)
                        do alpi = 1, ndim*nsinge
                          ii = iin + alpi
                            do alpj = 1, ndim*nsinge
                                jj = jjn + alpj
                                mmat(ii,jj) = mmat(ii,jj)+fk_escl(j,alpj,k)*fk_escl(i,alpi,l)*mm
                            enddo
                        enddo
                    end do
                end do
            end do
        end do
    endif
! --------------------- CALCUL DE [C] ----------------------------------
!
!-------------- SEULEUMENT EN METHODE PENALISEE ------------------------
!
    if (lpenac) then
        do i = 1, nnc
            do j = 1, nnc
                call xplma2(ndim, nne, nnes, ddles, i,&
                            nfhe, pli)
                call xplma2(ndim, nne, nnes, ddles, j,&
                            nfhe, plj)
                mmat(pli,plj) = mmat(pli,plj)-hpg*ffc(j)*ffc(i)*jacobi/coefcp
            end do
        end do
    endif
!
end subroutine
