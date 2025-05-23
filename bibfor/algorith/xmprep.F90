! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine xmprep(cface, contac, elref, elrefc, elc, &
                  ffc, ffp, fpg, iaint, ibasec, &
                  iptint, ifa, igeom, ipgf, jac, &
                  jlst, lact, nd, ndim, ninter, &
                  nlact, nno, nnos, nptf, nvit, &
                  rr, singu, tau1, tau2, ka, &
                  mu, jbaslo, jstno, jlsn, fk)
! aslint: disable=W1504
    implicit none
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/assert.h"
#include "asterfort/lteatt.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xjacf2.h"
#include "asterfort/xjacff.h"
#include "asterfort/xmoffc.h"
#include "asterfort/xxmmvd.h"
! PREPARATION MATRICE OU VECTEUR EL DE CONTACT DANS UN TE
! IN CFACE  : TABLEAU CONNECTIVITE FACETTES DE CONTACT
! IN CONTAC : TYPE DE CONTACT (1: P1P1, 3: P2P1)
! IN ELREF  : TYPE ELEMENT PARENT
! IN ELREFC : EQUIVALENT LINEAIRE DE ELREF (POUR LE P2P1)
! OUT FFC   : FONCTIONS DE FORME DE CONTACT
! OUT FFP   : FONCTIONS DE FORME ELEMENT PARENT
! IN  IAINT : ADRESSE TOPOFAC.AI FISSURE COURANTE
! IN  IFA   : NUMERO DE FACETTE
! IN IGEOM  : ADRESSE COORDONNEES NOEUD EL PARENT
! IN IPGF   : NUMERO PT GAUSS FACETTE
! OUT JAC   : PRODUIT JACOBIEN*POIDS
! IN JLSN   : ADRESSE LSN
! IN JLST   : ADRESSE LST
! IN LACT   : LAGRANGES ACTIFS
! OUT ND    : NORMALE A SURFACE DE CONTACT AU PG
! IN  NDIM
! IN NINTER
! IN NLACT
! IN NNO
! IN NNOS
! OUT NVIT  : ARETE VITALE OU NON
! OUT RR    : RACINE DU RAYON A LE POINTE DE FISSURE
! IN SINGU  : PRESENCE OU NON DE LA SINGULARITE
! OUT TAU1  : 1ERE TANGENTE A LE SURFACE DE CONTACT AU PG
! OUT TAU2  : 2EME TANGENTE (EN 3D)
    integer :: cface(30, 6), contac
    integer :: i, iaint, ibasec, ifa, igeom, ipgf, iptint, j, n
    integer :: jlst, k, lact(8), ndim, ninter, nlact, nno
    integer :: nnos, nptf, nvit, singu, zxain
    integer :: alp
    real(kind=8) :: dfbid(nno, 3), ffc(8), ffp(27), ffpc(27), g(3), jac
    real(kind=8) :: nd(3), lst, r, rr, tau1(3), tau2(3), x(4)
    real(kind=8) :: fk_mait(27, 3, 3), fk_escl(27, 3, 3)
    real(kind=8), optional :: ka, mu, fk(27, 3, 3)
    integer, optional :: jstno, jlsn, jbaslo
    character(len=8) :: elc, elref, elrefc, fpg
    aster_logical :: axi
    aster_logical :: largs_fk
!
!
! --- CALCUL DE JAC (PRODUIT DU JACOBIEN ET DU POIDS)
! --- ET DES FF DE L'ELEMENT PARENT AU POINT DE GAUSS
! --- ET LA NORMALE ND ORIENTÉE DE ESCL -> MAIT
!
    if (ndim .eq. 3) then
        call xjacff(elref, elrefc, elc, ndim, fpg, &
                    iptint, ifa, cface, ipgf, nno, &
                    nnos, igeom, ibasec, g, jac, &
                    ffp, ffpc, dfbid, nd, tau1, &
                    tau2)
    else if (ndim .eq. 2) then
        call xjacf2(elref, elrefc, elc, ndim, fpg, &
                    iptint, ifa, cface, nptf, ipgf, &
                    nno, nnos, igeom, ibasec, g, &
                    jac, ffp, ffpc, dfbid, nd, &
                    tau1)
    end if
!Check for axi-symmetric elements
    axi = lteatt('AXIS', 'OUI')
    if (axi) then
        r = 0.d0
        do n = 1, nno
            r = r+ffp(n)*zr(igeom-1+2*(n-1)+1)
        end do
        if (r .eq. 0d0) r = 1.d-7
        ASSERT(r .gt. 0d0)
        jac = jac*r
!
    end if
!
! --- CALCUL DES FONCTIONS DE FORMES DE CONTACT
!
    if ((contac .eq. 1) .or. (contac .eq. 2)) then
        call xmoffc(lact, nlact, nno, ffp, ffc)
    else if (contac .eq. 3) then
        call xmoffc(lact, nlact, nnos, ffpc, ffc)
    end if
!
! --- CE POINT DE GAUSS EST-IL SUR UNE ARETE?
!
    k = 0
    do i = 1, ninter
        if (k .eq. 0) then
            x(4) = 0.d0
            do j = 1, ndim
                x(j) = zr(iptint-1+ndim*(i-1)+j)
            end do
            do j = 1, ndim
                x(4) = x(4)+(x(j)-g(j))*(x(j)-g(j))
            end do
            x(4) = sqrt(x(4))
            if (x(4) .lt. 1.d-12) then
                k = i
                goto 17
            end if
        end if
17      continue
    end do
!
    zxain = xxmmvd('ZXAIN')
    if (k .ne. 0) then
        nvit = nint(zr(iaint-1+zxain*(k-1)+5))
    else
        nvit = 0
    end if
!      IL NE FAUT PAS UTILISER NVIT SI LE SCHEMA D'INTEGRATION
!      NE CONTIENT PAS DE NOEUDS
    if ((fpg(1:3) .eq. 'FPG') .or. (fpg .eq. 'GAUSS') .or. (fpg .eq. 'XCON')) nvit = 1
!
! --- CALCUL DE RR = SQRT(DISTANCE AU FOND DE FISSURE)
!
    if (singu .eq. 1) then
        fk(:, :, :) = 0.d0
        lst = 0.d0
        do i = 1, nno
            lst = lst+zr(jlst-1+i)*ffp(i)
        end do
        r = abs(lst)
        rr = sqrt(r)
! --- CALCUL FORMEL DU SAUT MOYEN POUR LES FCTS SINGULIERES
        largs_fk = present(ka) .and. &
                   present(mu) .and. &
                   present(jbaslo) .and. &
                   present(fk) .and. &
                   present(jlsn) .and. &
                   present(jstno)
        ASSERT(largs_fk)
        call xcalfev_wrap(ndim, nno, zr(jbaslo), zi(jstno), +1.d0, &
                          zr(jlsn), zr(jlst), zr(igeom), ka, mu, &
                          ffp, fk_mait, face='MAIT')
        call xcalfev_wrap(ndim, nno, zr(jbaslo), zi(jstno), -1.d0, &
                          zr(jlsn), zr(jlst), zr(igeom), ka, mu, &
                          ffp, fk_escl, face='ESCL')
        do i = 1, nno
            do alp = 1, ndim
                do j = 1, ndim
                    fk(i, alp, j) = (fk_mait(i, alp, j)-fk_escl(i, alp, j))/2.d0
                end do
            end do
        end do
!
    end if
!
end subroutine
