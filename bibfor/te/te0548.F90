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
subroutine te0548(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/elelin.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tecael.h"
#include "asterfort/xjacf2.h"
#include "asterfort/xjacff.h"
#include "asterfort/xlacti.h"
#include "asterfort/xminte.h"
#include "asterfort/xmoffc.h"
#include "asterfort/xplmat.h"
#include "asterfort/xteini.h"
!
    character(len=16) :: option, nomte
!
! person_in_charge: samuel.geniaut at edf.fr
!.......................................................................
!
!                CONTACT X-FEM METHODE CONTINUE :
!             MISE A JOUR DU SEUIL DE FROTTEMENT
!             MISE A JOUR DE LA COHESION DANS LE CAS COHESIF
!
!
!  OPTION : 'XREACL' (X-FEM MISE À JOUR DU SEUIL DE FROTTEMENT)
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!......................................................................
!
!
    integer :: i, j, ifa, ipgf, isspg, ni, pli
    integer :: ideppl, jaint, jcface, jlonch, jseuil, ipoidf, ivff, idfdef
    integer :: iadzi, iazk24, ipoids, ivf, idfde, jgano, jdonco
    integer :: ndim, nno, nnos, npg, nfh, ddlc, nnom, integ, ninter, nfe
    integer :: nface, cface(30, 6), ibid, nnof, npgf, jptint
    integer :: singu, jbasec, nptf
    integer :: ddls, nddl, nnol, lact(8), nlact, igeom, ddlm
    integer :: contac
    character(len=8) :: elref, typma, fpg, elc, elrefc
    real(kind=8) :: seuil, ffi, g(3), rbid, ffp(27), ffc(8), nd(3)
    real(kind=8) :: ffpc(27), dfbid(27, 3), r3bid(3)
!......................................................................
!
!
    call elref1(elref)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!     INITIALISATION DES DIMENSIONS DES DDLS X-FEM
    call xteini(nomte, nfh, nfe, singu, ddlc, &
                nnom, ddls, nddl, ddlm, ibid, &
                contac)
!
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3) (1:8)
!
! --- ROUTINE SPECIFIQUE P2P1
!
    call elelin(contac, elref, elrefc, ibid, ibid)
!
!     DEP ACTUEL (DEPPLU) : 'PDEPL_P'
    call jevech('PDEPL_P', 'L', ideppl)
    call jevech('PDONCO', 'L', jdonco)
    call jevech('PAINTER', 'L', jaint)
    call jevech('PCFACE', 'L', jcface)
    call jevech('PLONGCO', 'L', jlonch)
    call jevech('PPINTER', 'L', jptint)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PSEUIL', 'E', jseuil)
    call jevech('PBASECO', 'L', jbasec)
!
!     RECUPERATIONS DES DONNEES SUR LE CONTACT ET
!     SUR LA TOPOLOGIE DES FACETTES
    ninter = zi(jlonch-1+1)
    nface = zi(jlonch-1+2)
    nptf = zi(jlonch-1+3)
!
    do i = 1, 30
        do j = 1, 6
            cface(i, j) = 0
        end do
    end do
!
    do i = 1, nface
        do j = 1, nptf
            cface(i, j) = zi(jcface-1+nptf*(i-1)+j)
        end do
    end do
!
!     SCHEMA D'INTEGRATION NUMERIQUE ET ELEMENT DE REFERENCE DE CONTACT
    integ = nint(zr(jdonco-1+4))
    call xminte(ndim, integ, fpg)
!
    if (ndim .eq. 3) then
        if (contac .le. 2) then
            elc = 'TR3'
        else
            elc = 'TR3'
        end if
    else if (ndim .eq. 2) then
        if (contac .le. 2) then
            elc = 'SE2'
        else
            elc = 'SE3'
        end if
    end if
!
    call elrefe_info(elrefe=elc, fami=fpg, nno=nnof, npg=npgf, jpoids=ipoidf, &
                     jvf=ivff, jdfde=idfdef)
!
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3) (1:8)
!
!     LISTE DES LAMBDAS ACTIFS
!
    call xlacti(typma, ninter, jaint, lact, nlact)
    if (contac .eq. 1) nnol = nno
    if (contac .eq. 3) nnol = nnos
!
!     BOUCLE SUR LES FACETTES
    do ifa = 1, nface
!
!       BOUCLE SUR LES POINTS DE GAUSS DES FACETTES
        do ipgf = 1, npgf
!
!         INDICE DE CE POINT DE GAUSS DANS PSEUIL
            isspg = npgf*(ifa-1)+ipgf
!
!         CALCUL DE JAC (PRODUIT DU JACOBIEN ET DU POIDS)
!         ET DES FF DE L'ELEMENT PARENT AU POINT DE GAUSS
!         ET LA NORMALE ND ORIENTÉE DE ESCL -> MAIT
            if (ndim .eq. 3) then
                call xjacff(elref, elrefc, elc, ndim, fpg, &
                            jptint, ifa, cface, ipgf, nno, &
                            nnos, igeom, jbasec, g, rbid, &
                            ffp, ffpc, dfbid, nd, r3bid, &
                            r3bid)
            else if (ndim .eq. 2) then
                call xjacf2(elref, elrefc, elc, ndim, fpg, &
                            jptint, ifa, cface, nptf, ipgf, &
                            nno, nnos, igeom, jbasec, g, &
                            rbid, ffp, ffpc, dfbid, nd, &
                            r3bid)
            end if
!
!        CALCUL DES FONCTIONS DE FORMES DE CONTACT DANS LE CAS LINEAIRE
            if (contac .eq. 1) then
                call xmoffc(lact, nlact, nno, ffp, ffc)
            else if (contac .eq. 3) then
                call xmoffc(lact, nlact, nnos, ffpc, ffc)
            else if (contac .eq. 4) then
                call xmoffc(lact, nlact, nno, ffp, ffc)
            end if
!
!         CALCUL DU NOUVEAU SEUIL A PARTIR DES LAMBDA DE DEPPLU
            seuil = 0.d0
            do i = 1, nnol
                ffi = ffc(i)
                ni = i
                call xplmat(ddls, ddlc, ddlm, nno, nnom, &
                            ni, pli)
                seuil = seuil+ffi*zr(ideppl-1+pli)
            end do
!
!         LORS D'UNE CONVERGENCE FORCEE, IL SE PEUT QUE LES REACTIONS
!         SOIENT TROP PETITES. LE POINT DOIT ETRE CONSIDERE GLISSANT.
            if (abs(seuil) .lt. 1.d-11) then
                seuil = 0.d0
            end if
!
            zr(jseuil-1+isspg) = seuil
!
        end do
    end do
!
end subroutine
