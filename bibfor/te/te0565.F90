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

subroutine te0565(nomopt, nomte)
!
    implicit none
!
!     BUTS: .CALCUL DE L'ENERGIES DE DEFORMATION ELASTIQUE - CAS X-FEM
!
! -----------------------------------------------------------------
!
!  OPTION ENEL_ELEM : CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!  ================   DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!  EN HPP
!   ENELAS =  SOMME_VOLUME((SIG_T*(1/D)*SIG).DV)
!
!        OU  .SIG       EST LE TENSEUR DES CONTRAINTES DE CAUCHY
!            .D         EST LE TENSEUR DE HOOKE
!
!  EN GRANDES DEFORMATIONS SIMO MIEHE POUR ELAS OU VMIS_ISOT
!   ENERLAS = ENERGIE ELASTIQUE SPECIFIQUE
!           = K(0.5(J^2-1)-lnJ)+0.5mu(tr(J^(-2/3)be)-3)
!           SI PRESENCE DE THERMIQUE, ON AJOUTE UNE CORRECTION
!           SPECIFIQUE PRESENTEE DANS LA DOC R
!  EN GRANDES DEFORMATIONS GDEF_LOG
!   ENERELAS = SOMME_VOLUME((T_T*(1/D)*T).DV)
!        OU  .T       EST LE TENSEUR DES CONTRAINTES DU FORMALISME
!            .D         EST LE TENSEUR DE HOOKE
!
!   REMARQUE : EN GRANDE DEFORMATION ON INTEGRE SUR LE VOLUME INITIALE
! -----------------------------------------------------------------
!          ELEMENTS ISOPARAMETRIQUES 2D ET 3D
!
!          OPTIONS : 'ENEL_ELEM'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    character(len=16) :: nomte, nomopt
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/enelpg.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/reeref.h"
#include "asterfort/tecach.h"
#include "asterfort/teattr.h"
#include "asterfort/ltequa.h"
!-----------------------------------------------------------------------
    integer(kind=8) :: idene1
    integer(kind=8) :: idfde, idsig, idvari, igeom, imate, itemps
    integer(kind=8) :: ipoids, ivf
    integer(kind=8) :: nbsgm, nbsig, nbvari, ndim, nno
    integer(kind=8) :: npg, iret, i, jtab(7)
    parameter(nbsgm=6)
    real(kind=8) :: enelas
    real(kind=8) :: deux, trois
    real(kind=8) :: un, undemi, untier, welas, wtotal
    real(kind=8) :: zero
    real(kind=8) :: sigma(nbsgm)
    real(kind=8) :: angl_naut(3), instan
    real(kind=8) :: f(3, 3), r
    character(len=16) :: compor(3)
    integer(kind=8) :: jpintt, jpmilt, jcnset, jlonch
    integer(kind=8) :: nse, nnop
    character(len=8) :: elrefp, elrese(6), fami(6), enr
    real(kind=8) :: coorse(81), xg(3), xe(3), ff(27)
    real(kind=8) :: jac
    integer(kind=8) :: irese, idecpg, idebs, idebv
    integer(kind=8) :: j, in, ino, ise, kpg, n
    aster_logical :: grand, axi
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
!-----------------------------------------------------------------------
!
!
! ---- INITIALISATIONS :
!
    zero = 0.0d0
    undemi = 0.5d0
    un = 1.0d0
    deux = 2.0d0
    trois = 3.0d0
    untier = 1.0d0/3.0d0
    enelas = zero
    welas = zero
    wtotal = zero
    instan = zero
!
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
!
!
! --- ELEMENT DE REFERENCE PARENT : RECUP DE NDIM ET NNOP
    call elref1(elrefp)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop)
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NNO, NPG ET IVF
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! --- TYPE DE MODELISATION
!
    axi = lteatt('AXIS', 'OUI')
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT :
!
    nbsig = nbsigm()
!
! ---- RECUPERATION DES COORDONNEES DES CONNECTIVITES :
!
    call jevech('PGEOMER', 'L', igeom)
!
! ---- RECUPERATION DU MATERIAU :
!
    call jevech('PMATERC', 'L', imate)
!
! ---- RECUPERATION DES CHAMPS X-FEM :
!
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PLONCHA', 'L', jlonch)
!     PROPRES AUX ELEMENTS 1D ET 2D (QUADRATIQUES)
    call teattr('S', 'XFEM', enr, iret)
    if ((iret .eq. 0) .and. ltequa(elrefp, enr)) &
        call jevech('PPMILTO', 'L', jpmilt)
!
! ---- RECUPERATION  DES DONNEEES RELATIVES AU REPERE D'ORTHOTROPIE :
    call getElemOrientation(ndim, nnop, igeom, angl_naut)
!
! ---- RECUPERATION DU CHAMP DE DEPLACEMENTS AUX NOEUDS  :
!
!    RQ: - on suppose qu'on est en h.p.p. et qu'il n'est donc pas
!          necessaire de disposer du champ de deplacement pour
!          calculer l'energie ; on suppose que la variation du volume
!          de l'element est negligeable
!        - cette hypothese permet d'eviter d'avoir a creer un champ
!          de deplacement avec suffisament de composantes (H1X, etc)
!    call jevech('PDEPLR', 'L', idepl)
!
! ---- RECUPERATION DU CHAMP DE CONTRAINTES AUX POINTS D'INTEGRATION :
!
    call jevech('PCONTPR', 'L', idsig)
!
! ---- RECUPERATION DE L'INSTANT DE CALCUL
!      -----------------------------------
    call tecach('NNO', 'PINSTR', 'L', iret, iad=itemps)
    if (itemps .ne. 0) instan = zr(itemps)
!
! ----RECUPERATION DU TYPE DE COMPORTEMENT  :
!     N'EXISTE PAS EN LINEAIRE
    call tecach('NNO', 'PCOMPOR', 'L', iret, nval=7, &
                itab=jtab)
    compor(1) = 'ELAS'
    compor(2) = ' '
    compor(3) = 'PETIT'
    if (iret .eq. 0) then
        compor(1) = zk16(jtab(1))
        compor(3) = zk16(jtab(1)+2)
    end if
!
!     GRANDES DEFORMATIONS
!
    if ((compor(3) .eq. 'SIMO_MIEHE') .or. (compor(3) .eq. 'GDEF_LOG')) then
        grand = .true.
    else
        grand = .false.
    end if
!
! ----   RECUPERATION DU CHAMP DE VARIABLES INTERNES  :
!        N'EXISTE PAS EN LINEAIRE
    call tecach('ONO', 'PVARIPR', 'L', iret, nval=7, &
                itab=jtab)
    if (iret .eq. 0) then
        idvari = jtab(1)
        nbvari = max(jtab(6), 1)*jtab(7)
    else
        idvari = 1
        nbvari = 0
    end if
!
! --- RECUPERATION DU CHAMP POUR STOCKER L'ENERGIE ELSTIQUE
    call jevech('PENERD1', 'E', idene1)
!
! --- RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = zi(jlonch-1+1)
!
! --- BOUCLE D'INTEGRATION SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!       BOUCLE SUR LES SOMMETS DU SOUS-TRIA (DU SOUS-SEG)
        do in = 1, nno
            ino = zi(jcnset-1+nno*(ise-1)+in)
            do j = 1, ndim
                if (ino .lt. 1000) then
                    coorse(ndim*(in-1)+j) = zr(igeom-1+ndim*(ino-1)+j)
                else if (ino .gt. 1000 .and. ino .lt. 2000) then
                    coorse(ndim*(in-1)+j) = zr(jpintt-1+ndim*(ino-1000-1)+j)
                else if (ino .gt. 2000 .and. ino .lt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-2000-1)+j)
                else if (ino .gt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-3000-1)+j)
                end if
            end do
        end do
!
!-----------------------------------------------------------------------
!         BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELT
!-----------------------------------------------------------------------
!
        do kpg = 1, npg
!
!     --- COORDONNÉES DU PT DE GAUSS DANS LE REPÈRE RÉEL : XG
            xg(:) = 0.d0
            do i = 1, ndim
                do n = 1, nno
                    xg(i) = xg(i)+zr(ivf-1+nno*(kpg-1)+n)*coorse(ndim*(n-1)+i)
                end do
            end do
!
!           JUSTE POUR CALCULER LES FF
!
            call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                        xe, ff)
!
!           CALCULER LE JACOBIEN DE LA TRANSFO SSTET->SSTET REF
!           AVEC LES COORDONNEES DU SOUS-ELEMENT
            if (ndim .eq. 2) then
                call dfdm2d(nno, kpg, ipoids, idfde, coorse, &
                            jac)
            else if (ndim .eq. 3) then
                call dfdm3d(nno, kpg, ipoids, idfde, coorse, &
                            jac)
            end if
!
! -         CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE):
            if (axi) then
                r = 0.d0
                do n = 1, nnop
                    r = r+ff(n)*zr(igeom-1+2*(n-1)+1)
                end do
!
                ASSERT(r .gt. 0d0)
!               MODIFICATION DU JACOBIEN
                jac = jac*r
            end if

!           DEBUT DE LA ZONE MEMOIRE DE SIG ET VI CORRESPONDANTE
            idecpg = npg*(ise-1)
            idebs = nbsig*idecpg
            idebv = nbvari*idecpg

!     --- TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION COURANT :
!
            do i = 1, nbsig
                sigma(i) = zr(idsig+idebs+(kpg-1)*nbsig+i-1)
            end do
!
!           CALCUL DU GRADIENT DE LA TRANSFORMATION
!
!           Hypothese : on est en H.P.P.
            ASSERT(.not. grand)
!
            f(1:3, 1:3) = 0.d0
            do i = 1, 3
                f(i, i) = 1.d0
            end do
!
            call enelpg('XFEM', zi(imate), instan, kpg, angl_naut, &
                        compor, f, sigma, nbvari, &
                        zr(idvari+idebv+(kpg-1)*nbvari), enelas)
!
            welas = welas+enelas*jac
!
        end do
!
    end do
!
    zr(idene1) = welas
!
end subroutine
