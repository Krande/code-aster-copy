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

subroutine te0576(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/enelpg.h"
#include "asterfort/eps1mc.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/nmgeom.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!.......................................................................
!
! FONCTIONS REALISEES:
!
!      CALCUL DE LA DENSITE D'ENERGIE POTENTIELLE THERMOELASTIQUE
!      A L'EQUILIBRE POUR LES ELEMENTS ISOPARAMETRIQUES 3D
!      .SOIT AUX POINTS D'INTEGRATION : OPTION 'ENEL_ELGA'
!
!      OPTIONS : 'ENEL_ELGA'
!
!      CALCUL DE LA DENSITE D'ENERGIE TOTALE
!      A L'EQUILIBRE POUR LES ELEMENTS ISOPARAMETRIQUES 3D
!      .SOIT AUX POINTS D'INTEGRATION : OPTION 'ETOT_ELGA'
!
!      OPTIONS : 'ETOT_ELGA'
!
!
! ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    integer(kind=8) :: jgano, ndim, nno, i, nnos, ipoids, ivf, nbnomx, nbcont, npg1
    integer(kind=8) :: nbsig, igau, isig, igeom, itemps, nbvari, imate, idener
    integer(kind=8) :: idfde, idepl, ideplm, idepmm, idvari, idsig, idsigm, mxcmel, iret
    integer(kind=8) :: idenem, jtab(7)
    parameter(nbnomx=27)
    parameter(nbcont=6)
    parameter(mxcmel=162)
    real(kind=8) :: epsi(nbcont), angl_naut(3), instan, zero, undemi, enelem
    real(kind=8) :: enerpg(nbnomx)
    real(kind=8) :: nharm, deux, integ1, integ2, integ, r
    real(kind=8) :: epsim(nbcont), delta(nbcont), epss(mxcmel)
    real(kind=8) :: epssm(mxcmel), sigmm(nbcont), sigma(nbcont), f(3, 3)
    real(kind=8) :: poids, epsbid(6), dfdbid(27*3)
    character(len=4) :: fami
    character(len=16) :: compor(3)
    aster_logical :: grand
! DEB ------------------------------------------------------------------
!
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
!      ------------------------
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!      -----------------------------------------
    nbsig = nbsigm()
!
! --- INITIALISATIONS :
!     -----------------
    zero = 0.0d0
    undemi = 0.5d0
    deux = 2.0d0
    nharm = zero
    enelem = zero
    instan = 0.d0
!
    do i = 1, nbnomx
        enerpg(i) = zero
    end do
!
! ---- RECUPERATION DES COORDONNEES DES CONNECTIVITES
!      ----------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
    if (option(1:4) .eq. 'ENEL') then
!
! ----   RECUPERATION DU MATERIAU
!        ------------------------
        call jevech('PMATERC', 'L', imate)
!
! ----   RECUPERATION  DES DONNEES RELATIVES AU REPERE D'ORTHOTROPIE
!        -----------------------------------------------------------
        call getElemOrientation(ndim, nno, igeom, angl_naut)
!
! ---    RECUPERATION DU CHAMP DE DEPLACEMENT A L'INSTANT COURANT :
!        --------------------------------------------------------
        call jevech('PDEPLAR', 'L', idepl)
!
! ----   RECUPERATION DU CHAMP DE CONTRAINTES AUX POINTS D'INTEGRATION
!        -------------------------------------------------------------
        call jevech('PCONTRR', 'L', idsig)
!
! ----    RECUPERATION DE L'INSTANT DE CALCUL
!         -----------------------------------
        call tecach('ONO', 'PINSTR', 'L', iret, iad=itemps)
        if (itemps .ne. 0) instan = zr(itemps)
!
! ----   RECUPERATION DU CHAMP DE VARIABLES INTERNES  :
!        N'EXISTE PAS EN LINEAIRE
        call tecach('ONO', 'PVARIGR', 'L', iret, nval=7, &
                    itab=jtab)
        if (iret .eq. 0) then
            idvari = jtab(1)
            nbvari = max(jtab(6), 1)*jtab(7)
        else
            idvari = 1
            nbvari = 0
        end if
!
    end if
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
!
! --- CAS DU CALCUL DE LA DENSITE D'ENERGIE TOTALE :
!     ============================================
    if (option(1:4) .eq. 'ETOT') then
!
        if (grand) then
            call utmess('F', 'COMPOR1_79', sk=compor(3))
        end if
!
! ---    RECUPERATION DU CHAMP DE DEPLACEMENT A L'INSTANT COURANT :
!        --------------------------------------------------------
        call jevech('PDEPLR', 'L', idepl)
!
! ---    RECUPERATION EVENTUELLE DU CHAMP DE DEPLACEMENT A
! ---    L'INSTANT PRECEDENT :
!        -------------------
        call tecach('NNO', 'PDEPLM', 'L', iret, iad=ideplm)
        if (ideplm .ne. 0) then
            call jevech('PDEPLM', 'L', idepmm)
        end if
!
! ---    RECUPERATION DU CHAMP DE CONTRAINTES AUX POINTS D'INTEGRATION
! ---    A L'INSTANT COURANT :
!        -------------------
        call jevech('PCONTPR', 'L', idsig)
!
! ---    RECUPERATION EVENTUELLE DU CHAMP DE CONTRAINTES A
! ---    L'INSTANT PRECEDENT :
!        -------------------
        call tecach('NNO', 'PCONTMR', 'L', iret, iad=idsigm)
        if (idsigm .ne. 0) then
            call jevech('PCONTMR', 'L', idsigm)
        end if
!
! ---    CALCUL DU CHAMP DE DEFORMATIONS AU PREMIER ORDRE
! ---    CORRESPONDANT AU CHAMP DE DEPLACEMENT COURANT :
!        ---------------------------------------------
        call eps1mc(nno, ndim, nbsig, npg1, ipoids, &
                    ivf, idfde, zr(igeom), zr(idepl), nharm, &
                    epss)
!
! ---    CALCUL EVENTUEL DU CHAMP DE DEFORMATIONS AU PREMIER ORDRE
! ---    CORRESPONDANT AU CHAMP DE DEPLACEMENT A L'INSTANT PRECEDENT :
!        -----------------------------------------------------------
        if (ideplm .ne. 0) then
            call eps1mc(nno, ndim, nbsig, npg1, ipoids, &
                        ivf, idfde, zr(igeom), zr(idepmm), nharm, &
                        epssm)
        end if
!
    end if
!
! ---- BOUCLE SUR LES POINTS D'INTEGRATION :
!      ===================================
    do igau = 1, npg1
!
!  --    CALCUL DU JACOBIEN AU POINT D'INTEGRATION COURANT :
!        -------------------------------------------------
        call dfdm3d(nno, igau, ipoids, idfde, zr(igeom), &
                    poids)
!
        do isig = 1, nbsig
            epsi(isig) = zero
        end do
!
!  --    CALCUL DE LA DENSITE D'ENERGIE POTENTIELLE THERMOELASTIQUE :
!        ==========================================================
        if (option(1:4) .eq. 'ENEL') then
!
! --- TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION COURANT :
!
            do i = 1, nbsig
                sigma(i) = zr(idsig+(igau-1)*nbsig+i-1)
            end do
!
! --- CALCUL DU JACOBIEN AU POINT D'INTEGRATION COURANT :
            call nmgeom(3, nno, .false._1, grand, zr(igeom), &
                        igau, ipoids, ivf, idfde, zr(idepl), &
                        .true._1, poids, dfdbid, f, epsbid, &
                        r)
!
! ---     CALCUL DE L'ENERGIE ELASTIQUE AU POINT D'INTEGRATION COURANT
!
            call enelpg(fami, zi(imate), instan, igau, angl_naut, &
                        compor, f, sigma, nbvari, &
                        zr(idvari+(igau-1)*nbvari), enerpg(igau))
!
!
!  --    CALCUL DE LA DENSITE D'ENERGIE TOTALE :
!        =====================================
        else if (option(1:4) .eq. 'ETOT') then
!
!  --      TENSEURS DES DEFORMATIONS  ET DES CONTRAINTES AU PAS DE
!  --      TEMPS COURANT ET AU PAS DE TEMPS PRECEDENT S'IL Y A LIEU,
!  --      AU POINT D'INTEGRATION COURANT :
!          ------------------------------
            do i = 1, nbsig
                epsi(i) = epss(i+(igau-1)*nbsig)
                if (ideplm .ne. 0) then
                    epsim(i) = epssm(i+(igau-1)*nbsig)
                end if
                sigma(i) = zr(idsig+(igau-1)*nbsig+i-1)
                if (idsigm .ne. 0) then
                    sigmm(i) = zr(idsigm+(igau-1)*nbsig+i-1)
                end if
            end do
!
            if (ideplm .ne. 0) then
                do i = 1, nbsig
                    delta(i) = epsi(i)-epsim(i)
                end do
            else
                do i = 1, nbsig
                    delta(i) = 0.d0
                end do
            end if
!
!  --      CALCUL DES TERMES A SOMMER POUR OBTENIR LA DENSITE
!  --      D'ENERGIE TOTALE :
!          ----------------
            integ1 = sigma(1)*delta(1)+sigma(2)*delta(2)+sigma(3)*delta(3)+deux*sigma(4)*d&
                     &elta(4)+deux*sigma(5)*delta(5)+deux*sigma(6)*delta(6)
!
            if (ideplm .ne. 0 .and. idsigm .ne. 0) then
                integ2 = sigmm(1)*delta(1)+sigmm(2)*delta(2)+sigmm(3)*delta(3)+deux*sigmm(4&
                         &)*delta(4)+deux*sigmm(5)*delta(5)+deux*sigmm(6)*delta(6)
!
                enerpg(igau) = undemi*(integ1+integ2)
            else
!
!  --        CAS D'ORDRE NUMERO 1 :
!            --------------------
                integ = sigma(1)*epsi(1)+sigma(2)*epsi(2)+sigma(3)*epsi(3)+deux*sigma(4)*e&
                        &psi(4)+deux*sigma(5)*epsi(5)+deux*sigma(6)*epsi(6)
!
                enerpg(igau) = undemi*integ
!
            end if
!
            enelem = enelem+enerpg(igau)*poids
!
        end if
!
    end do
!
! ---- RECUPERATION DU CHAMP DES DENSITES D'ENERGIE DE DEFORMATION
! ---- ELASTIQUE EN SORTIE
!      -------------------
    call jevech('PENERDR', 'E', idener)
!
! ---- OPTIONS ENEL_* ET ETOT_*
!      ==============================
    if (option(1:4) .eq. 'ETOT') then
        call jevech('PENERDM', 'L', idenem)
!
! ----   OPTION ETOT_ELGA
!        ================
        if (option .eq. 'ETOT_ELGA') then
            do igau = 1, npg1
                zr(idener+igau-1) = zr(idenem+igau-1)+enerpg(igau)
            end do
!
! ----   OPTION ETOT_ELEM
!        ================
        else if (option .eq. 'ETOT_ELEM') then
            zr(idener) = zr(idenem)+enelem
        end if
    else
!
! ----   OPTION ENEL_ELGA
!        ================
        if (option .eq. 'ENEL_ELGA') then
            do igau = 1, npg1
                zr(idener+igau-1) = enerpg(igau)
            end do
        end if
    end if
end subroutine
