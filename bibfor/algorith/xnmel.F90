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
! aslint: disable=W1306,W1504
! person_in_charge: samuel.geniaut at edf.fr
!
subroutine xnmel(nnop, nfh, nfe, ddlc, &
                 ddlm, igeom, typmod, option, imate, &
                 compor, lgpg, carcri, jpintt, cnset, &
                 heavt, lonch, basloc, instam, instap, idepl, lsn, &
                 lst, sig, vi, matuu, ivectu, &
                 codret, jpmilt, nfiss, jheavn, jstno, &
                 l_line, l_nonlin, lMatr, lVect, lSigm)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/nbsigm.h"
#include "asterfort/tecach.h"
#include "asterfort/xxnmel.h"
!
    integer(kind=8) :: nnop, imate, lgpg, codret, igeom, nfiss, jheavn
    integer(kind=8) :: cnset(4*32), heavt(*), lonch(10), ndim
    integer(kind=8) :: nfh, nfe, ddlc, ddlm
    integer(kind=8) :: ivectu, idepl, jpintt, jpmilt
    integer(kind=8) :: jstno
    character(len=8) :: typmod(*)
    character(len=16) :: option, compor(*)
    real(kind=8) :: instam, instap
    real(kind=8) :: carcri(*), vi(*), crit2(1), vi2(1), sig2(1)
    real(kind=8) :: lsn(nnop)
    real(kind=8) :: lst(nnop), matuu(*), sig(*), basloc(*)
    aster_logical, intent(in) :: l_line, l_nonlin, lMatr, lVect, lSigm
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:  PRÉLIMINAIRES AU CALCUL DES OPTIONS RIGI_MECA_TANG,
!           RAPH_MECA ET FULL_MECA  EN HYPER-ELASTICITE AVEC X-FEM
!
! --------------------------------------------------------------------------------------------------
!
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  IVF     : VALEUR  DES FONCTIONS DE FORME
! IN  NFH     : NOMBRE DE DDL HEAVISIDE (PAR NOEUD)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  DDLM    : NOMBRE DE DDL PAR NOEUD MILIEU
! IN  IGEOM   : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  OPTION  : OPTION DE CALCUL
! IN  NOMTE   : NOM DU TE
! IN  IMATE   : MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT
! IN  LGPG  : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!              CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  PINTT   : COORDONNÉES DES POINTS D'INTERSECTION
! IN  CNSET   : CONNECTIVITE DES SOUS-ELEMENTS
! IN  HEAVT   : VALEURS DE L'HEAVISIDE SUR LES SS-ELTS
! IN  LONCH   : LONGUEURS DES CHAMPS UTILISÉES
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE
! IN  IDEPL   : DEPLACEMENT A PARTIR DE LA CONF DE REF
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  PMILT   : COORDONNEES DES POINTS MILIEUX
! IN  NFISS   : NOMBRE DE FISSURES "VUES" PAR L'ÉLÉMENT
! IN  JHEAVN  : POINTEUR VERS LA DEFINITION HEAVISIDE
!
! OUT SIG     : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA)
! OUT VI      : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA)
! OUT MATUU   : MATRICE DE RIGIDITE PROFIL (RIGI_MECA_TANG ET FULL_MECA)
! OUT IVECTU  : VECTEUR FORCES NODALES (RAPH_MECA ET FULL_MECA)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elrefp, fami_se
    real(kind=8) :: coorse(81), he(nfiss)
    integer(kind=8) :: nse, npg
    integer(kind=8) :: nnops, ibid, ibid2
    integer(kind=8) :: j, ise, in, ino, idebs, idebv
    integer(kind=8) :: nbsig, idecpg, jtab(7), ncomp, iret
    integer(kind=8) :: ncompn, heavn(nnop, 5)
    integer(kind=8) :: irese, nno, ig, ifiss
    character(len=8), parameter :: elrese(6) = (/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/)
    character(len=8), parameter :: fami(6) = (/'BID ', 'XINT', 'XINT', 'BID ', 'XINT', 'XINT'/)
!
! --------------------------------------------------------------------------------------------------
!
!     ATTENTION, DEPL ET VECTU SONT ICI DIMENSIONNÉS DE TELLE SORTE
!     QU'ILS NE PRENNENT PAS EN COMPTE LES DDL SUR LES NOEUDS MILIEU
!
    call elref1(elrefp)
!
!     NOMBRE DE COMPOSANTES DE PHEAVTO (DANS LE CATALOGUE)
    call tecach('OOO', 'PHEAVTO', 'L', iret, nval=2, &
                itab=jtab)
    ncomp = jtab(2)
!
!     ELEMENT DE REFERENCE PARENT : RECUP DE NDIM
    call elrefe_info(fami='RIGI', ndim=ndim, nnos=nnops)
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NPG
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
!
!     ADRESSE DES COORD DU SOUS ELT EN QUESTION
    fami_se = fami(ndim+irese)
    if (nfe .gt. 0) then
        if (ndim .eq. 3 .and. &
            (count(zi((jstno-1+1):(jstno-1+nnop)) .eq. 2)+ &
             count(zi((jstno-1+1):(jstno-1+nnop)) .eq. 0)) .eq. nnop) then
            fami_se = 'XGEO'
        end if
    end if
!
! - Get element parameters
!
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami_se, nno=nno, npg=npg)
!
!     NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
    nbsig = nbsigm()
!    RECUPERATION DE LA DEFINITION DES DDL HEAVISIDES
    if (nfh .gt. 0) then
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, itab=jtab)
        ncompn = jtab(2)/jtab(3)
        ASSERT(ncompn .eq. 5)
        do ino = 1, nnop
            do ig = 1, ncompn
                heavn(ino, ig) = zi(jheavn-1+ncompn*(ino-1)+ig)
            end do
        end do
    end if
!     RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = lonch(1)
!
!       BOUCLE D'INTEGRATION SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!       BOUCLE SUR LES 4/3 SOMMETS DU SOUS-TETRA/TRIA
        do in = 1, nno
            ino = cnset(nno*(ise-1)+in)
            do j = 1, ndim
                if (ino .lt. 1000) then
                    coorse(ndim*(in-1)+j) = zr(igeom-1+ndim*(ino-1)+j)
                else if (ino .gt. 1000 .and. ino .lt. 2000) then
                    coorse(ndim*(in-1)+j) = zr(jpintt-1+ndim*(ino-1000- &
                                                              1)+j)
                else if (ino .gt. 2000 .and. ino .lt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-2000- &
                                                              1)+j)
                else if (ino .gt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-3000- &
                                                              1)+j)
                end if
            end do
        end do
!
!       FONCTION HEAVISIDE CSTE POUR CHAQUE FISSURE SUR LE SS-ELT
        do ifiss = 1, nfiss
            he(ifiss) = heavt(ncomp*(ifiss-1)+ise)
        end do
!
!       DEBUT DE LA ZONE MEMOIRE DE SIG ET VI CORRESPONDANTE
        idecpg = npg*(ise-1)
        idebs = nbsig*idecpg
        idebv = lgpg*idecpg
!
        if (ndim .eq. 3) then
            ASSERT(nbsig .eq. 6)
        else if (ndim .eq. 2) then
            ASSERT(nbsig .eq. 4)
        end if
!
        if (l_line) then
            call xxnmel(elrefp, elrese(ndim+irese), ndim, coorse, &
                        igeom, he, nfh, ddlc, ddlm, &
                        nnops, nfe, basloc, nnop, npg, &
                        typmod, option, imate, compor, lgpg, &
                        crit2, instam, instap, ibid, lsn, lst, idecpg, &
                        sig2, vi2, matuu, ibid2, codret, &
                        nfiss, heavn, jstno, &
                        l_line, l_nonlin, lMatr, lVect, lSigm)
        elseif (l_nonlin) then
            call xxnmel(elrefp, elrese(ndim+irese), ndim, coorse, &
                        igeom, he, nfh, ddlc, ddlm, &
                        nnops, nfe, basloc, nnop, npg, &
                        typmod, option, imate, compor, lgpg, &
                        carcri, instam, instap, idepl, lsn, lst, idecpg, &
                        sig(idebs+1), vi(idebv+1), matuu, ivectu, codret, &
                        nfiss, heavn, jstno, &
                        l_line, l_nonlin, lMatr, lVect, lSigm)
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
end subroutine
