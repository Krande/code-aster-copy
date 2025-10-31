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
!
subroutine xside3(elrefp, ndim, coorse, elrese, igeom, &
                  he, nfh, ddlc, ddlm, nfe, &
                  basloc, nnop, npg, idecpg, jvMaterCode, &
                  idepl, lsn, lst, nfiss, &
                  heavn, jstno, sig)
!
    use BehaviourStrain_module
    use BehaviourStrain_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dmatmc.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epstmc.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/reeref.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xcinem.h"
#include "asterfort/xkamat.h"
#include "asterfort/xnbddl.h"
#include "jeveux.h"
!
    integer(kind=8) :: ndim, igeom, jvMaterCode, nnop, npg
    integer(kind=8) :: nfh, ddlc, nfe, idecpg
    integer(kind=8) :: nfiss, heavn(nnop, 5), jstno
    character(len=8) :: elrefp, elrese
    real(kind=8) :: basloc(9*nnop), he(nfiss), coorse(*)
    real(kind=8) :: lsn(nnop), lst(nnop), sig(6, npg)
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:  CALCUL DE L'OPTION SIEF_ELGA AVEC X-FEM EN 3D
!
! --------------------------------------------------------------------------------------------------
!
! IN  ELREFP  : ELEMENT DE REFERENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNEES DES SOMMETS DU SOUS-ELEMENT
! IN  IGEOM   : COORDONNEES DES NOEUDS DE L'ELEMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ELT
! IN  NFH     : NOMBRE DE FONCTIONS HEAVISIDE
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  DDLS    : NOMBRE DE DDL PAR NOEUD MILIEU
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE AUX NOEUDS
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS DU SOUS-ELEMENT
! IN  IMATE   : MATERIAU CODE
! IN  DEPL    : DEPLACEMENT A PARTIR DE LA CONF DE REF
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  NFISS   : NOMBRE DE FISSURES "VUES" PAR L'ELEMENT
! IN  JHEAVN  : POINTEUR VERS LA DEFINITION HEAVISIDE
!
! OUT SIG     : CONTRAINTES (SIEF_ELGA)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    character(len=4), parameter :: fami = 'XFEM'
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    aster_logical, parameter :: axi = ASTER_FALSE
    integer(kind=8) :: iNode, iSigm, iEpsi, iDim
    integer(kind=8) :: kpgXFEM, kpg, hea_se
    integer(kind=8) :: ddls, ddld, ddlm
    integer(kind=8) :: ivf, idepl
    integer(kind=8) :: nno, npgbis, nbsig, nnops, ndimb
    integer(kind=8) :: singu
    real(kind=8) :: f(3, 3), epsi(6), epsiTher(6)
    real(kind=8) :: time
    real(kind=8) :: xg(ndim), xe(ndim), ff(nnop)
    real(kind=8) :: anglNaut(3), dfdi(nnop, ndim)
    real(kind=8) :: fk(27, 3, 3), dkdgl(27, 3, 3, 3)
    real(kind=8) :: grad(3, 3)
    real(kind=8) :: s, sigmTher, d(6, 6)
    real(kind=8) :: ka, mu
    integer(kind=8) :: elasID
    character(len=16) :: elasKeyword
    type(All_Varc_Strain) :: allVarcStrain
!
!     ATTENTION, DEPL ET VECTU SONT ICI DIMENSIONNES DE TELLE SORTE
!     QU'ILS NE PRENNENT PAS EN COMPTE LES DDL SUR LES NOEUDS MILIEU
! --------------------------------------------------------------------------------------------------
!
    call get_elas_id(jvMaterCode, elasID, elasKeyword)
    ASSERT(elasID .eq. ELAS_ISOT)
    anglNaut = 0.d0
    nbsig = nbsigm()
    ASSERT(ndim .eq. 3)

! - RECUPERATION DU CHAMP DE DEPLACEMENT SUR L'ELEMENT
    call jevech('PDEPLAR', 'L', idepl)

! - Get time
    time = r8vide()
    allVarcStrain%hasTime = ASTER_FALSE

!   NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD
    call xnbddl(ndim, nfh, nfe, ddlc, ddld, ddls, singu)
    call elrefe_info(fami='RIGI', nnos=nnops)

! - INITIALISATION
    call elrefe_info(elrefe=elrese, fami='XINT', ndim=ndimb, nno=nno, npg=npgbis, &
                     jvf=ivf)
    ASSERT(npg .eq. npgbis .and. ndim .eq. ndimb)

! CALCUL DE L IDENTIFIANT DU SS ELEMENT
    hea_se = xcalc_code(nfiss, he_real=[he])

! - Loop on XFEM Gauss points
    do kpgXFEM = 1, npg
        kpg = idecpg+kpgXFEM

!       COORDONNÉES DU PT DE GAUSS DANS LE REPÈRE RÉEL : XG
        xg(:) = 0.d0
        do iDim = 1, ndim
            do iNode = 1, nno
                xg(iDim) = xg(iDim)+zr(ivf-1+nno*(kpgXFEM-1)+iNode)*coorse(ndim*(iNode-1)+iDim)
            end do
        end do

!       CALCUL DES FF
        call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                    xe, ff, dfdi=dfdi)

!       FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS ET LEURS DÉRIVÉES
        if (singu .gt. 0) then
            call xkamat(jvMaterCode, ndim, axi, ka, mu)
            call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he(1), &
                              lsn, lst, zr(igeom), ka, mu, ff, fk, dfdi, dkdgl)
        end if

!       CALCUL DES DEFORMATIONS EPS
        call xcinem(axi, igeom, nnop, nnops, idepl, &
                    ndim, he, &
                    nfiss, nfh, singu, ddls, ddlm, &
                    fk, dkdgl, ff, dfdi, f, &
                    epsi, grad, heavn)

! ----- Compute thermal strains
        epsiTher = 0.d0
        call epstmc(fami, '+', kpg, ksp, ndim, &
                    time, anglNaut, jvMaterCode, &
                    VARC_STRAIN_TEMP, allVarcStrain, &
                    epsiTher)

! ----- Hooke matrix
        call dmatmc(fami, jvMaterCode, time, '+', kpg, &
                    ksp, anglNaut, nbsig, d)

! ----- Compute stress
        do iSigm = 1, nbsig
            s = 0.d0
            sigmTher = 0.d0
            do iEpsi = 1, nbsig
                s = s+epsi(iEpsi)*d(iSigm, iEpsi)
                sigmTher = sigmTher+epsiTher(iEpsi)*d(iSigm, iEpsi)
            end do
            sig(iSigm, kpgXFEM) = s-sigmTher
        end do
        sig(4, kpgXFEM) = sig(4, kpgXFEM)*rac2
        sig(5, kpgXFEM) = sig(5, kpgXFEM)*rac2
        sig(6, kpgXFEM) = sig(6, kpgXFEM)*rac2
    end do
!
end subroutine
