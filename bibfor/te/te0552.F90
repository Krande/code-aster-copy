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

subroutine te0552(option, nomte)
    implicit none
!
    character(len=16) :: option, nomte
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES DEPLACEMENTS AUX POINTS DE GAUSS
!
!     OPTION TRAITEE  ==>  DEPL_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter  :: mxnoeu = 4, mxnpg = 4, nbcompo = 6, nbdepl = 3
    integer(kind=8), parameter  :: sp_couche_dkt = 3, sp_couche_gri = 1
!
    integer(kind=8)             :: icmp, indga, indno, indx, ino, ipg, iptc, nbcou, icou
    integer(kind=8)             :: ndim, nno, nnos, npg, kdec
    integer(kind=8)             :: icacoq, jdeplga, jvf
    integer(kind=8)             :: jdepg, jgeom, ipoids, idfdx, jgano, jnbspi
!
    logical             :: elem_dkt, elem_gri
!
    real(kind=8)        :: pgl(3, 3)
    real(kind=8)        :: deplno(nbcompo*mxnoeu), deplga(nbcompo*mxnpg), deplsp(nbdepl)
    real(kind=8)        :: epcoqu, hcouche, zic, zmin, excent
!
    character(len=4)    :: fami
!
! --------------------------------------------------------------------------------------------------
!
    if (option .ne. 'DEPL_ELGA') then
        ASSERT(.false.)
    end if
    !
    ! Les éléments traités
    elem_dkt = (nomte .eq. 'MEDKTR3') .or. (nomte .eq. 'MEDKQU4') .or. &
               (nomte .eq. 'MEDSQU4') .or. (nomte .eq. 'MEDSTR3') .or. &
               (nomte .eq. 'MEQ4QU4') .or. (nomte .eq. 'MET3TR3')
    elem_gri = (nomte .eq. 'MEGCQU4') .or. (nomte .eq. 'MEGCTR3')
    ASSERT(elem_dkt .or. elem_gri)
    !
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=jvf, jdfde=idfdx, jgano=jgano)
    !
    ASSERT(nno .le. mxnoeu)
    ASSERT(npg .le. mxnpg)
    ! La Géométrie
    call jevech('PGEOMER', 'L', jgeom)
    ! Les déplacements aux noeuds
    call jevech('PDEPLAR', 'L', jdepg)
    ! Matrice de passage Global-->Local
    if (nno .eq. 3) then
        call dxtpgl(zr(jgeom), pgl)
    else if (nno .eq. 4) then
        call dxqpgl(zr(jgeom), pgl)
    end if
    ! Passage des déplacements dans le repère local
    call utpvgl(nno, nbcompo, pgl, zr(jdepg), deplno)
    !
    call jevech('PDEPLGA', 'E', jdeplga)
    !
    ! Calcul des déplacements aux points de gauss de l'élément support
    deplga(:) = 0.0
    do ipg = 1, npg
        kdec = (ipg-1)*nno
        do ino = 1, nno
            do icmp = 1, nbcompo
                indga = (ipg-1)*nbcompo+icmp
                indno = (ino-1)*nbcompo+icmp
                deplga(indga) = deplga(indga)+deplno(indno)*zr(jvf+kdec+ino-1)
            end do
        end do
    end do
    !
    ! Calcul des déplacements aux sous-points
    ! Caractéristiques des coques
    call jevech('PCACOQU', 'L', icacoq)
    epcoqu = zr(icacoq)
    !
    if (elem_dkt) then
        excent = zr(icacoq+4)
        call jevech('PNBSP_I', 'L', jnbspi)
        nbcou = zi(jnbspi)
        if (nbcou .le. 0) then
            call utmess('F', 'ELEMENTS_46')
        end if
        !
        hcouche = epcoqu/nbcou
        zmin = excent-epcoqu*0.50
        !
        do icou = 1, nbcou
            do iptc = 1, sp_couche_dkt
                ! Altitude des points dans la couche
                if (iptc .eq. 1) then
                    zic = zmin+(icou-1)*hcouche
                else if (iptc .eq. 2) then
                    zic = zmin+(icou-1)*hcouche+hcouche*0.50
                else
                    zic = zmin+(icou-1)*hcouche+hcouche
                end if
                !
                do ipg = 1, npg
                    indga = (ipg-1)*nbcompo
                    deplsp(1) = deplga(indga+1)+deplga(indga+5)*zic
                    deplsp(2) = deplga(indga+2)-deplga(indga+4)*zic
                    deplsp(3) = deplga(indga+3)
                    ! Adresse du sous-point
                    indx = jdeplga+nbdepl*(sp_couche_dkt*(nbcou*(ipg-1)+icou-1)+iptc-1)
                    ! Passage des déplacements dans le repère global
                    call utpvlg(1, nbdepl, pgl, deplsp, zr(indx))
                end do
            end do
        end do
    else if (elem_gri) then
        excent = zr(icacoq+3)
        nbcou = 1
        icou = 1
        iptc = 1
        !
        zic = excent
        do ipg = 1, npg
            indga = (ipg-1)*nbcompo
            deplsp(1) = deplga(indga+1)+deplga(indga+5)*zic
            deplsp(2) = deplga(indga+2)-deplga(indga+4)*zic
            deplsp(3) = deplga(indga+3)
            ! Adresse du sous-point
            indx = jdeplga+nbdepl*(sp_couche_gri*(nbcou*(ipg-1)+icou-1)+iptc-1)
            ! Passage des déplacements dans le repère global
            call utpvlg(1, nbdepl, pgl, deplsp, zr(indx))
        end do
    end if
!
end subroutine
