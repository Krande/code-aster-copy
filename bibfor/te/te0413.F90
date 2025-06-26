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

subroutine te0413(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/glrc_recup_mate.h"
#include "asterfort/gquad4.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/Behaviour_type.h"
    character(len=16) :: option, nomte
!
!
! FONCTIONS REALISEES:
!
!      CALCUL DE LA DENSITE DE DISSIPATION
!      A L'EQUILIBRE POUR LES ELEMENTS DKTG ET LA LOI GLRC_DM
!      .SOIT AUX POINTS D'INTEGRATION : OPTION 'DISS_ELGA'
!      .SOIT L INTEGRALE PAR ELEMENT  : OPTION 'DISS_ELEM'
!
!      OPTIONS : 'DISS_ELGA'
!                'DISS_ELEM'
!
! ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    integer(kind=8) :: npgmx
    parameter(npgmx=4)
!
    real(kind=8) :: pgl(3, 3)
    real(kind=8) :: qsi, eta, xyzl(3, 4), jacob(5), poids, cara(25)
    real(kind=8) :: disse(npgmx), dse
    real(kind=8) :: ep, seuil
!
    integer(kind=8) :: ndim, nno, nnoel, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: jgeom, ipg, idener, imate
    integer(kind=8) :: icacoq, jvari, nbvar
!
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: valk(2)
    aster_logical :: dkq, lkit
!
    if (nomte .eq. 'MEDKQG4') then
        dkq = .true.
    else if (nomte .eq. 'MEDKTG3') then
        dkq = .false.
    else
        call utmess('F', 'ELEMENTS_34', sk=nomte)
    end if
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnoel, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                     jgano=jgano)
!
    call jevech('PGEOMER', 'L', jgeom)
    call jevech('PCOMPOR', 'L', vk16=compor)
    if (nno .eq. 3) then
        call dxtpgl(zr(jgeom), pgl)
    else if (nno .eq. 4) then
        call dxqpgl(zr(jgeom), pgl)
    end if
!
    lkit = compor(RELA_NAME) (1:7) .eq. 'KIT_DDI'
!
    if ((compor(RELA_NAME) (1:7) .eq. 'GLRC_DM') .or. &
        (lkit .and. (compor(CREEP_NAME) (1:7) .eq. 'GLRC_DM'))) then
!
        call jevech('PCACOQU', 'L', icacoq)
!
        call utpvgl(nno, 3, pgl, zr(jgeom), xyzl)
!
        if (dkq) then
            call gquad4(xyzl, cara)
        else
            call gtria3(xyzl, cara)
        end if
!
        read (compor(NVAR), '(I16)') nbvar
        ep = zr(icacoq)
!
        if (option .eq. 'DISS_ELGA') then
            call jevech('PVARIGR', 'L', jvari)
        else if (option .eq. 'DISS_ELEM') then
            call jevech('PVARIPR', 'L', jvari)
        end if
!
        call r8inir(npgmx, 0.d0, disse, 1)
        dse = 0.0d0
!
! ---- BOUCLE SUR LES POINTS D'INTEGRATION :
!      ===================================
        do ipg = 1, npg
!
            qsi = zr(icoopg-1+ndim*(ipg-1)+1)
            eta = zr(icoopg-1+ndim*(ipg-1)+2)
            if (dkq) then
                call jquad4(xyzl, qsi, eta, jacob)
                poids = zr(ipoids+ipg-1)*jacob(1)
            else
                poids = zr(ipoids+ipg-1)*cara(7)
            end if
!
            call jevech('PMATERC', 'L', imate)
!
            call glrc_recup_mate(zi(imate), compor(RELA_NAME), .false._1, ep, seuil=seuil)
!
!  --    CALCUL DE LA DENSITE D'ENERGIE POTENTIELLE ELASTIQUE :
!        ==========================================================
            if ((option .eq. 'DISS_ELGA') .or. (option .eq. 'DISS_ELEM')) then
!
                disse(ipg) = (zr(jvari-1+(ipg-1)*nbvar+1)+zr(jvari-1+(ipg-1)*nbvar+2))*seuil
                dse = dse+disse(ipg)*poids
!
            end if
        end do
!
! ---- RECUPERATION DU CHAMP DES DENSITES D'ENERGIE DE DEFORMATION
! ---- ELASTIQUE EN SORTIE
!      -------------------
        if (option .eq. 'DISS_ELGA') then
            call jevech('PDISSPG', 'E', idener)
        else if (option .eq. 'DISS_ELEM') then
            call jevech('PDISSD1', 'E', idener)
        end if
!
! --- OPTIONS DISS_ELGA
!     ==============================
        if (option .eq. 'DISS_ELGA') then
            do ipg = 1, npg
                zr(idener-1+(ipg-1)*1+1) = disse(ipg)
            end do
!
! --- OPTION DISS_ELEM
!     ================
        else if (option .eq. 'DISS_ELEM') then
            zr(idener-1+1) = dse
        end if
    elseif (compor(RELA_NAME) (1:4) .eq. 'DHRC') then
        call jevech('PCACOQU', 'L', icacoq)
        call utpvgl(nno, 3, pgl, zr(jgeom), xyzl)
        if (dkq) then
            call gquad4(xyzl, cara)
        else
            call gtria3(xyzl, cara)
        end if
!
        read (compor(NVAR), '(I16)') nbvar
!
        if (option .eq. 'DISS_ELGA') then
            call jevech('PVARIGR', 'L', jvari)
        else if (option .eq. 'DISS_ELEM') then
            call jevech('PVARIPR', 'L', jvari)
        end if
!
        call r8inir(npgmx, 0.d0, disse, 1)
        dse = 0.0d0
!
! ---- BOUCLE SUR LES POINTS D'INTEGRATION :
!      ===================================
        do ipg = 1, npg
!
            qsi = zr(icoopg-1+ndim*(ipg-1)+1)
            eta = zr(icoopg-1+ndim*(ipg-1)+2)
            if (dkq) then
                call jquad4(xyzl, qsi, eta, jacob)
                poids = zr(ipoids+ipg-1)*jacob(1)
            else
                poids = zr(ipoids+ipg-1)*cara(7)
            end if
!
!  --    CALCUL DE LA DENSITE D'ENERGIE POTENTIELLE ELASTIQUE :
!        ==========================================================
            if ((option .eq. 'DISS_ELGA') .or. (option .eq. 'DISS_ELEM')) then
!
                disse(ipg) = zr(jvari-1+(ipg-1)*nbvar+9)
                dse = dse+disse(ipg)*poids
!
            end if
        end do
!
! ---- RECUPERATION DU CHAMP DES DENSITES D'ENERGIE DE DEFORMATION
! ---- ELASTIQUE EN SORTIE
!      -------------------
        if (option .eq. 'DISS_ELGA') then
            call jevech('PDISSPG', 'E', idener)
        else if (option .eq. 'DISS_ELEM') then
            call jevech('PDISSD1', 'E', idener)
        end if
!
! --- OPTIONS DISS_ELGA
!     ==============================
        if (option .eq. 'DISS_ELGA') then
            do ipg = 1, npg
                zr(idener-1+(ipg-1)*1+1) = disse(ipg)
            end do
!
! --- OPTION DISS_ELEM
!     ================
        else if (option .eq. 'DISS_ELEM') then
            zr(idener-1+1) = dse
        end if
    else
!      RELATION NON PROGRAMMEE
        valk(1) = option
        valk(2) = compor(RELA_NAME) (1:7)
        call utmess('A', 'ELEMENTS4_63', nk=2, valk=valk)
    end if
!
end subroutine
