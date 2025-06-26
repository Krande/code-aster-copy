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
!
subroutine caliel(valeTypeZ, loadZ, modelZ)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterfort/aflrch.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/caarle.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/nueffe.h"
#include "asterfort/raco3d.h"
#include "asterfort/rapo2d.h"
#include "asterfort/rapo3d.h"
#include "asterfort/rapoco.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: loadZ, valeTypeZ, modelZ
!
! --------------------------------------------------------------------------------------------------
!
!     MODELISATION DU RACCORD ENTRE DES ELEMENTS
!     AYANT DES MODELISATIONS DIFFERENTES PAR DES RELATIONS
!     LINEAIRES ENTRE DDLS.
!     CES RELATIONS SONT AFFECTEES A LA CHARGE CHARGZ.
!     TYPES DES RACCORDS TRAITES :
!       1) RACCORD POUTRE-3D PAR DES RELATIONS LINEAIRES
!          ENTRE LES NOEUDS DES MAILLES DE SURFACE MODELISANT
!          LA TRACE DE LA SECTION DE LA POUTRE SUR LE MASSIF 3D
!          ET LE NOEUD DE LA POUTRE DONNE PAR L'UTILISATEUR
!
!       2) RACCORD POUTRE-COQUE PAR DES RELATIONS LINEAIRES
!          ENTRE LES NOEUDS DES MAILLES DE BORD DE COQUE MODELISANT
!          LA TRACE DE LA SECTION DE LA POUTRE SUR A COQUE
!          ET LE NOEUD DE LA POUTRE DONNE PAR L'UTILISATEUR
! -------------------------------------------------------
!  FONREZ        - IN    - K4   - : 'REEL' OU 'FONC'
!  CHARGZ        - IN    - K8   - : NOM DE LA SD CHARGE
!                - JXVAR -      -
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: renumSans = "SANS"
    character(len=16), parameter :: factorKeyword = "LIAISON_ELEM"
    character(len=8) :: model, load
    character(len=14), parameter :: numeDof = '&&CALIEL.NUMED'
    character(len=16) :: option
    character(len=19) :: modelLigrel
    character(len=19), parameter :: lisrel = '&&CALIEL.RLLISTE'
    integer(kind=8) :: iocc, nbOcc, iop, nbLigr
    character(len=24), pointer :: listLigr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    load = loadZ
    model = modelZ
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
!
    call getfac(factorKeyword, nbOcc)
    if (nbOcc .ne. 0) then
! ----- Create numbering on model
        nbLigr = 1
        AS_ALLOCATE(vk24=listLigr, size=nbLigr)
        listLigr(1) = modelLigrel
        call nueffe(nbLigr, listLigr, 'VV', numeDof, renumSans, model)
        AS_DEALLOCATE(vk24=listLigr)

! ----- Create relation
        do iocc = 1, nbOcc
            call getvtx(factorKeyword, 'OPTION', iocc=iocc, scal=option, nbret=iop)
            if (option .eq. '3D_POU') then
                call rapo3d(numeDof, iocc, valeTypeZ, lisrel, loadZ)
            else if (option .eq. 'COQ_3D') then
                call raco3d(numeDof, iocc, valeTypeZ, lisrel, loadZ)
            else if (option .eq. '3D_POU_ARLEQUIN') then
                call caarle(numeDof, iocc, lisrel, loadZ)
            else if (option .eq. '2D_POU') then
                call rapo2d(numeDof, iocc, valeTypeZ, lisrel, loadZ)
            else if (option .eq. '3D_TUYAU') then
                call rapo3d(numeDof, iocc, valeTypeZ, lisrel, loadZ)
            else if (option .eq. 'PLAQ_POUT_ORTH') then
                call rapo3d(numeDof, iocc, valeTypeZ, lisrel, loadZ)
            else if (option .eq. 'COQ_POU') then
                call rapoco(numeDof, iocc, valeTypeZ, lisrel, loadZ)
            else if (option .eq. 'COQ_TUYAU') then
                call rapoco(numeDof, iocc, valeTypeZ, lisrel, loadZ)
            else
                ASSERT(ASTER_FALSE)
            end if
        end do

! ----- Affect relations to load
        call aflrch(lisrel, load, 'NLIN')

! ----- Clean
        call jedetc('V', lisrel, 1)
        call jedetr('&&CALIEL.LIGRMO')
        call jedetr('&&CALIEL.NUMED')
!
    end if

    call jedema()
end subroutine
