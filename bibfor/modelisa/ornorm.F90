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
subroutine ornorm(mesh, listCellNume, nbCell, reorie, norien, nconex, &
                  onlySkin1D_)
!
    use mesh_module, only: checkCellsAreSkin
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/indiis.h"
#include "asterfort/infniv.h"
#include "asterfort/iorim1.h"
#include "asterfort/iorim2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/utmavo.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nbCell
    integer(kind=8), pointer :: listCellNume(:)
    aster_logical, intent(in) :: reorie
    integer(kind=8), intent(out) :: norien, nconex
    aster_logical, intent(in), optional :: onlySkin1D_
!.======================================================================
!
!   ORNORM  --  LE BUT EST QUE TOUTES LES MAILLES DE LA LISTE SOIENT
!               ORIENTEES COMME LA PREMIERE MAILLE DE LA LISTE.
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NOMA           IN    K8      NOM DU MAILLAGE
!    LISTMA         IN    I       LISTE DES MAILLES A REORIENTER
!    NBMAIL         IN    I       NB DE MAILLES DE LA LISTE
!    NORIEN        VAR            NOMBRE DE MAILLES REORIENTEES
!.========================= DEBUT DES DECLARATIONS ====================
! -----  VARIABLES LOCALES
    integer(kind=8) :: iliste
    integer(kind=8) :: iCell, cellNume, norieg, lliste
    integer(kind=8) :: im1, im2, ico, ibid(1)
    integer(kind=8) :: p1, p2, ifm, niv, p3, p4
    integer(kind=8) :: jdesm1, jdesm2
    integer(kind=8) :: nbmavo, indi, im3
    integer(kind=8), parameter :: zero = 0
    aster_logical :: hasSkin1D, hasSkin2D, onlySkin1D
    character(len=1) :: lect
    character(len=2) :: kdim
    character(len=8) :: cellName
    character(len=24) :: nomavo
    integer(kind=8), pointer :: ori1(:) => null()
    integer(kind=8), pointer :: ori2(:) => null()
    integer(kind=8), pointer :: ori3(:) => null()
    integer(kind=8), pointer :: ori4(:) => null()
    character(len=8), pointer :: ori5(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
!
#define pasori(iCell) ori1(iCell).eq.0
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
    if (nbCell .eq. 0) goto 999
!
    call infniv(ifm, niv)

    onlySkin1D = ASTER_FALSE
    if (present(onlySkin1D_)) then
        onlySkin1D = onlySkin1D_
    end if
    nconex = 0

! - Options
    lect = 'L'
    if (reorie) then
        lect = 'E'
    end if

! - Access to mesh datastructures
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', p2)
    call jeveuo(mesh//'.CONNEX', lect, p1)

! - Working vectors
    AS_ALLOCATE(vi=ori1, size=nbCell)
    AS_ALLOCATE(vi=ori2, size=nbCell)
    AS_ALLOCATE(vi=ori3, size=nbCell)
    AS_ALLOCATE(vi=ori4, size=nbCell)
    AS_ALLOCATE(vk8=ori5, size=nbCell)

! - Check type of cells (only skin)
    call checkCellsAreSkin(mesh, &
                           nbCell, listCellNume, &
                           onlySkin1D, &
                           ori5, &
                           hasSkin1D, hasSkin2D)

! - Prepare lists
    do iCell = 1, nbCell
        cellNume = listCellNume(iCell)
        ori3(iCell) = zi(p2+cellNume)-zi(p2-1+cellNume)
        ori4(iCell) = zi(p2+cellNume-1)
    end do
!
!
! --- RECUPERATION DES MAILLES VOISINES DU GROUP_MA :
!     ---------------------------------------------
    kdim = '  '
    if (hasSkin1D) kdim = '1D'
    if (hasSkin2D) kdim = '2D'
    nomavo = '&&ORNORM.MAILLE_VOISINE '
    call utmavo(mesh, kdim, listCellNume, nbCell, 'V', &
                nomavo, zero, ibid)
    call jeveuo(jexatr(nomavo, 'LONCUM'), 'L', p4)
    call jeveuo(nomavo, 'L', p3)
!
    norieg = 0
!
! --- LA BOUCLE 100 DEFINIT LES CONNEXES
!
    do iCell = 1, nbCell
        cellNume = listCellNume(iCell)
! ----- SI LA MAILLE N'EST PAS ORIENTEE ON L'ORIENTE
        if (pasori(iCell)) then
            if (niv .eq. 2) then
                cellName = int_to_char8(cellNume)
                call utmess('I', 'MESH3_9', sk=cellName)
            end if
            nconex = nconex+1
            ori1(iCell) = 1
            lliste = 0
            iliste = 0
            ori2(lliste+1) = iCell
!
! ------- ON ORIENTE TOUTES LES MAILLES DU CONNEXE
!
200         continue
!
            im1 = ori2(iliste+1)
            jdesm1 = ori4(im1)
! ------- ON ESSAYE D'ORIENTER LES MAILLES VOISINES
            nbmavo = zi(p4+im1)-zi(p4-1+im1)
            do im3 = 1, nbmavo
                indi = zi(p3+zi(p4+im1-1)-1+im3-1)
                im2 = indiis(listCellNume, indi, 1, nbCell)
                if (im2 .eq. 0) goto 210
                cellNume = listCellNume(im2)
                if (pasori(im2)) then
                    jdesm2 = ori4(im2)
!             VERIFICATION DE LA CONNEXITE ET REORIENTATION EVENTUELLE
                    if (hasSkin1D) ico = iorim1(zi(p1+jdesm1-1), zi(p1+jdesm2-1), reorie)
                    if (hasSkin2D) ico = iorim2( &
                                         zi(p1+jdesm1-1), ori3(im1), zi(p1+jdesm2-1), ori3(im2), &
                                         reorie &
                                         )
!             SI MAILLES CONNEXES
                    if (ico .ne. 0) then
                        ori1(im2) = 1
                        lliste = lliste+1
                        ori2(lliste+1) = im2
                        if (reorie .and. niv .eq. 2) then
                            cellName = int_to_char8(cellNume)
                            if (ico .lt. 0) then
                                call utmess('I', 'MESH3_7', sk=cellName)
                            else
                                call utmess('I', 'MESH3_8', sk=cellName)
                            end if
                        end if
                    end if
!
!             SI ORIENTATIONS CONTRAIRES
                    if (ico .lt. 0) norieg = norieg+1
!
                end if
210             continue
            end do
            iliste = iliste+1
            if (iliste .le. lliste) goto 200
        end if
    end do
!
    norien = norien+norieg
!
    AS_DEALLOCATE(vi=ori1)
    AS_DEALLOCATE(vi=ori2)
    AS_DEALLOCATE(vi=ori3)
    AS_DEALLOCATE(vi=ori4)
    AS_DEALLOCATE(vk8=ori5)
    call jedetr(nomavo)
!
999 continue
    call jedema()
end subroutine
