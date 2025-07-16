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
subroutine ajellt(ligretZ, meshZ, nbCell, listCell, &
                  phenomZ, modelisaZ)
!
    implicit none
!
#include "asterfort/crelgt.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/ligretDebug.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: ligretZ, meshZ
    integer(kind=8), intent(in) :: nbCell
    integer(kind=8), intent(in) :: listCell(nbCell)
    character(len=*), intent(in) :: phenomZ, modelisaZ
!
! --------------------------------------------------------------------------------------------------
!
! Create LIGRET of virtual elements for contact - To use with lgtlgr/crelgt suboutines
!
! In  ligret           : name of LIGRET
! In  mesh             : name of mesh
! In  nbCell           : number of cells to apply
! In  listCell         : list of cells
! In  phenom           : phenomenon
! In  modelisa         : modelisation
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: debug = ASTER_FALSE
    integer(kind=8) :: iCell, idlima, idlino, idlity
    integer(kind=8) :: idpoma, idpono, modelisaNume
    integer(kind=8) :: iret, elemTypeNume, jdpm
    integer(kind=8) :: jdtm, nbCellAffe, nbCellMaxi, lopomx
    integer(kind=8) :: matard, nbCellAvail, nbmail, nbmax
    integer(kind=8) :: nlolim, cellNume, cellTypeNume
    parameter(nbmail=10000)
    character(len=8) :: mesh
    character(len=16) :: phenom, modelisa
    character(len=19) :: ligret
    character(len=24) :: typmai
    integer(kind=8), pointer :: apma(:) => null()
    character(len=16), pointer :: phen(:) => null()
    integer(kind=8), pointer :: apno(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    integer(kind=8), pointer :: vnbma(:) => null()
    integer(kind=8), pointer :: mata(:) => null()
    character(len=16), pointer :: mode(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    mesh = meshZ
    phenom = phenomZ
    modelisa = modelisaZ
    ligret = ligretZ
    typmai = mesh//'.TYPMAIL'
    matard = 0
!
! --- ON VERIFIE SI LE LIGRET EXISTE ,S'IL N'EXISTE PAS, ON LE CREE :
!     -------------------------------------------------------------
    call jeexin(ligret//'.LGRF', iret)
    if (iret .eq. 0) then
        call crelgt('V', ligret)
    end if
!
! --- VECTEUR DE LA LISTE DES MAILLES CUMULEES DU LIGRET :
!     --------------------------------------------------
    call jeveuo(ligret//'.LIMA', 'E', idlima)
!
! --- VECTEUR DES TYPES DES MAILLES CUMULEES DU LIGRET :
!     ------------------------------------------------
    call jeveuo(ligret//'.LITY', 'E', idlity)
!
! --- NOM DE LA MODELISATION :
!     ----------------------
    call jeveuo(ligret//'.MODE', 'E', vk16=mode)
!
! --- NOM DU PHENOMENE :
!     ----------------
    call jeveuo(ligret//'.PHEN', 'E', vk16=phen)
!
! --- TABLEAU DE POINTEURS DANS LA LISTE DES MAILLES :
!     ----------------------------------------------
    call jeveuo(ligret//'.POMA', 'E', idpoma)
!
! --- TABLEAU DE POINTEURS DANS LA LISTE DES NOEUDS :
!     ---------------------------------------------
    call jeveuo(ligret//'.PONO', 'E', idpono)
!
! --- NOM DU MAILLAGE :
!     ---------------
    call jeveuo(ligret//'.LGRF', 'E', vk8=lgrf)
!
! --- NOMBRE DE MAILLES TARDIVES :
!     --------------------------
    call jeveuo(ligret//'.MATA', 'E', vi=mata)
!
! --- VECTEUR DE LA LISTE DES NOEUDS CUMULES DU LIGRET :
!     ------------------------------------------------
    call jeveuo(ligret//'.LINO', 'E', idlino)
!
! --- NOMBRE D'AFFECTATIONS DE MAILLES :
!     --------------------------------
    call jeveuo(ligret//'.APMA', 'E', vi=apma)
!
! --- NOMBRE D'AFFECTATIONS DE NOEUDS :
!     -------------------------------
    call jeveuo(ligret//'.APNO', 'E', vi=apno)
!
! --- NOMBRE D'AFFECTATIONS DE MAILLES :
!     --------------------------------
    call jeveuo(ligret//'.NBMA', 'E', vi=vnbma)
!
    vnbma(1) = vnbma(1)+nbCell
!
! --- ON AFFECTE UNE FOIS POUR TOUTES LE NOM DU MAILLAGE :
!     --------------------------------------------------
    if (iret .eq. 0) then
        lgrf(1) = mesh
    end if
!
! --- VERIFICATION DE L'ADEQUATION DE L'AFFECTATION DES MAILLES
! --- A LA LISTE DES MAILLES CUMULEES :
!     ===============================
    if (listCell(1) .gt. 0 .and. nbCell .ge. 1) then
!
! ---   NOMBRE DE MAILLES DEJA AFFECTEES :
!       --------------------------------
        call jelira(ligret//'.LIMA', 'LONUTI', nbCellAffe)
!
! ---   LONGUEUR DU VECTEUR LIGRET.LIMA :
!       -------------------------------
        call jelira(ligret//'.LIMA', 'LONMAX', nbCellMaxi)
!
! ---   NOMBRE DE MAILLES DISPONIBLES :
!       -----------------------------
        nbCellAvail = nbCellMaxi-nbCellAffe
!
! ---   REAJUSTEMENT EVENTUEL DES VECTEURS LIMA ET LITY :
!       -----------------------------------------------
        if (nbCell .gt. nbCellAvail) then
            nlolim = nbCell-nbCellAvail
            nbmax = nbCellMaxi+max(nlolim, nbmail)
            call juveca(ligret//'.LIMA', nbmax)
            call jeveuo(ligret//'.LIMA', 'E', idlima)
            call juveca(ligret//'.LITY', nbmax)
            call jeveuo(ligret//'.LITY', 'E', idlity)
        end if
!
! ---   VERIFICATION DE L'ADEQUATION DE LA TAILLE DU VECTEUR
! ---   DES POINTEURS DANS LA LISTE DE MAILLES :
!       --------------------------------------
!
! ---   NOMBRE D'AFFECTATIONS DE MAILLES :
!       --------------------------------
        apma(1) = apma(1)+1
!
! ---   LONGUEUR DU VECTEUR LIGRET.POMA :
!       -------------------------------
        call jelira(ligret//'.POMA', 'LONMAX', lopomx)
!
! ---   REAJUSTEMENT EVENTUEL DU VECTEUR POMA :
!       -------------------------------------
        if (apma(1) .ge. lopomx) then
            call juveca(ligret//'.POMA', 2*lopomx)
            call jeveuo(ligret//'.POMA', 'E', idpoma)
        end if
    end if

! - AFFECTATION DU TYPE DES MAILLES
    call jenonu(jexnom('&CATA.'//phenom(1:13)//'.MODL', modelisa), modelisaNume)
    call jeveuo(jexnum('&CATA.'//phenom, modelisaNume), 'L', jdpm)
    call jeveuo(typmai, 'L', jdtm)

! - AFFECTATION DE LA LISTE DES MAILLES CUMULEES
    do iCell = 1, nbCell
        zi(idlima+zi(idpoma+apma(1)-1)+iCell-1) = listCell(iCell)
        cellNume = listCell(iCell)
        cellTypeNume = zi(jdtm+cellNume-1)
        elemTypeNume = zi(jdpm+cellTypeNume-1)
        zi(idlity+zi(idpoma+apma(1)-1)+iCell-1) = elemTypeNume
    end do
!
    vnbma(1) = vnbma(1)+nbCell
!
! ---   VECTEUR DE POINTEURS DANS LE VECTEUR DES MAILLES :
!       ------------------------------------------------
    zi(idpoma+apma(1)) = zi(idpoma+apma(1)-1)+nbCell
!
    call jeecra(ligret//'.LIMA', 'LONUTI', zi(idpoma+apma(1)))

! - AFFECTATION DE LA MODELISATION AU LIGRET
    mode(1) = modelisa

! - AFFECTATION DU PHENOMENE AU LIGRET
    phen(1) = phenom

! - For debug
    if (debug) then
        call ligretDebug(ligret)
    end if
!
    call jedema()
!
end subroutine
