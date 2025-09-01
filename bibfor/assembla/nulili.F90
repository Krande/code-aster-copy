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
subroutine nulili(nbLigr, listLigr, lili, base, gran_name, &
                  igds, mesh, nec, nlili, modeLocZ_)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/nbgrel.h"
#include "asterfort/typele.h"
!
    integer(kind=8), intent(in) :: nbLigr
    character(len=24), pointer :: listLigr(:)
    character(len=24), intent(in):: lili
    character(len=1), intent(in):: base
    character(len=8), intent(out) :: gran_name
    integer(kind=8), intent(out) :: igds
    character(len=8), intent(out) :: mesh
    integer(kind=8), intent(out) :: nec
    integer(kind=8), intent(out) :: nlili
    character(len=*), optional, intent(in) :: modeLocZ_
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering - Create LILI objects
!
! --------------------------------------------------------------------------------------------------
!
! In  listLigr       : pointer to list of LIGREL
! In  nbLigr         : number of LIGREL in list
! In  lili           : name of LILI object to create
! In  base           : JEVEUX base to create LILI object
! Out gran_name      : name of GRANDEUR to number
! Out igds           : index of GRANDEUR to number
! Out mesh           : name of mesh
! Out nec            : number of coding integers for GRANDEUR
! Out nlili          : size of LILI object
! In  modeLoc        : local mode
!
!----------------------------------------------------------------------
! ATTENTION : NE PAS FAIRE JEMARQ/JEDEMA CAR CETTE ROUTINE
!             RECOPIE DES ADRESSES JEVEUX DANS .ADNE ET .ADLI
!----------------------------------------------------------------------
!
!
!    --- DESCRIPTION DES OBJETS ADNE ET ADLI ---
!     ADNE (1          ) = NBRE DE MAILLES DU MAILLAGE
!     ADNE (2          ) = 0
!     ADNE (3          ) = 0
!     ADLI (1          ) = 0
!     ADLI (2          ) = 0
!     ADLI (3          ) = 0
!     POUR 2<=ILI<=NLILI
!     ADNE (3*(ILI-1)+1) = NBRE MAX D'OBJETS DE LA COLLECTION
!                            LILI(ILI).NEMA
!     ADNE (3*(ILI-1)+2) = ADRESSE DE L'OBJET LILI(ILI).NEMA
!     ADNE (3*(ILI-1)+3) = ADRESSE DU VECTEUR DES LONG. CUMULEES DE
!                            LILI(ILI).NEMA
!     ADLI (3*(ILI-1)+1) = NBRE MAX D'OBJETS DE LA COLLECTION
!                            LILI(ILI).LIEL
!     ADLI (3*(ILI-1)+2) = ADRESSE DE L'OBJET LILI(ILI).LIEL
!     ADLI (3*(ILI-1)+3) = ADRESSE DU VECTEUR DES LONG. CUMULEES DE
!                            LILI(ILI).LIEL
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: ligrelMesh, answer
    character(len=16) :: phenom, ligrelPhenom, nomte
    character(len=19) :: prefix
    character(len=24) :: ligrName
    integer(kind=8) :: iad, iGrel, nbCell
    integer(kind=8) :: iligr, iret
    integer(kind=8) :: nbgr, nbsup, jmoloc, imode, cellTypeNume
    character(len=8) :: modeLoc
    integer(kind=8), pointer :: adli(:) => null()
    integer(kind=8), pointer :: adne(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
!    call jemarq() FORBIDDEN !
!
    prefix = lili(1:14)
    nlili = nbLigr+1
    ASSERT(nlili .gt. 1)

! - Local mode
    modeLoc = ' '
    if (present(modeLocZ_)) then
        modeLoc = modeLocZ_
    end if

! - Create main LILI repertory object
    call jecreo(lili, base//' N  K24')
    call jeecra(lili, 'NOMMAX', nlili)
    call jecroc(jexnom(lili, '&MAILLA'))

! - Creation of temporary objects ADNE et ADLI
    call jecreo(prefix//'.ADNE', 'V V I')
    call jeecra(prefix//'.ADNE', 'LONMAX', 3*nlili)
    call jeveuo(prefix//'.ADNE', 'E', vi=adne)
    call jecreo(prefix//'.ADLI', 'V V I')
    call jeecra(prefix//'.ADLI', 'LONMAX', 3*nlili)
    call jeveuo(prefix//'.ADLI', 'E', vi=adli)

! - Save temporary objects ADNE et ADLI
    do iligr = 1, nbLigr
        ligrName = listLigr(iligr)

! ----- Only one phenomenon
        call dismoi('PHENOMENE', ligrName, 'LIGREL', repk=ligrelPhenom)
        if (iligr .eq. 1) then
            phenom = ligrelPhenom
        else
            ASSERT(phenom .eq. ligrelPhenom)
        end if

! ----- Same mesh
        call dismoi('NOM_MAILLA', ligrName, 'LIGREL', repk=ligrelMesh)
        if (iligr .eq. 1) then
            mesh = ligrelMesh
        end if
        ASSERT(mesh .eq. ligrelMesh)

! ----- Create object in collection
        call jecroc(jexnom(lili, ligrName))

! ----- Set ADNE/ADLI objects
        call dismoi('EXI_ELEM', ligrName, 'LIGREL', repk=answer)
        if (answer(1:3) .ne. 'NON') then
            call jeexin(ligrName(1:19)//'.NEMA', iret)
            if (iret .ne. 0) then
!
!---- ADNE(3*(iligr)+1)=NBRE DE MAILLES SUP DU LIGREL NOMLI
!
                call jelira(ligrName(1:19)//'.NEMA', 'NUTIOC', nbsup)
                adne(1+3*(iligr)) = nbsup
                call jeveut(ligrName(1:19)//'.NEMA', 'L', iad)
                adne(1+3*(iligr)+1) = iad
                call jeveut(jexatr(ligrName(1:19)//'.NEMA', 'LONCUM'), 'L', iad)
                adne(1+3*(iligr)+2) = iad
            else
                adne(1+3*(iligr)) = 0
                adne(1+3*(iligr)+1) = 2**30
                adne(1+3*(iligr)+2) = 2**30
            end if
!
!---- ADLI(3*(iligr)+1)=NBRE DE MAILLES DU LIGREL NOMLI
!
            call jelira(ligrName(1:19)//'.LIEL', 'NUTIOC', nbgr)
            adli(1+3*(iligr)) = nbgr
            call jeveut(ligrName(1:19)//'.LIEL', 'L', iad)
            adli(1+3*(iligr)+1) = iad
            call jeveut(jexatr(ligrName(1:19)//'.LIEL', 'LONCUM'), 'L', iad)
            adli(1+3*(iligr)+2) = iad
        end if
    end do
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCell)
    adne(1) = nbCell

! - Information about GRANDEUR
    if (modeLoc .eq. ' ') then
        call dismoi('NOM_GD', phenom, 'PHENOMENE', repk=gran_name)
    else
        ligrName = listLigr(1) (1:19)
        do iGrel = 1, nbgrel(ligrName)
            cellTypeNume = typele(ligrName, iGrel)
            call jenuno(jexnum('&CATA.TE.NOMTE', cellTypeNume), nomte)
            call jenonu(jexnom('&CATA.TE.NOMMOLOC', nomte//modeLoc), imode)
            if (imode .gt. 0) then
                call jeveuo(jexnum('&CATA.TE.MODELOC', imode), 'L', jmoloc)
                call jenuno(jexnum('&CATA.GD.NOMGD', zi(jmoloc-1+2)), gran_name)
                goto 30
            end if
        end do
        ASSERT(ASTER_FALSE)
30      continue

! --- Pour certains RESU_ELEM, il faut changer le nom
! Il en manque sûrement (voir creprn.F90)
        if (gran_name(1:4) == "MDEP" .or. gran_name(1:4) == "VDEP" &
            .or. gran_name(1:4) == "MDNS") then
            gran_name = "DEPL_"//gran_name(6:6)
        elseif (gran_name(1:4) == "MTEM" .or. gran_name(1:4) == "VTEM" &
                .or. gran_name(1:4) == "MTNS") then
            gran_name = "TEMP_"//gran_name(6:6)
        elseif (gran_name(1:4) == "MPRE" .or. gran_name(1:4) == "VPRE") then
            gran_name = "PRES_"//gran_name(6:6)
        elseif (gran_name(1:4) == "MSIZ" .or. gran_name(1:4) == "VSIZ") then
            gran_name = "SIZZ_"//gran_name(6:6)
        elseif (gran_name(1:4) == "MZNS" .or. gran_name(1:4) == "VNEU") then
            gran_name = "NEUT_"//gran_name(6:6)
        end if
    end if
    call jenonu(jexnom('&CATA.GD.NOMGD', gran_name), igds)
    ASSERT(igds .ne. 0)
!
    nec = nbec(igds)
!
!    call jedema() FORBIDDEN !
!
end subroutine
