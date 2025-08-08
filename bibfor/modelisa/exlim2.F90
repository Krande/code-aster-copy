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

subroutine exlim2(lismai, nbmail, ligrmoz, basez, ligrez)
    implicit none
#include "jeveux.h"
#include "asterfort/adalig.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cormgi.h"
#include "asterfort/dismoi.h"
#include "asterfort/initel.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: lismai(*), nbmail
    character(len=*) :: ligrmoz, basez, ligrez

! But : Creer le ligrel "reduit" correspondant au ligrel donné
!       mais en ne conservant que certaines mailles.
!
! in  : lismai : liste des numeros de mailles constituant le
!                ligrel a creer.
!                Remarque : les mailles de numero 0 sont ignorees.
! in  : nbmail : longueur de la liste des mailles
! in  : ligrmo : nom du ligrel referencant les mailles de lismai
!                des grels
! in  : basez  : base sur laquelle on cree le ligrel
! out : ligrez : ligrel a creer
!----------------------------------------------------------------------
!
!
    character(len=1) :: base
    character(len=8) :: noma, nomail
    character(len=19) :: ligrel, ligrmo
    character(len=24) :: cptlie
    integer(kind=8) :: i, j, lont, numvec, numail, igrel, nbmam, k
    integer(kind=8) :: lcliel, jdnb, iadm, jdli, nbCell
    integer(kind=8) :: jtyp, jnel, typele, typel1, nbtyel, itype, nmel, nbmail_nz
    integer(kind=8), pointer :: lismai_nz(:) => null()
    integer(kind=8), pointer :: liel(:) => null()
    integer(kind=8), pointer :: repe(:) => null()
    integer(kind=8), pointer :: cell(:) => null()
    integer(kind=8), pointer :: mocell(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()

    base = basez
    ligrmo = ligrmoz
    ligrel = ligrez
    call dismoi('NOM_MAILLA', ligrmo, 'LIGREL', repk=noma)
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbCell)

    call jeveuo(ligrmo//'.REPE', 'L', vi=repe)
    call jeveuo(jexatr(ligrmo//'.LIEL', 'LONCUM'), 'L', lcliel)
    call jeveuo(ligrmo//'.LIEL', 'L', vi=liel)

!   -- objet .NBNO
!   --------------
    call wkvect(ligrel//'.NBNO', base//' V I', 1, jdnb)
    zi(jdnb) = 0

!   -- objet .LGRF
!   --------------
    call jedupo(ligrmo//'.LGRF', base, ligrel//'.LGRF', .false._1)

!   -- objet .TYFE
!   --------------
    call wkvect(ligrel//'.TYFE', base//' V I', nbCell, vi=cell)
    call jeveuo(ligrmo//'.TYFE', 'L', vi=mocell)
    cell = 0

!   -- on retire les mailles 0 de lismai ainsi que les mailles qui ne
!      font pas partie du ligrel :
!   ------------------------------------------------------------------
    AS_ALLOCATE(vi=lismai_nz, size=nbmail)
    nbmail_nz = 0
    do k = 1, nbmail
        numail = lismai(k)
        if (numail .gt. 0) then
            igrel = repe(1+2*(numail-1))
            if (igrel .eq. 0) then
                nomail = int_to_char8(numail)
                call utmess('A', 'MODELISA4_50', sk=nomail)
                cycle
            end if
            nbmail_nz = nbmail_nz+1
            lismai_nz(nbmail_nz) = numail
            cell(numail) = mocell(numail)
        end if
    end do

    if (nbmail_nz .eq. 0) then
        call utmess('F', 'MODELISA4_51', sk=nomail)
    end if

!   -- On compte les types d'element et le nombre de mailles par type.
!      En realite, on ne "groupe" dans un grel que les mailles successives
!      portant des elements de meme type.
!      Mais le ligrel sera optimise plus tard (appel a adalig).
!   -------------------------------------------------------------------
    call wkvect('&&EXLIM1.TYPE_NOMBRE', 'V V I', 2*nbmail_nz, jtyp)
    jnel = jtyp+nbmail_nz

    typel1 = 0
    nbtyel = 0
    itype = 0
    do i = 1, nbmail_nz
        numail = lismai_nz(i)
        igrel = repe(1+2*(numail-1))
        ASSERT(igrel .gt. 0)
        iadm = zi(lcliel+igrel)
        typele = liel(iadm-1)
        if (typele .eq. typel1) then
            zi(jnel-1+itype) = zi(jnel-1+itype)+1
        else
            nbtyel = nbtyel+1
            itype = nbtyel
            typel1 = typele
            zi(jnel-1+itype) = 1
            zi(jtyp-1+nbtyel) = typele
        end if
    end do

    nbmam = 0
    do i = 1, nbtyel
        nbmam = max(nbmam, zi(jnel-1+i))
    end do

!   -- objet .LIEL
!   ==============
    cptlie = ligrel//'.LIEL'
    lont = nbtyel*(nbmam+1)
    call jecrec(cptlie, base//' V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbtyel)
    call jeecra(cptlie, 'LONT', lont)
    call jeveuo(cptlie, 'E', jdli)

!   -- stockage des groupes elements dans liel
!   ------------------------------------------
    numvec = 0
    numail = 0
    do i = 1, nbtyel
        nmel = zi(jnel-1+i)
!
        call jecroc(jexnum(cptlie, i))
        call jeecra(jexnum(cptlie, i), 'LONMAX', nmel+1)
!
        do j = 1, nmel
            numvec = numvec+1
            numail = numail+1
            zi(jdli+numvec-1) = lismai_nz(numail)
        end do
!
        numvec = numvec+1
        zi(jdli+numvec-1) = zi(jtyp-1+i)
!
    end do
!
    call jedetr('&&EXLIM1.TYPE_NOMBRE')

!   --  adaptation de la taille des grels et initialisation des elements
!   ---------------------------------------------------------------------
    call adalig(ligrel)
    call cormgi(base, ligrel)
    call initel(ligrel)
!
    AS_DEALLOCATE(vi=lismai_nz)
    call jedema()
end subroutine
