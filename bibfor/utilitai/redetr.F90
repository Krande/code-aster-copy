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
subroutine redetr(matelz)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/zerosd.h"
!
    character(len=*) :: matelz
! person_in_charge: jacques.pellet at edf.fr
!
!      BUT: DETRUIRE DANS LE MATR_ELEM  MATELZ LES RESUELEM NULS
!           MAIS EN PRENANT GARDE QU'IL RESTE QUELQUE CHOSE
!
!     IN/OUT  : MATELZ = NOM DE LA SD MATR_ELEM A NETTOYER
!
!
    integer(kind=8) :: iret1, iexi, iexiav
    integer(kind=8) :: izero, ico, k, nb1, nbdet, nb1av
    aster_logical :: ldetr
    character(len=3) :: kret
    character(len=8) :: noma
    character(len=19) :: matele, resuel
    character(len=24), pointer :: relr(:) => null()
    integer(kind=8), pointer :: adetr(:) => null()
    character(len=24), pointer :: tempor(:) => null()
!
    call jemarq()
!
    matele = matelz
    ldetr = .false.
    call dismoi('NOM_MAILLA', matele, 'MATR_ELEM', repk=noma)
    call dismoi('PARALLEL_MESH', noma, 'MAILLAGE', repk=kret)
!
!     -- SI LE MATR_ELEM NE CONTIENT QUE DES MACRO-ELEMENTS,
!        L'OBJET .RELR N'EXISTE PAS ET IL N'Y A RIEN A FAIRE :
    call jeexin(matele//'.RELR', iexi)
    iexi = min(1, abs(iexi))
    if (kret .eq. 'NON') then
        iexiav = iexi
        call asmpi_comm_vect('MPI_MAX', 'I', sci=iexi)
        iexi = min(1, abs(iexi))
        ASSERT(iexi .eq. iexiav)
    end if
    if (iexi .eq. 0) goto 60
!
    call jeveuo(matele//'.RELR', 'E', vk24=relr)
    call jelira(matele//'.RELR', 'LONUTI', nb1)
!
!     -- LE MATR_ELEM DOIT CONTENIR LE MEME NOMBRE DE RESUELEM
!        SUR TOUS LES PROCESSEURS :
    if (kret .eq. 'NON') then
        nb1av = nb1
        call asmpi_comm_vect('MPI_MAX', 'I', sci=nb1)
        ASSERT(nb1 .eq. nb1av)
    end if
!
!     -- SI LE MATR_ELEM NE CONTIENT QU'1 RESUELEM OU AUCUN,
!        IL NE FAUT RIEN DETRUIRE
    if (nb1 .eq. 1 .or. nb1 .eq. 0) goto 60
!
!     -- CREATION DES OBJETS TEMPORAIRES DE TRAVAIL
!        ET DU BOOLEEN POUR DESTRUCTION A LA SORTIE
    ldetr = .true.
    AS_ALLOCATE(vk24=tempor, size=nb1)
    AS_ALLOCATE(vi=adetr, size=nb1)
!
!     -- ON EXAMINE LES RESUELEM CANDIDATS A LA DESTRUCTION :
!        ADETR(K)=1 : LE NOM EST ' '
!        ADETR(K)=2 : LE NOM EST NON ' ' MAIS LA SD N'EXISTE PAS
!        ADETR(K)=3 : LA SD EXISTE ET ELLE EST NULLE
!        ADETR(K)=0 : LA SD EXISTE ET ELLE EST NON NULLE
!     REMARQUE : LES CAS 1 ET 2 N'EXISTENT PAS ENCORE
!                J'ESPERE QU'ILS N'ARRIVERONT JAMAIS
    do k = 1, nb1
        adetr(k) = 0
        resuel = relr(k) (1:19)
        if (resuel .eq. ' ') then
            ASSERT(.false.)
            adetr(k) = 1
            goto 10
        end if
!
!       -- EXISTENCE DU RESU_ELEM ?
        call exisd('RESUELEM', resuel, iret1)
        if (iret1 .eq. 0) then
            adetr(k) = 2
            ASSERT(.false.)
            goto 10
        end if
!
!
!       -- SI LE RESU_ELEM EST NUL SUR TOUS LES PROCS,
!          ON PEUT LE DETRUIRE:
        izero = 1
        if (zerosd('RESUELEM', resuel)) izero = 0
        if (kret .eq. 'NON') then
            call asmpi_comm_vect('MPI_MAX', 'I', sci=izero)
        end if
        if (izero .eq. 0) then
            adetr(k) = 3
        else
            adetr(k) = 0
        end if
10      continue
    end do
!
!
!     -- ON COMPTE LES RESUELEM A DETRUIRE :
    nbdet = 0
    do k = 1, nb1
        if (adetr(k) .eq. 3) nbdet = nbdet+1
    end do
    if (nbdet .eq. 0) goto 60
!
!     -- ON DETRUIT LES RESULEM NULS (ON EN GARDE AU MOINS 1) :
!        ON PART DE LA FIN CAR LA MATRICE NON SYMETRIQUE EST
!        EN GENERAL STOCKEE APRES LA SYMETRIQUE
    nbdet = min(nbdet, nb1-1)
    ico = 0
    do k = nb1, 1, -1
        if (adetr(k) .eq. 3) then
            ico = ico+1
            if (ico .gt. nbdet) goto 31
            resuel = relr(k) (1:19)
            call detrsd('RESUELEM', resuel)
            relr(k) = ' '
        end if
    end do
31  continue
!
!     -- ON COMPACTE LE MATR_ELEM POUR QUE TOUS SES RESUELEM
!        SOIENT "VRAIS"
    ico = 0
    do k = 1, nb1
        resuel = relr(k) (1:19)
        if (resuel .ne. ' ') then
            ico = ico+1
            tempor(ico) = resuel
        end if
    end do
    ASSERT(ico .gt. 0)
!
    call jeecra(matele//'.RELR', 'LONUTI', ico)
    do k = 1, ico
        relr(k) = tempor(k)
    end do
!
60  continue
!
!     -- DESTRUCTION DES OBJETS TEMPORAIRES SI BESOIN
    if (ldetr) then
        AS_DEALLOCATE(vk24=tempor)
        AS_DEALLOCATE(vi=adetr)
    end if
!
    call jedema()
!
end subroutine
