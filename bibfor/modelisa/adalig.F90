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

subroutine adalig(ligrz, partsdz)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/adalig_sd.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jevtbl.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/gettco.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=*), intent(in) :: ligrz
    character(len=*), intent(in), optional :: partsdz
!----------------------------------------------------------------------
! But: Reorganiser la collection .LIEL de ligrz afin de regrouper
!      les elements de meme TYPE_ELEM dans un meme GREL.
!
! On veut :
!   * Limiter la taille des GRELS (pas plus de nelmx elements)
!   * Faire en sorte que l'equilibrage soit bon pour DISTRIBUTION / METHODE='GROUP_ELEM' :
!     * Pour chaque TYPE_ELEM :
!       On decoupe le paquet d'elements en un nombre de grels multiple de nbproc.
!     * Si partsd est n'est pas fourni :
!       L'equilibrage est presque parfait :
!       Les GRELS ont tous le meme nombre d'elements (a 1 pres)
!
!     * Si partsd est fourni:
!       * On ajoute une nouvelle contrainte pour les GRELS :
!         * le GREL kgrel ne contient que des elements des sous-domaines affectes au
!           processeur kproc [0, ..., nbproc-1] avec : mod(kgrel,nbproc)=kproc
!
!
! Arguments d'entree:
!     ligrz  (o) : nom du ligrel
!     partsd (f) : nom de la partsd
!----------------------------------------------------------------------

    character(len=19) :: ligr, partsd
    character(len=1) :: clas
    character(len=8) :: noma
    character(len=16) :: typsd
    character(len=24) :: liel, tliel
    integer(kind=8) :: i, iret, nbtg, iad, iadp, iadt, iadtp, jliel, jtlie2
    integer(kind=8) ::  jtliel, igrel, itype, j, jtype
    integer(kind=8) :: nbel, nbelem, nbg, nbgrel, nbtype, nel, nelem, ntot
    integer(kind=8) :: nbelmx, rang, nbproc, np1, nspaq, nbelgr, igre2
    integer(kind=8) :: ktype, k, nbelgv, lont
    integer(kind=8), pointer :: gteut(:) => null()
    integer(kind=8), pointer :: nteut(:) => null()
    integer(kind=8), pointer :: teut(:) => null()
    mpi_int :: mrank, msize
    aster_logical :: lhpc
!----------------------------------------------------------------------

    call jemarq()

    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    ligr = ligrz
    liel = ligr//'.LIEL'
    typsd = '****'
    call jeexin(liel, iret)
    if (iret .eq. 0) then
        goto 999
    end if
    call jelira(liel, 'NUTIOC', nbgrel)
    if (nbgrel .eq. 0) then
        goto 999
    end if

!   -- le maillage est-il distribué pour être en mode HPC
!   ----------------------------------------------------
    lhpc = ASTER_FALSE
    call dismoi('NOM_MAILLA', ligr, 'LIGREL', repk=noma)
    call gettco(noma, typsd)
    if (typsd .eq. 'MAILLAGE_P') then
        lhpc = ASTER_TRUE
    end if

!   -- recopie de liel dans tliel et destruction de liel
!   ----------------------------------------------------
    tliel = '&&ADALIG.LIEL'
    call jelira(liel, 'CLAS', cval=clas)
    call jedupo(liel, 'V', tliel, .true._1)
    call jedetr(liel)
    call jelira(tliel, 'NMAXOC', nbtg)

!   -- Calcul de 2 vecteurs de travail (sur-dimensionnes) :
!     teut  : liste des type_elem utilises dans le ligrel
!     nteut : nombre total d'elements du ligrel (par type_elem)
!   ----------------------------------------------------
    AS_ALLOCATE(vi=teut, size=nbtg)
    AS_ALLOCATE(vi=nteut, size=nbtg)

    call jeveuo(tliel, 'L', jtliel)
    call jeveuo(jexatr(tliel, 'LONCUM'), 'L', jtlie2)
    iad = zi(jtlie2)
    nbtype = 0
    do i = 1, nbtg
        iadp = zi(jtlie2+i)
        nbelem = iadp-iad-1
        iad = iadp
        if (nbelem .gt. 0) then
            itype = zi(jtliel-1+iadp-1)
            do j = 1, nbtype
                if (itype .eq. teut(j)) then
                    nteut(j) = nteut(j)+nbelem
                    goto 1
                end if
            end do
            nbtype = nbtype+1
            teut(nbtype) = itype
            nteut(nbtype) = nbelem
        end if
1       continue
    end do

!   -- si partsd est fourni, il faut utiliser un autre algorithme :
!   ---------------------------------------------------------------
    if (present(partsdz)) then
        partsd = partsdz
        call adalig_sd(ligr, partsd, tliel, nbtype, clas, teut, nteut)
        goto 998
    end if

!   -- Calcul du nombre de grels du nouveau .LIEL
!      et de la dimension totale de la collection
!   -----------------------------------------------
!   gteut : nombre de grels du ligrel (par type_elem)
    AS_ALLOCATE(vi=gteut, size=nbtg)

    lont = 0
    nbgrel = 0
    nbelmx = int(jevtbl('TAILLE_GROUP_ELEM'))
    do ktype = 1, nbtype
        nbel = nteut(ktype)
        if (lhpc) then
            nspaq = nbel/nbelmx
            ASSERT((nspaq*nbelmx .le. nbel))
            if (nspaq*nbelmx .lt. nbel) nspaq = nspaq+1
            nbg = nspaq
        else
            nspaq = (nbel/nbproc)/nbelmx
            ASSERT((nspaq*nbproc*nbelmx .le. nbel))
            if (nspaq*nbproc*nbelmx .lt. nbel) nspaq = nspaq+1
            nbg = nspaq*nbproc
        end if
!
        gteut(ktype) = nbg
        nbgrel = nbgrel+nbg
        lont = lont+nbel+nbg
    end do
    if (.not. lhpc) then
        ASSERT((nbgrel/nbproc)*nbproc .eq. nbgrel)
    end if

!   -- Allocation du nouveau .LIEL
!   ------------------------------------
    call jecrec(liel, clas//' V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbgrel)
    call jeecra(liel, 'LONT', lont)
    igrel = 0
    do ktype = 1, nbtype
        itype = teut(ktype)
        ntot = nteut(ktype)
        nbg = gteut(ktype)
        nbelgr = ntot/nbg
!       -- attention : pour les petits GREL, il peut arriver que
!          nbelgr=0 (si nbproc > ntot)
!          il y aura alors des GREL vides

!       -- le nombre d'elements par GREL sera nbelgr ou nbelgr+1
!       -- les np1 1ers GREL de ktype auront nbelgr+1 elements
!          les autres auront nbelgr elements
        np1 = ntot-nbg*nbelgr
        ASSERT(np1 .lt. nbg)
        do k = 1, nbg
            nbelgv = nbelgr
            if (k .le. np1) nbelgv = nbelgv+1
            call jecroc(jexnum(liel, igrel+k))
            call jeecra(jexnum(liel, igrel+k), 'LONMAX', nbelgv+1)
            call jeveuo(jexnum(liel, igrel+k), 'E', jliel)
            zi(jliel+nbelgv) = itype
        end do
        igrel = igrel+nbg
    end do
    ASSERT(nbgrel .eq. igrel)

!   -- Remplissage des nouveaux GREL
!   ----------------------------------
    igrel = 0
    do ktype = 1, nbtype
        itype = teut(ktype)
        ntot = nteut(ktype)
        nbg = gteut(ktype)
        nbelgr = ntot/nbg
        np1 = ntot-nbg*nbelgr
        ASSERT(np1 .lt. nbg)

        igre2 = 1
        nbelgv = nbelgr
        if (igre2 .le. np1) nbelgv = nbelgv+1
        call jeveuo(jexnum(liel, igrel+igre2), 'E', jliel)
        nelem = 0

!       -- On remplit les nouveaux GREL avec les elements du bon type
        iadt = zi(jtlie2)
        do j = 1, nbtg
            iadtp = zi(jtlie2+j)
            jtype = zi(jtliel-2+iadtp)
            if (jtype .eq. itype) then
                nel = iadtp-iadt-1
                do k = 1, nel
!                   -- Il faut changer de GREL :
                    if (nelem .ge. nbelgv) then
                        igre2 = igre2+1
                        nbelgv = nbelgr
                        if (igre2 .le. np1) nbelgv = nbelgv+1
                        call jeveuo(jexnum(liel, igrel+igre2), 'E', jliel)
                        nelem = 0
                    end if
                    nelem = nelem+1
                    zi(jliel-1+nelem) = zi(jtliel-1+iadt+k-1)
                end do
            end if
            iadt = iadtp
        end do
        ASSERT(igre2 .le. nbg)
        ASSERT(nelem .eq. nbelgv)
        igrel = igrel+nbg
    end do
    ASSERT(igrel .eq. nbgrel)

!   -- Destruction des objets de travail
    AS_DEALLOCATE(vi=gteut)

998 continue
    AS_DEALLOCATE(vi=teut)
    AS_DEALLOCATE(vi=nteut)
    call jedetr(tliel)

999 continue
    call jedema()
end subroutine
