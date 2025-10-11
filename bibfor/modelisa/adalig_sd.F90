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

subroutine adalig_sd(ligr, partsd, ntliel, nbtype, clas, teut, nteut)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jevtbl.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/sdpart.h"
#include "asterfort/wkvect.h"
!
    character(len=19), intent(in) :: ligr
    character(len=8), intent(in) :: partsd
    character(len=24), intent(in) :: ntliel
    integer(kind=8), intent(in) :: nbtype
    character(len=1), intent(in) :: clas
    integer(kind=8), pointer :: nteut(:)
    integer(kind=8), pointer :: teut(:)
!----------------------------------------------------------------------
! But: Reorganiser la collection .LIEL de ligrz afin de regrouper
!      les elements de meme TYPE_ELEM dans un meme GREL
!      et en respectant la partition partsd.
!
! De plus, on veut :
!   * Limiter la taille des GRELS (pas plus de nelmx elements)
!   * Faire en sorte que les GRELS respectent partsd :
!     (tous les elements d'un GREL appartiennent aux sous-domaines traites par le
!      processeur qui traitera ce GREL)
!     * Pour chaque TYPE_ELEM :
!       * On decoupe le paquet d'elements en un nombre de grels multiple de nbproc.
!       * le GREL kgrel ne contient que des elements des sous-domaines affectes au
!           processeur kproc [0, ..., nbproc-1] avec : mod(kgrel,nbproc)=kproc
!
!    Remarques :
!      * L'equilibrage n'est pas toujours optimal : cela depend de la "qualite" de
!        la partition. Il est preferable que nbSD soit un multiple de nbproc,
!        mais ca ne suffit pas pour assurer un equlibrage "parfait"
!      * Le desequilibrage conduit parfois a creer des GRELS vides.
!
!
! Arguments d'entree:
!     ligr   (o) : nom du ligrel
!     partsd (o) : nom de la partsd
!----------------------------------------------------------------------
    character(len=24) :: liel
    character(len=8) :: noma
    integer(kind=8) :: i, jliel
    integer(kind=8) ::  igrel, itype, nbma, ksd, kproc
    integer(kind=8) :: nbel, nbgrel, kgre1, nbgre1
    integer(kind=8) :: nbelmx, npaq, nbelgr, nel1, nel2, nel3, decal
    integer(kind=8) :: ktype, lont, nbgrel_av, rang, nbproc
    integer(kind=8) :: nbsd, igr, ktyp1, ktyp2, iel, numa, ngrmx, ngr1, ipaq
    integer(kind=8), pointer :: traite_par(:) => null()
    integer(kind=8), pointer :: nbel1(:) => null()
    integer(kind=8), pointer :: gteut(:) => null()
    integer(kind=8), pointer :: nbelgrel(:) => null()
    integer(kind=8), pointer :: utilise(:) => null()
    integer(kind=8), pointer :: utilise_1(:) => null()
    integer(kind=8), pointer :: utilise_2(:) => null()
    integer(kind=8), pointer :: ordre_stockage(:) => null()
    mpi_int :: mrank, msize

    integer(kind=8), pointer :: feta(:) => null()
    integer(kind=8), pointer :: tliel(:) => null()
    integer(kind=8), pointer :: lctliel(:) => null()
    integer(kind=8), pointer :: sdloc(:) => null()
#define numail(igr,iel) tliel(lctliel(igr)+iel-1)
!----------------------------------------------------------------------

    call jemarq()

    liel = ligr//'.LIEL'
    call dismoi('NOM_MAILLA', ligr, 'LIGREL', repk=noma)
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
    call jelira(ntliel, 'NMAXOC', nbgrel_av)
    call jeveuo(ntliel, 'L', vi=tliel)
    call jeveuo(jexatr(ntliel, 'LONCUM'), 'L', vi=lctliel)
    call jelira(partsd//'.FETA', 'NMAXOC', nbsd)

    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)

!   -- Calcul du vecteur traite_par:
!      traite_par(ima) = kproc
!   -----------------------------------------------
!   -- on repartit les sous-domaines
    call wkvect('&&ADALIG_SD.PART.SD', 'V V I', nbsd, vi=sdloc)
    call sdpart(nbsd, 0, sdloc)

    AS_ALLOCATE(vi=traite_par, size=nbma)
    traite_par(:) = -99
    do ksd = 1, nbsd
        if (sdloc(ksd) .eq. 1) then
            call jeveuo(jexnum(partsd//'.FETA', ksd), 'L', vi=feta)
            if (rang .eq. 0) then
                do i = 1, size(feta)
                    traite_par(feta(i)) = nbproc-1
                end do
            else
                do i = 1, size(feta)
                    traite_par(feta(i)) = rang-1
                end do
            end if
        end if
    end do
    call asmpi_comm_vect('MPI_MAX', 'I', nbval=nbma, vi=traite_par)
    call jedetr('&&ADALIG_SD.PART.SD')

!   -- Calcul du vecteur nbel1:
!      nbel1((ktype-1)*nbproc+kproc+1) = nel
!      nel est le nombre d'elements de type ktype qui seront
!      calcules par le proc kproc
!   -----------------------------------------------
    AS_ALLOCATE(vi=nbel1, size=nbtype*nbproc)
    nbel1(:) = 0
    do igr = 1, nbgrel_av
        nbelgr = lctliel(igr+1)-lctliel(igr)-1
        ktyp1 = tliel(lctliel(igr)-1+nbelgr+1)
        ktype = 0
        do ktyp2 = 1, nbtype
            if (ktyp1 .eq. teut(ktyp2)) then
                ktype = ktyp2
                exit
            end if
        end do
        ASSERT(ktype .gt. 0)
        do iel = 1, nbelgr
            numa = numail(igr, iel)
            kproc = traite_par(numa)
            ASSERT(kproc .ge. 0)
            nbel1((ktype-1)*nbproc+kproc+1) = nbel1((ktype-1)*nbproc+kproc+1)+1
        end do
    end do

!   -- Calcul du nombre de grels du nouveau .LIEL
!      et de la dimension totale de la collection
!   -----------------------------------------------
!   gteut : nombre de grels du ligrel (par type_elem)
    AS_ALLOCATE(vi=gteut, size=nbtype)
    lont = 0
    nbgrel = 0
    nbelmx = int(jevtbl('TAILLE_GROUP_ELEM'))
    do ktype = 1, nbtype
        ngrmx = 0
        do kproc = 0, nbproc-1
            nel1 = nbel1((ktype-1)*nbproc+kproc+1)
            ngr1 = nel1/nbelmx
            if (mod(nel1, nbelmx) .gt. 0) ngr1 = ngr1+1
            ngrmx = max(ngrmx, ngr1)
        end do
        nbgre1 = nbproc*ngrmx
        gteut(ktype) = nbgre1
        nbgrel = nbgrel+nbgre1
        nbel = nteut(ktype)
        lont = lont+nbel+nbgre1
    end do

!   -- Calcul du nombre d'elements des GRELS du nouveau .LIEL
!   -------------------------------------------------------------
    AS_ALLOCATE(vi=nbelgrel, size=nbgrel)
    igrel = 0
    do ktype = 1, nbtype
        nbgre1 = gteut(ktype)
        npaq = nbgre1/nbproc
        do ipaq = 1, npaq
            do kproc = 0, nbproc-1
                igrel = igrel+1
                nel1 = nbel1((ktype-1)*nbproc+kproc+1)
                nel2 = ipaq*nbelmx
                nel3 = (ipaq-1)*nbelmx
                if (nel2 .le. nel1) then
                    nbelgrel(igrel) = nbelmx
                else
                    if (nel3 .gt. nel1) then
                        nbelgrel(igrel) = 0
                    else
                        nbelgrel(igrel) = nel1-nel3
                    end if
                end if
            end do
        end do
    end do

!   -- Allocation du nouveau .LIEL
!   ------------------------------------
    call jecrec(liel, clas//' V I', 'NU', 'CONTIG', 'VARIABLE', nbgrel)
    call jeecra(liel, 'LONT', lont)
    igrel = 0
    do ktype = 1, nbtype
        itype = teut(ktype)
        nbgre1 = gteut(ktype)
        do kgre1 = 1, nbgre1
            igrel = igrel+1
            nbelgr = nbelgrel(igrel)
            call jecroc(jexnum(liel, igrel))
            call jeecra(jexnum(liel, igrel), 'LONMAX', nbelgr+1)
            call jeveuo(jexnum(liel, igrel), 'E', jliel)
            zi(jliel+nbelgr) = itype
        end do
    end do
    ASSERT(nbgrel .eq. igrel)

!   -- Remplissage des nouveaux GREL
!   ----------------------------------
    AS_ALLOCATE(vi=ordre_stockage, size=nbma)
    AS_ALLOCATE(vi=utilise, size=nbproc)
    AS_ALLOCATE(vi=utilise_2, size=nbproc+1)
    AS_ALLOCATE(vi=utilise_1, size=nbproc)
    igrel = 0
    do ktyp2 = 1, nbtype
        nbgre1 = gteut(ktyp2)
        npaq = nbgre1/nbproc

!       -- on remplit des objets qui facilitent le remplissage des grels :
        utilise(:) = 0
        do igr = 1, nbgrel_av
            nbelgr = lctliel(igr+1)-lctliel(igr)-1
            ktyp1 = tliel(lctliel(igr)-1+nbelgr+1)
            if (ktyp1 .ne. teut(ktyp2)) cycle
            do iel = 1, nbelgr
                numa = numail(igr, iel)
                kproc = traite_par(numa)
                utilise(kproc+1) = utilise(kproc+1)+1
            end do
        end do
        utilise_2(1) = 0
        do kproc = 0, nbproc-1
            utilise_2(kproc+2) = utilise_2(kproc+1)+utilise(kproc+1)
        end do

        ordre_stockage(:) = 0
        utilise_1(:) = 0
        do igr = 1, nbgrel_av
            nbelgr = lctliel(igr+1)-lctliel(igr)-1
            ktyp1 = tliel(lctliel(igr)-1+nbelgr+1)
            if (ktyp1 .ne. teut(ktyp2)) cycle
            do iel = 1, nbelgr
                numa = numail(igr, iel)
                kproc = traite_par(numa)
                utilise_1(kproc+1) = utilise_1(kproc+1)+1
                ordre_stockage(utilise_2(kproc+1)+utilise_1(kproc+1)) = numa
            end do
        end do

        do ipaq = 1, npaq
            do kproc = 0, nbproc-1
                igrel = igrel+1
                call jeveuo(jexnum(liel, igrel), 'E', jliel)

                nbelgr = nbelgrel(igrel)
                decal = max(0, (ipaq-1)*nbelmx)
                do iel = 1, nbelgr
                    numa = ordre_stockage(utilise_2(kproc+1)+decal+iel)
                    zi(jliel-1+iel) = numa
                end do
            end do
        end do
    end do
    ASSERT(nbgrel .eq. igrel)

    AS_DEALLOCATE(vi=traite_par)
    AS_DEALLOCATE(vi=nbel1)
    AS_DEALLOCATE(vi=gteut)
    AS_DEALLOCATE(vi=nbelgrel)
    AS_DEALLOCATE(vi=ordre_stockage)
    AS_DEALLOCATE(vi=utilise)
    AS_DEALLOCATE(vi=utilise_1)
    AS_DEALLOCATE(vi=utilise_2)

    call jedema()
end subroutine
