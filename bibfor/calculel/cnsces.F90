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

subroutine cnsces(cnsz, typces, cesmoz, mnogaz, base, &
                  cesz)
! person_in_charge: jacques.pellet at edf.fr

    use compensated_ops_module, only: sum, dot_product

    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "MeshTypes_type.h"
!
    character(len=*) :: cnsz, cesz, base, cesmoz, typces, mnogaz
! ------------------------------------------------------------------
! But: transformer un cham_no_s en cham_elem_s
! ------------------------------------------------------------------
! Arguments:
! cnsz  in/jxin  k19 : sd cham_no_s a transformer
! typces in       k4  : type voulu pour le cham_elem_s
!                      /'ELEM' /'ELGA' /'ELNO'
! cesmoz in/jxin  k19 :  sd cham_elem_s "modele" pour cesz
!       si typces = 'ELEM' : cesmoz n'est pas utilise
!       si typces  ='ELGA' on se sert de cesmoz pour determiner
!          le nombre de points et de sous-points  du cham_elem_s
!       si typces  ='ELNO' on se sert de cesmoz pour determiner
!          le nombre de sous-points  du cham_elem_s
!
! mnogaz in/jxin  k19 : sd cham_elem_s contenant les matrices
!                       de passage noeud -> gauss.
!                       mnogaz n'est utilise que si elno->elga
! attention :  mnogaz est un cham_elem_s avec une convention
!              tres particuliere  (maille de reference)
!              (voir routine manopg.f)
!
! cesz   in/jxout k19 : sd cham_elem_s resultat
! base    in      k1  : base de creation pour cesz : g/v/l
!-----------------------------------------------------------------------
!
!  Principes retenus pour la conversion :
!  --------------------------------------
!
!  1) On ne traite que les cham_no_s reels (r8)
!  2)
!    si  typces='ELEM'
!       on affecte a la maille la moyenne arithmetique des noeuds
!    si  typces='ELNO'
!       on recopie la valeur du noeud global sur le noeud local
!    si  typces='ELGA'
!       on utilise la matrice de passage noeud -> gauss.
!
!  3) les eventuels sous-points portent tous les memes valeurs
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ima, ncmp, icmp, jcnsl
    integer(kind=8) :: jcesd, jcesl, nbma, iret, nbsp, nbno, ico
    integer(kind=8) :: iad, nbpt, ipt, ino, nuno, isp, nbpg2, nbno2, iad1
    integer(kind=8) :: ilcnx1, nbpg, ipg, imaref, sz
    integer(kind=8) :: mnogal, mnogad
    character(len=8) :: ma, nomgd
    character(len=3) :: tsca
    character(len=19) :: ces, cesmod, cns, mnoga
    real(kind=8) :: v
    real(kind=8), pointer :: cesv(:) => null()
    real(kind=8), pointer :: nmnogav(:) => null()
    real(kind=8) :: liv1(MT_NNOMAX), liv2(MT_NNOMAX)
    character(len=8), pointer :: cnsc(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
    character(len=8), pointer :: cesk(:) => null()
    integer(kind=8), pointer :: cemd(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
    integer(kind=8), pointer :: vnbpt(:) => null()
    integer(kind=8), pointer :: vnbsp(:) => null()
!------------------------------------------------------------------
    call jemarq()
!
!
    ces = cesz
    cesmod = cesmoz
    cns = cnsz
!
!
!   1- Recuperation d'informations dans cns :
!      ma    : nom du maillage
!      nomgd : nom de la grandeur
!      ncmp  : nombre de cmps dans cns
!   ------------------------------------------
    call jeveuo(cns//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(cns//'.CNSC', 'L', vk8=cnsc)
    call jeveuo(cns//'.CNSV', 'L', vr=cnsv)
    call jeveuo(cns//'.CNSL', 'L', jcnsl)
!
    ma = cnsk(1)
    nomgd = cnsk(2)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    ASSERT(tsca .eq. 'R')
    call jeveuo(ma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    call jelira(cns//'.CNSC', 'LONMAX', ncmp)
!
!
!   2. CALCUL DES OBJETS  '&&CNSCES.NBPT','&CNSCES.NBSP'
!   -----------------------------------------------------------------
    AS_ALLOCATE(vi=vnbpt, size=nbma)
    AS_ALLOCATE(vi=vnbsp, size=nbma)
!
!
!   -- PAR DEFAUT : NBSP=1
    do ima = 1, nbma
        vnbsp(ima) = 1
    end do
!
    call exisd('CHAM_ELEM_S', cesmod, iret)
    ASSERT((typces .ne. 'ELGA') .or. (iret .gt. 0))
!
    if (iret .gt. 0) then
        call jeveuo(cesmod//'.CESD', 'L', vi=cemd)
        do ima = 1, nbma
            vnbpt(ima) = cemd(5+4*(ima-1)+1)
            vnbsp(ima) = cemd(5+4*(ima-1)+2)
        end do
    end if
!
    if (typces .eq. 'ELEM') then
        do ima = 1, nbma
            vnbpt(ima) = 1
        end do
    else if (typces .eq. 'ELNO') then
        do ima = 1, nbma
            vnbpt(ima) = zi(ilcnx1+ima)-zi(ilcnx1+ima-1)
        end do
    else if (typces .eq. 'ELGA') then
!       DEJA FAIT GRACE A CESMOD
    else
        ASSERT(.false.)
    end if
!
!
!   5- CREATION DE CES :
!   ---------------------------------------
    call cescre(base, ces, typces, ma, nomgd, &
                ncmp, cnsc, vnbpt, vnbsp, [-ncmp])
!
    call jeveuo(ces//'.CESD', 'L', jcesd)
    call jeveuo(ces//'.CESV', 'E', vr=cesv)
    call jeveuo(ces//'.CESL', 'E', jcesl)
!
!
!
!   6- REMPLISSAGE DES OBJETS .CESL ET .CESV :
!   ------------------------------------------
!
!
    if (typces .eq. 'ELEM') then
!     --------------------------
        do ima = 1, nbma
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            nbsp = zi(jcesd-1+5+4*(ima-1)+2)
            nbno = zi(ilcnx1+ima)-zi(ilcnx1-1+ima)
            do icmp = 1, ncmp
!
!           - ON VERIFIE QUE TOUS LES NOEUDS PORTENT BIEN LA CMP :
                ico = 0
                do ino = 1, nbno
                    nuno = connex(1+zi(ilcnx1-1+ima)-2+ino)
                    if (zl(jcnsl-1+(nuno-1)*ncmp+icmp)) ico = ico+1
                end do
                if (ico .ne. nbno) goto 90
!
!         -- CALCUL DE LA MOYENNE ARITHMETIQUE :
                liv1(:) = 0.d0
                sz = 1
                do ino = 1, nbno
                    nuno = connex(1+zi(ilcnx1-1+ima)-2+ino)
                    if (zl(jcnsl-1+(nuno-1)*ncmp+icmp)) then
                        liv1(sz) = cnsv((nuno-1)*ncmp+icmp)
                        sz = sz+1
                    end if
                end do
!
                v = sum(liv1)/nbno
!
                do ipt = 1, nbpt
                    do isp = 1, nbsp
                        call cesexi('C', jcesd, jcesl, ima, ipt, &
                                    isp, icmp, iad)
                        ASSERT(iad .lt. 0)
                        zl(jcesl-1-iad) = .true.
                        cesv(1-1-iad) = v
                    end do
                end do
90              continue
            end do
        end do
!
!
    else if (typces .eq. 'ELNO') then
!   --------------------------------
        do ima = 1, nbma
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            nbsp = zi(jcesd-1+5+4*(ima-1)+2)
            nbno = zi(ilcnx1+ima)-zi(ilcnx1-1+ima)
            ASSERT(nbno .eq. nbpt)
!
            do icmp = 1, ncmp
!
!           - ON VERIFIE QUE TOUS LES NOEUDS PORTENT BIEN LA CMP :
                ico = 0
                do ino = 1, nbno
                    nuno = connex(1+zi(ilcnx1-1+ima)-2+ino)
                    if (zl(jcnsl-1+(nuno-1)*ncmp+icmp)) ico = ico+1
                end do
                if (ico .ne. nbno) goto 140
!
                do ino = 1, nbno
                    nuno = connex(1+zi(ilcnx1-1+ima)-2+ino)
                    if (.not. zl(jcnsl-1+(nuno-1)*ncmp+icmp)) goto 130
                    v = cnsv((nuno-1)*ncmp+icmp)
                    do isp = 1, nbsp
                        call cesexi('C', jcesd, jcesl, ima, ino, &
                                    isp, icmp, iad)
                        ASSERT(iad .lt. 0)
                        zl(jcesl-1-iad) = .true.
                        cesv(1-1-iad) = v
                    end do
130                 continue
                end do
140             continue
            end do
        end do
!
!
    else if (typces .eq. 'ELGA') then
!     --------------------------
        mnoga = mnogaz
        call jeveuo(mnoga//'.CESK', 'L', vk8=cesk)
        call jeveuo(mnoga//'.CESD', 'L', mnogad)
        call jeveuo(mnoga//'.CESL', 'L', mnogal)
        call jeveuo(mnoga//'.CESV', 'L', vr=nmnogav)
        ASSERT(cesk(1) .eq. ma)
!
        do ima = 1, nbma
            call cesexi('C', mnogad, mnogal, ima, 1, &
                        1, 1, iad)
            if (iad .le. 0) goto 210
            if (nint(nmnogav(iad)) .gt. 0) then
                imaref = ima
            else
                imaref = -nint(nmnogav(iad))
            end if
            call cesexi('C', mnogad, mnogal, imaref, 1, &
                        1, 1, iad)
            if (iad .le. 0) goto 210
!
            nbno2 = nint(nmnogav(iad))
            nbpg2 = nint(nmnogav(iad+1))
!
            nbpg = zi(jcesd-1+5+4*(ima-1)+1)
            nbsp = zi(jcesd-1+5+4*(ima-1)+2)
            nbno = zi(ilcnx1+ima)-zi(ilcnx1-1+ima)
            if (nbno .ne. nbno2 .and. cnsz .eq. '&&VRCIN1.CNS1') nbno = nbno2
            ASSERT(nbno .eq. nbno2)
            ASSERT(nbpg .eq. nbpg2)
            liv1(:) = 0.d0
            liv2(:) = 0.d0
!
            do icmp = 1, ncmp
!
!           - ON VERIFIE QUE TOUS LES NOEUDS PORTENT BIEN LA CMP :
                ico = 0
                do ino = 1, nbno
                    nuno = connex(1+zi(ilcnx1-1+ima)-2+ino)
                    if (zl(jcnsl-1+(nuno-1)*ncmp+icmp)) ico = ico+1
                end do
                if (ico .ne. nbno) goto 200
!
                do ipg = 1, nbpg
                    do ino = 1, nbno
                        nuno = connex(1+zi(ilcnx1-1+ima)-2+ino)
                        liv1(ino) = cnsv((nuno-1)*ncmp+icmp)
                        liv2(ino) = nmnogav(iad+1+nbno*(ipg-1)+ino)
                    end do
!
                    v = dot_product(liv1, liv2)
!
                    do isp = 1, nbsp
                        call cesexi('C', jcesd, jcesl, ima, ipg, &
                                    isp, icmp, iad1)
                        ASSERT(iad1 .lt. 0)
                        zl(jcesl-1-iad1) = .true.
                        cesv(1-1-iad1) = v
                    end do
                end do
!
200             continue
            end do
210         continue
        end do
!
    end if
!
!
!     7- MENAGE :
!     -----------
    AS_DEALLOCATE(vi=vnbpt)
    AS_DEALLOCATE(vi=vnbsp)
!
    call jedema()
end subroutine
