! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine peenca2(champ, long, vr, nbmail, nummai, ligrel, nbgr,&
                  ztot, option, nbproc, rang)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/digdel.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nbelem.h"
#include "asterfort/utmess.h"
#include "asterfort/vecint.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: champ
    integer :: long, nbmail, nummai(*), nbgr, option
    real(kind=8) :: vr(long), ztot
    character(len=19) :: ligrel
!     FAIRE DES OPERATIONS SUR UN CHAM_ELEM DE TYPE ENERGIE
!            (NOTION D'INTEGRALE DU CHAMP SUR LE MODELE)
!     ------------------------------------------------------------------
! IN  : CHAMP  : NOM DU CHAM_ELEM
! IN  : LONG   : LONGUEUR DU VECTEUR VR
! OUT : VR     : VECTEUR CONTENANT LES RESULATTS GLOBAUX
! IN  : NBMAIL : = 0 , CALCUL SUR TOUT LE CHAM_ELEM
!                SINON CALCUL SUR UN NOMBRE DE MAILLES
! IN  : NUMMAI : NUMEROS DES MAILLES
!
!-----------------------------------------------------------------------
    integer :: icoef, idecgr, iel, im, inum, j, k, nel, longt, mode
    integer :: ind, iaux, nbproc, rang, nbpas2, jldist2, iaux1, numim, ind1
    integer :: jad, jnbgr
    aster_logical :: ldist2, lluck
    real(kind=8), pointer :: celv(:) => null()
    integer, pointer :: liel(:) => null()
    integer, pointer :: celd(:) => null()
!-----------------------------------------------------------------------
! Init.
    call jemarq()
    call jeveuo(champ//'.CELD', 'L', vi=celd)
    icoef=max(1,celd(4))
    call jeveuo(champ//'.CELV', 'L', vr=celv)
    do k=1,long
      vr(k)=0.d0
    enddo
!
! Activation du parallelisme en espace pour construire la connectivite inverse
    ldist2=.true.
    if (nbmail.lt.nbproc) ldist2=.false.
!
!-----------------------------------------------------------------------------
! PREMIER PASSAGE: CALCUL CONNECTIVITE INVERSE ET CALCUL ENERGIE
!                  SI LDIST2=.TRUE. CALCUL CONNECTIVITE EN PARALLELE MPI
!-----------------------------------------------------------------------------
!
    if (option.lt.0) then
! Init.
      ind=-option
      call jeveuo(ligrel//'.LIEL', 'L', vi=liel)
! Objet pour stocker et mutualiser la connectivite inverse
      call jeveuo('&&PEEPOT_peenca', 'E', jad)
! Remplissage objet nbelem pour gagner le temps des jelira ds nbelem
      call wkvect('&&PEECA2_nbelem','V V I',nbgr,jnbgr)
      do j=1,nbgr
        zi(jnbgr+j-1)=nbelem(ligrel,j)
      enddo
!
      if (.not.ldist2) then
!
!-----------------------------------------------------------------------------
!     CALCUL CONNECTIVITE INVERSE (PART I) + ENERGIE (PART II) EN SEQUENTIEL
!-----------------------------------------------------------------------------
!PART I + II
        do 40 im = 1, nbmail
          inum = 0
          numim=nummai(im)
          do 42 j = 1, nbgr
            iaux=celd(4+j)+2
            mode=celd(iaux)
            nel=zi(jnbgr+j-1)
            if (mode .eq. 0) then
              inum = inum + nel + 1
              goto 42
            endif
            longt = digdel(mode)
            longt = longt * icoef
            idecgr=celd(iaux+6)
            do 44 k = 1, nel
              iel = liel(1+inum+k-1)
              if (iel .ne. numim) goto 44
              zi(jad+2*(ind-1))=j
              zi(jad+2*(ind-1)+1)=k
              ind=ind+1
              vr(1) = vr(1)+ celv(idecgr+(k-1)*longt)
              goto 40
44          continue
            inum = inum + nel + 1
42        continue
40      continue
!
      else
!
!-----------------------------------------------------------------------------
! PART I: CALCUL CONNECTIVITE INVERSE EN PARALLELE MPI (FILTRE VIA &PEECA2_vldist)
! PART II: CALCUL ENERGIE EN SEQUENTIEL GRACE A LA PART I
!-----------------------------------------------------------------------------
!PART I
! Filtre MPI type distribution de carte. Le reliquat est fait par tout le monde
! donc on ne le communique pas: pour les mailles de nbpas2*proc+1 jusqua nbmail
        call wkvect('&&PEECA2_vldist','V V I',nbmail,jldist2)
        call vecint(nbmail,rang,zi(jldist2))
! Determination du nbre de pas paralleles mpi: nbpas2
        nbpas2=nbmail/nbproc
        iaux1=0
        do k=1,nbpas2*nbproc
          if (iaux1.gt.(nbproc-1)) iaux1=0
          zi(jldist2+k-1)=iaux1
          iaux1=iaux1+1
        enddo
        lluck=.false.
        do 240 im = 1, nbmail
! Filtre MPI
          if (zi(jldist2+im-1).eq.rang) then
            numim=nummai(im)
            ind1=ind+im-1
            if (lluck) then
              k=k+1
              iel = liel(1+inum+k-1)
              if (iel .eq. numim) then
                zi(jad+2*(ind1-1))=j
                zi(jad+2*(ind1-1)+1)=k
                if (k.eq.nel) then
                  lluck=.false.
                  inum = 0
                else
                  lluck=.true.
                endif
                goto 239
              else
                lluck=.false.
                inum = 0
              endif
            else
              inum = 0
            endif
            do 242 j = 1, nbgr
              iaux=celd(4+j)+2
              mode=celd(iaux)
              nel=zi(jnbgr+j-1)
              if (mode .eq. 0) then
                inum = inum + nel + 1
                goto 242
              endif
              do 244 k = 1, nel
                iel = liel(1+inum+k-1)
                if (iel .ne. numim) goto 244
                if (k.eq.nel) then
                  lluck=.false.
                else
                  lluck=.true.
                endif
                zi(jad+2*(ind1-1))=j
                zi(jad+2*(ind1-1)+1)=k
                goto 239
244            continue
              inum = inum + nel + 1
242          continue
239          continue
! Fin du filtre
           endif
240     continue
! Nettoyage filtre et communication MPI associee (sauf le reliquat)
        call jedetr('&&PEECA2_vldist')
        iaux1=jad+2*(ind-1)
        call asmpi_comm_vect('MPI_SUM','I',nbval=2*nbpas2*nbproc,vi=zi(iaux1))
!
!PART II
! Calcul energie proprement dit via la connectivite nouvellement cree
        do im = 1, nbmail
          j=zi(jad+2*(ind-1))
          k=zi(jad+2*(ind-1)+1)
! On teste au cas ou afin de ne pas oublier des mailles.
! Visiblement cest acceptee pour cette option (cf. MA99 avec ssla200a)
          if (j <=0 .or. k<=0) then
              cycle
          endif
          ind=ind+1
          iaux=celd(4+j)+2
          mode=celd(iaux)
          longt = digdel(mode) * icoef
          idecgr=celd(iaux+6)
          vr(1) = vr(1)+ celv(idecgr+(k-1)*longt)
        enddo
      endif
! Nettoyage objet pour mutualisation nbelem
      call jedetr('&&PEECA2_nbelem')
! Renvoi a la routine appelante de l'increment
      option=-ind
!
!-----------------------------------------------------------------------------
! PASSAGES SUIVANTS: QUE CALCUL ENERGIE GRACE A LA CONNECTIVITE CALCULEE AU
!                  PASSAGE PRECEDENT
!-----------------------------------------------------------------------------
!
    else
! Init.
      call jeveuo('&&PEEPOT_peenca', 'L', jad)
      ind=option
! Calcul energie proprement dit via la connectivite cree au premier pas
      do 140 im = 1, nbmail
        j=zi(jad+2*(ind-1))
        k=zi(jad+2*(ind-1)+1)
! On teste au cas ou afin de ne pas oublier des mailles.
! Visiblement cest acceptee pour cette option (cf. MA99 avec ssla200a)
        if (j <=0 .or. k<=0) then
            cycle
        endif
        ind=ind+1
        iaux=celd(4+j)+2
        mode=celd(iaux)
        longt = digdel(mode) * icoef
        idecgr=celd(iaux+6)
        vr(1) = vr(1)+ celv(idecgr+(k-1)*longt)
140   continue
      option=ind
    endif
!
!-----------------------------------------------------------------------------
! MISE AU BON FORMAT DU CALCUL ENERGIE
!-----------------------------------------------------------------------------
!
    if (( vr(1).lt.r8prem() ) .and. ( ztot.lt.r8prem() )) then
      vr(2) = 0.0d0
    else
      vr(2) = 100.0d0 * vr(1) / ztot
    endif
!
    call jedema()
end subroutine
