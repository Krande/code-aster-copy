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

subroutine cgnoxf(mofaz, iocc, nomaz, lisnoz, nbno)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cnocns.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: iocc, nbno
    character(len=*) :: mofaz, nomaz, lisnoz
!
!       CGNOXF -- TRAITEMENT DE L'OPTION FISS_XFEM
!                 DU MOT FACTEUR CREA_GROUP_NO DE
!                 LA COMMANDE DEFI_GROUP
!
!      CETTE FONCTIONNALITE PERMET DE CREER UN GROUP_NO CONSTITUE
!      DE TOUS LES NOEUDS DE TYPE XFEM DEFINI PAR L'UTILISATEUR.
!
! -------------------------------------------------------
!  MOFAZ         - IN    - K16  - : MOT FACTEUR 'CREA_GROUP_NO'
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  LISNOZ        - JXVAR - K24  - : NOM DE LA LISTE DE NOEUDS
!                                   DU TYPE XFEM DEMANDE PAR
!                                   L'UTILISATEUR
!  NBNO          - OUT   -  I   - : LONGUEUR DE CETTE LISTE
! -------------------------------------------------------
!
!.========================= DEBUT DES DECLARATIONS ====================
!
!
! --------- VARIABLES LOCALES ---------------------------
    integer(kind=8) :: ibid
    integer(kind=8) :: n1, ifiss, nfiss
    integer(kind=8) :: ino, valeno, nbnot
    integer(kind=8) :: idlist, jstno
    character(len=8) :: noma, nomnoe, fiss, nomafi
    character(len=8) :: nomagr, valk(2), ma
    character(len=16) :: motfac, typgrp
    character(len=19) :: stno, cnslt, cnsln
    character(len=24) :: stnot
    character(len=24) :: lisnoe
    aster_logical :: grille
    real(kind=8) :: rayon, dist
    character(len=8), pointer :: vfiss(:) => null()
    integer(kind=8), pointer :: noeu(:) => null()
    real(kind=8), pointer :: lsn(:) => null()
    real(kind=8), pointer :: lst(:) => null()
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
!
! --- INITIALISATIONS :
!     ================
    motfac = mofaz
    noma = nomaz
    lisnoe = lisnoz
    nbno = 0
!
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnot)
    AS_ALLOCATE(vi=noeu, size=nbnot)
!
! --  RECUPERATION DU TYPE GROUPE :
!     ============================
    call getvtx(motfac, 'TYPE_GROUP', iocc=iocc, scal=typgrp, nbret=n1)
!
! --  RECUPERATION DES NOMS DES FISSURES :
!     ===================================
    call getvid(motfac, 'FISSURE', iocc=iocc, nbval=0, nbret=nfiss)
    nfiss = -nfiss
    AS_ALLOCATE(vk8=vfiss, size=nfiss)
    call getvid(motfac, 'FISSURE', iocc=iocc, nbval=nfiss, vect=vfiss, &
                nbret=ibid)
!
! --- TYPE DE NOEUD = 'HEAVISIDE'
!     ============================
    if (typgrp .eq. 'HEAVISIDE') then
        do ifiss = 1, nfiss
            fiss = vfiss(ifiss)
            stno = fiss//'.STNO'
            call jeveuo(stno//'.VALE', 'L', jstno)
            do ino = 1, nbnot
                valeno = zi(jstno+ino-1)
                if (valeno .eq. 1) then
                    nbno = nbno+1
                    noeu(nbno) = ino
                end if
            end do
            call jedetr(stno)
        end do
!
! --- TYPE DE NOEUD = 'CRACKTIP'
!     ============================
    else if (typgrp .eq. 'CRACKTIP') then
        do ifiss = 1, nfiss
            fiss = vfiss(ifiss)
            stno = fiss//'.STNO'
            call jeveuo(stno//'.VALE', 'L', jstno)
            do ino = 1, nbnot
                valeno = abs(zi(jstno+ino-1))
                if (valeno .eq. 2) then
                    nbno = nbno+1
                    noeu(nbno) = ino
                end if
            end do
            call jedetr(stno)
        end do
!
! --- TYPE DE NOEUD = 'MIXTE'
!     ============================
    else if (typgrp .eq. 'MIXTE') then
        do ifiss = 1, nfiss
            fiss = vfiss(ifiss)
            stno = fiss//'.STNO'
            call jeveuo(stno//'.VALE', 'L', jstno)
            do ino = 1, nbnot
                valeno = zi(jstno+ino-1)
                if (valeno .eq. 3) then
                    nbno = nbno+1
                    noeu(nbno) = ino
                end if
            end do
        end do
!
! --- TYPE DE NOEUD = 'XFEM'
!     ============================
    else if (typgrp .eq. 'XFEM') then
        do ifiss = 1, nfiss
            fiss = vfiss(ifiss)
            stno = fiss//'.STNO'
            call jeveuo(stno//'.VALE', 'L', jstno)
            do ino = 1, nbnot
                valeno = zi(jstno+ino-1)
                if (valeno .ne. 0) then
                    nbno = nbno+1
                    noeu(nbno) = ino
                end if
            end do
            call jedetr(stno)
        end do
!
!
! --- TYPE DE NOEUD = 'TORE'
!     ============================
    else if ((typgrp .eq. 'TORE') .or. (typgrp .eq. 'ZONE_MAJ')) then
!
        cnslt = '&&CGNOXF.CNSLT'
        cnsln = '&&CGNOXF.CNSLN'
!
        do ifiss = 1, nfiss
            fiss = vfiss(ifiss)
!
!           CHECK IF THE LOCALISATION HAS BEEN USED
            stnot = fiss//'.PRO.RAYON_TORE'
            call jeexin(stnot, ibid)
            if ((ibid .gt. 0) .and. (typgrp .eq. 'TORE')) then
                typgrp = 'ZONE_MAJ'
                call utmess('A', 'XFEM2_92', sk=fiss)
            end if
!
            if (typgrp .eq. 'TORE') then
!
!              GET THE CRACK MESH
                call dismoi('NOM_MAILLA', fiss, 'FISS_XFEM', repk=nomafi)
!
                call getvr8(motfac, 'RAYON_TORE', iocc=1, scal=rayon, nbret=ibid)
                rayon = rayon**2
!
!              RETREIVE THE TWO LEVEL SETS
                call getvid(' ', 'MAILLAGE', scal=ma, nbret=ibid)
                if (ibid .eq. 0) then
                    call getvid(' ', 'GRILLE', scal=ma, nbret=ibid)
!                  CHECK FOR THE PRESENCE OF THE GRID
                    stnot = fiss//'.GRI.MAILLA'
                    call jeexin(stnot, ibid)
                    if (ibid .gt. 0) then
!                    GRID NAME
                        call jeveuo(stnot, 'L', ibid)
                        nomagr = zk8(ibid)
                        if (nomagr .ne. ma) then
                            call utmess('F', 'XFEM2_86')
                        end if
                    else
                        call utmess('F', 'XFEM2_86')
                    end if
                    call cnocns(fiss//'.GRI.LTNO', 'V', cnslt)
                    call cnocns(fiss//'.GRI.LNNO', 'V', cnsln)
                else
                    if (nomafi .ne. ma) then
                        call utmess('F', 'XFEM2_86')
                    end if
                    call cnocns(fiss//'.LTNO', 'V', cnslt)
                    call cnocns(fiss//'.LNNO', 'V', cnsln)
                end if
                call jeveuo(cnslt//'.CNSV', 'L', vr=lst)
                call jeveuo(cnsln//'.CNSV', 'L', vr=lsn)
!
                do ino = 1, nbnot
                    dist = lst(ino)**2+lsn(ino)**2
                    if (dist .le. rayon) then
                        nbno = nbno+1
                        noeu(nbno) = ino
                    end if
                end do
!
                call jedetr(cnslt)
                call jedetr(cnsln)
!
            end if
!
!
! --- TYPE DE NOEUD = 'ZONE_MAJ'
!     ============================
            if (typgrp .eq. 'ZONE_MAJ') then
!
!             GET THE CRACK MESH
                call dismoi('NOM_MAILLA', fiss, 'FISS_XFEM', repk=nomafi)
!
!             CHECK FOR THE PRESENCE OF THE GRID
                stnot = fiss//'.GRI.MAILLA'
                call jeexin(stnot, ibid)
                if (ibid .gt. 0) then
                    grille = .true.
!                GRID NAME
                    call jeveuo(stnot, 'L', ibid)
                    nomagr = zk8(ibid)
                else
                    grille = .false.
                end if
!
                if (noma .eq. nomafi) then
                    if (grille) then
                        stnot = fiss//'.PRO.NOEUD_PROJ'
                    else
                        stnot = fiss//'.PRO.NOEUD_TORE'
                    end if
                else if (grille .and. (noma .eq. nomagr)) then
                    stnot = fiss//'.PRO.NOEUD_TORE'
                else
                    valk(1) = nomafi
                    if (grille) then
                        valk(2) = nomagr
                    else
                        valk(2) = 'AUCUN'
                    end if
                    call utmess('F', 'XFEM2_96', nk=2, valk=valk)
                end if
!
                call jeexin(stnot, ibid)
                if (ibid .gt. 0) then
                    call jeveuo(stnot, 'L', jstno)
                    do ino = 1, nbnot
                        if (zl(jstno+ino-1)) then
                            nbno = nbno+1
                            noeu(nbno) = ino
                        end if
                    end do
                else
!                THE LOCALISATION HAS NOT BEEN USED. ZONE_MAJ IS
!                COINCIDENT WITH THE WHOLE MODEL.
                    do ino = 1, nbnot
                        nbno = nbno+1
                        noeu(nbno) = ino
                    end do
                end if
            end if
!
        end do
!
    else
        ASSERT(.false.)
    end if
!
    if (nbno .ne. 0) then
        call wkvect(lisnoe, 'V V I', nbno, idlist)
!
        do ino = 1, nbno
            zi(idlist+ino-1) = noeu(ino)
            nomnoe = int_to_char8(zi(idlist+ino-1))
        end do
    end if
!
! --- FIN
!     ===
!
! --- MENAGE
!
!
    AS_DEALLOCATE(vk8=vfiss)
    AS_DEALLOCATE(vi=noeu)
!
    call jedema()
!.============================ FIN DE LA ROUTINE ======================
end subroutine
