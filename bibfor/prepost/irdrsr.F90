! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine irdrsr(ifi, nbno, desc, nec, dg, &
                  ncmpmx, vale, nomcmp, titr, nomnoe, &
                  nomsd, nomsym, ir, numnoe, lmasu, &
                  nbcmp, ncmps, nocmpl)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/ecrtes.h"
#include "asterfort/exisdg.h"
#include "asterfort/irgags.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/lxliis.h"
#include "asterfort/wkvect.h"
    integer :: ifi, nbno, desc(*), nec, dg(*), ncmpmx
    integer :: ir, numnoe(*), nbcmp, ncmps(*), ncmp
    real(kind=8) :: vale(*)
    character(len=*) :: nomcmp(*), nocmpl(*)
    character(len=*) :: titr, nomnoe(*), nomsd, nomsym
    aster_logical :: lmasu
!
!        ECRITURE D'UN CHAM_NO A REPRESENTATION CONSTANTE
!        SUR FICHIER UNIVERSEL, DATASET TYPE 55 A VALEURS REELLES
!      ENTREE:
!         IFI   : UNITE LOGIQUE DU FICHIER UNIVERSEL
!         NBNO  : NOMBRE DE NOEUDS DU LIGREL ( DU MAILLAGE)
!         DESC  :
!         NEC   : NOMBRE D'ENTIERS-CODES
!         DG    : ENTIERS CODES
!         NCMPMX: NOMBRE MAXI DE CMP DE LA GRANDEUR
!         VALE  : VALEURS DU CHAM_NO
!         NOMCMP: NOMS DES CMP
!         TITR  : 1 LIGNE DE TITRE
!         NOMNOE: NOMS DES NOEUDS
!         NUMNOE: NUMEROS DES NOEUDS
!         NOMSD : NOM DU RESULTAT
!         NOMSYM: NOM SYMBOLIQUE
!         IR    : NUMERO D'ORDRE DU CHAMP
!         LMASU : INDIQUE SI MAILLAGE SUPERTAB  .TRUE. MAILLAGE SUPERTAB
!         NBCMP : NOMBRE DE COMPOSANTES DE LA SELECTION A IMPRIMER
!         NCMPS : NUMEROS DES COMPOSANTES DE LA SELECTION A IMPRIMER
!         NOCMPL: NOMS DES COMPOSANTES DE LA SELECTION A IMPRIMER
!
    aster_logical :: lcmp
!     ------------------------------------------------------------------
    character(len=8) :: nocmp, nomgs
    character(len=24) :: nomst
    character(len=80) :: entete(10), titre, texte
    integer :: nbchs
    integer :: iente, iutil
!
!  --- INITIALISATIONS ----
!
!-----------------------------------------------------------------------
    integer :: i, ibcmps, ic, ichs, icmp, icms
    integer :: icompt, icp, ida, idebu, iec, ier, ifin
    integer :: ilig, indats, inno, ino, inom
    integer :: ires, iret, irval, ival, j, jadm, jj
    integer :: jl, jmax, jpos, jtitr, k, l, ll
    integer :: nbcmpt, nbdats, ni
    integer, pointer :: ipcmps(:) => null()
    character(len=8), pointer :: nomchs(:) => null()
    character(len=8), pointer :: nomgds(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    AS_ALLOCATE(vk8=nomgds, size=ncmpmx)
    AS_ALLOCATE(vk8=nomchs, size=ncmpmx)
    AS_ALLOCATE(vi=ipcmps, size=ncmpmx*ncmpmx)
!
    call wkvect('&&IRDRSR.NBCMPS', 'V V I', ncmpmx, ibcmps)
!
    nomst = '&&IRECRI.SOUS_TITRE.TITR'
    call jeveuo(nomst, 'L', jtitr)
    titre = zk80(jtitr)
!
! --- ALLOCATION DES TABLEAUX DE TRAVAIL ---
!
    call jeexin('&&IRDRSR.VAL', iret)
    if (iret .ne. 0) call jedetr('&&IRDRSR.VAL')
    call wkvect('&&IRDRSR.VAL', 'V V R', ncmpmx, irval)
    call jeexin('&&IRDRSR.NOM', iret)
    if (iret .ne. 0) call jedetr('&&IRDRSR.NOM')
    call wkvect('&&IRDRSR.NOM', 'V V K16', ncmpmx, inom)
!
    ncmp = -desc(2)
    do iec = 1, nec
        dg(iec) = desc(3+iec-1)
    end do
    icompt = 0
    do icmp = 1, ncmpmx
        if (exisdg(dg, icmp)) then
            icompt = icompt+1
            zk16(inom-1+icompt) = nomcmp(icmp)
        end if
    end do
!
!     ---- RECHERCHE DES GRANDEURS SUPERTAB -----
!
    call irgags(icompt, zk16(inom), nomsym, nbchs, nomchs, &
                zi(ibcmps), nomgds, ipcmps)
!
!
!      ==================
! ---- PARTIE 1 : NBCMP=0
!      ==================
    if (nbcmp .eq. 0) then
!
!     ---- BOUCLE SUR LES DIVERSES GRANDEURS SUPERTAB ----
        do ichs = 1, nbchs
            iente = 1
            lcmp = .false.
            call ecrtes(nomsd, titr, nomgds(ichs), ir, 'NOEU', &
                        zi(ibcmps-1+ichs), 2, entete, lcmp)
            idebu = 1
            entete(4) = ' '
            texte = ' '
            do icp = 1, zi(ibcmps-1+ichs)
                nocmp = nomcmp(ipcmps((ichs-1)*ncmpmx+icp))
                iutil = lxlgut(nocmp)
                ifin = idebu+iutil
                texte(idebu:ifin) = nocmp(1:iutil)//' '
                idebu = ifin+1
            end do
            iutil = lxlgut(texte)
            jmax = lxlgut(titre)
            jmax = min(jmax, (80-iutil-2))
            entete(4) = titre(1:jmax)//' - '//texte(1:iutil)
            do inno = 1, nbno
                ino = numnoe(inno)
                ival = (ino-1)*ncmp
!
                do ic = 1, zi(ibcmps-1+ichs)
                    zr(irval-1+ic) = 0.0d0
                end do
                do icms = 1, zi(ibcmps-1+ichs)
                    zr(irval-1+icms) = vale(ival+icms)
                end do
                if (iente .eq. 1) then
                    write (ifi, '(A80)') (entete(i), i=1, 10)
                    iente = 0
                end if
                if (lmasu) then
                    call lxliis(nomnoe(inno) (2:8), ino, ier)
                end if
                write (ifi, '(I10,5X,A,A)') ino, '% NOEUD ', nomnoe(inno)
                write (ifi, '(6(1PE13.5E3))') (zr(irval-1+i), i=1, zi( &
                                               ibcmps-1+ichs))
            end do
            if (iente .eq. 0) write (ifi, '(A)') '    -1'
        end do
        call jedetr('&&IRDRSR.VAL')
        call jedetr('&&IRDRSR.NOM')
        call jedetr('&&IRDRSR.NBCMPS')
!
!      =====================
! ---- PARTIE 2 : NBCMP.NE.0
!      =====================
!
    else
!
! --- NOM DE LA GRANDEUR SUPERTAB
        do i = 1, nbchs
            do j = 1, zi(ibcmps+i-1)
                if (ncmps(1) .eq. ipcmps((i-1)*ncmpmx+j)) goto 899
            end do
        end do
899     continue
        nomgs = nomgds(i)
!
! --- NOMBRE DE DATASET
        call wkvect('&&IRDRSR.CMP_DATS', 'V V I', nbcmp, indats)
        nbcmpt = 6
        ilig = nbcmp/6
        ires = nbcmp-ilig*6
        ni = 0
        zi(indats) = ni
        if (ires .eq. 0) then
            nbdats = ilig
            do i = 1, nbdats
                zi(ibcmps+i-1) = 6
                ni = ni+6
                zi(indats+i) = ni
            end do
        else
            nbdats = ilig+1
            do i = 1, nbdats-1
                zi(ibcmps+i-1) = 6
                ni = ni+6
                zi(indats+i) = ni
            end do
            zi(ibcmps+nbdats-1) = ires
            zi(indats+nbdats) = ni+ires
        end if
!
! --- ECRITURE DE L'ENTETE SUPERTAB ----
        lcmp = .true.
        call ecrtes(nomsd, titr, nomgs, ir, 'NOEU', &
                    nbcmpt, 2, entete, lcmp)
!
!
! --- COMPOSANTES ADMISES
        call jedetr('&&IRDESR.CMP')
        call wkvect('&&IRDESR.CMP', 'V V I', ncmpmx, jadm)
        call jedetr('&&IRDESR.POS')
        call wkvect('&&IRDESR.POS', 'V V I', nbcmp, jpos)
        k = 0
        do icmp = 1, ncmpmx
            if (exisdg(dg, icmp)) then
                zi(jadm+k) = icmp
                k = k+1
            end if
        end do
!
! --- BOUCLES SUR LES DATASETS
! ----------------------------
        do ida = 1, nbdats
!
            iente = 1
            ifin = 1
            idebu = 1
            entete(4) = ' '
            texte = ' '
!
            do icp = 1, zi(ibcmps+ida-1)
                nocmp = nocmpl(icp+zi(indats+ida-1))
                iutil = lxlgut(nocmp)
                ifin = idebu+iutil
                texte(idebu:ifin) = nocmp(1:iutil)//' '
                idebu = ifin+1
            end do
!
            iutil = lxlgut(texte)
            jmax = lxlgut(titre)
            jmax = min(jmax, (80-iutil-2))
            entete(4) = titre(1:jmax)//' - '//texte(1:iutil)
!
!
! ---    POSITIONS DES COMPOSANTES SELECTIONNEES PARMI LES
!        COMPOSANTES ADMISES
            l = 0
            do j = 1, zi(ibcmps+ida-1)
                ll = 0
                do jl = 1, ncmp
                    ll = ll+1
                    if (zi(jadm+jl-1) .eq. ncmps(j+zi(indats+ida-1))) goto 780
                end do
780             continue
                zi(jpos+l) = ll
                l = l+1
            end do
!
! ---    BOUCLES SUR LES NOEUDS
!        ---------------------
            do inno = 1, nbno
!
                ino = numnoe(inno)
                jj = (ino-1)*ncmp
!
                do ic = 1, 6
                    zr(irval-1+ic) = 0.0d0
                end do
!
                do icms = 1, zi(ibcmps+ida-1)
                    zr(irval-1+icms) = vale(jj+zi(jpos+icms-1))
                end do
!
                if (iente .eq. 1) then
                    write (ifi, '(A80)') (entete(i), i=1, 10)
                    iente = 0
                end if
                if (lmasu) then
                    call lxliis(nomnoe(inno) (2:8), ino, ier)
                end if
                write (ifi, '(I10,5X,A,A)') ino, '% NOEUD ', nomnoe(inno)
                write (ifi, '(6(1PE13.5E3))') (zr(irval-1+i), i=1, 6)
            end do
            if (iente .eq. 0) write (ifi, '(A)') '    -1'
        end do
!
        call jedetr('&&IRDESR.CMP_DATS')
        call jedetr('&&IRDESR.CMP')
        call jedetr('&&IRDESR.POS')
    end if
!
    AS_DEALLOCATE(vk8=nomgds)
    AS_DEALLOCATE(vk8=nomchs)
    AS_DEALLOCATE(vi=ipcmps)
!
!
    call jedema()
end subroutine
