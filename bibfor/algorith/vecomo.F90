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

subroutine vecomo(modgen, sst1, sst2, intf1, intf2, &
                  nliais, option)
    implicit none
!  C. VARE     DATE 18/09/94
!-----------------------------------------------------------------------
!  BUT:  < VERIFIER LA COHERENCE DU MODELE GENERALISE >
!
!  ON VERIFIE QUE LA LIAISON DEFINIE DANS DEFI_MODELE_GENE EST
!  COMPATIBLE AVEC LES ORIENTATIONS ET LES TRANSLATIONS AFFECTEES AUX
!  SOUS-STRUCTURES. LES NOEUDS DES DEUX INTERFACES QUI FORMENT LA
!  LIAISON DOIVENT ETRE CONFONDUS 2 A 2
!
!-----------------------------------------------------------------------
!
! MODGEN  /I/ : NOM K8 DU MODELE GENERALISE
! SST1    /I/ : NOM K8 DE LA PREMIERE SOUS-STRUCTURE DE LA LIAISON
! SST2    /I/ : NOM K8 DE LA SECONDE SOUS-STRUCTURE DE LA LIAISON
! INTF1   /I/ : NOM K8 DE L'INTERFACE DE SST1
! INTF2   /I/ : NOM K8 DE L'INTERFACE DE SST2
! NLIAIS  /I/ : NUMERO DE L'INTERFACE COURANTE
!         /O/ : CODE RETOUR : 1 SI INTERFACE COMPATIBLE, 0 SINON
! OPTION  /I/ : INDIQUE SI ON EST EN CLASSIQUE OU EN REDUIT
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/intet0.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mgutdm.h"
#include "asterfort/pmppr.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
!
!
!   PARAMETRE REPRESENTANT LE NOMBRE MAX DE COMPOSANTES DE LA GRANDEUR
!   SOUS-JACENTE TRAITEE
!
    aster_logical :: saut, ordre
    integer(kind=8) :: nbno, icrit, nliais, ldlid, llint3, llint4, iret, nbec, nbcmpm
    integer(kind=8) :: ival, nusst1, nusst2, llrot1, llrot2, lltra1, lltra2, i, j, k, l
    integer(kind=8) :: ibid, nbno1, nbno2, llint1, llint2
    integer(kind=8) :: llistb, inu1, nuno1, inu2, nuno2, jnode, ip, inu
    integer(kind=8) :: nuno, ldac2
    parameter(nbcmpm=10)
    character(len=4) :: nliai
    character(len=8) :: modgen, lint1, lint2, criter, temp
    character(len=8) :: sst1, sst2, intf1, intf2, mail1, mail2, nomnoi
    character(len=8) :: nomnoj, nomnop, option
    character(len=24) :: repnom, int1, int2, famli, ordol, ordod
    character(len=24) :: valk(5)
    real(kind=8) :: x1(3), x2(3), xr1(3), xr2(3), rot1(3), rot2(3), dxr
    real(kind=8) :: mat1(nbcmpm, nbcmpm), mat2(nbcmpm, nbcmpm), tra1(3)
    real(kind=8) :: mat3(nbcmpm, nbcmpm), zero, dxrm, lcaram, tra2(3)
    real(kind=8) :: mattmp(nbcmpm, nbcmpm), difmax, lcara1, lcara2
    real(kind=8) :: matro1(nbcmpm, nbcmpm), matro2(nbcmpm, nbcmpm), seuil, dxrij
    integer(kind=8), pointer :: lista(:) => null()
    real(kind=8), pointer :: nllcoo1(:) => null()
    real(kind=8), pointer :: nllcoo2(:) => null()
    integer(kind=8), pointer :: nldesc1(:) => null()
    integer(kind=8), pointer :: nldesc2(:) => null()
!
!-----------------------------------------------------------------------
    data zero/0.0d+00/
!-----------------------------------------------------------------------
!
!-----SEUIL DE TOLERANCE ET CRITERE DE PRECISION
!
    call jemarq()
    difmax = 1.d-3
    call getvr8('VERIF', 'PRECISION', iocc=1, scal=seuil, nbret=ival)
    if (ival .ne. 0) difmax = seuil
    icrit = 1
    call getvtx('VERIF', 'CRITERE', iocc=1, scal=criter, nbret=ival)
    if (ival .ne. 0) then
        if (criter .eq. 'ABSOLU') icrit = 2
    end if
!
!-----RECUPERATION DES ROTATIONS ET DES TRANSLATIONS
!
    repnom = modgen//'      .MODG.SSNO'
    call jenonu(jexnom(repnom, sst1), nusst1)
    call jenonu(jexnom(repnom, sst2), nusst2)
    call jeveuo(jexnum(modgen//'      .MODG.SSOR', nusst1), 'L', llrot1)
    call jeveuo(jexnum(modgen//'      .MODG.SSOR', nusst2), 'L', llrot2)
    do i = 1, 3
        rot1(i) = zr(llrot1+i-1)
        rot2(i) = zr(llrot2+i-1)
    end do
    call jeveuo(jexnum(modgen//'      .MODG.SSTR', nusst1), 'L', lltra1)
    call jeveuo(jexnum(modgen//'      .MODG.SSTR', nusst2), 'L', lltra2)
    do i = 1, 3
        tra1(i) = zr(lltra1+i-1)
        tra2(i) = zr(lltra2+i-1)
    end do
!
!-----CALCUL DES MATRICES DE ROTATION
!
    call intet0(rot1(1), mat1, 3)
    call intet0(rot1(2), mat2, 2)
    call intet0(rot1(3), mat3, 1)
    call r8inir(nbcmpm*nbcmpm, zero, mattmp, 1)
    call pmppr(mat1, nbcmpm, nbcmpm, 1, mat2, &
               nbcmpm, nbcmpm, 1, mattmp, nbcmpm, &
               nbcmpm)
    call r8inir(nbcmpm*nbcmpm, zero, matro1, 1)
    call pmppr(mattmp, nbcmpm, nbcmpm, 1, mat3, &
               nbcmpm, nbcmpm, 1, matro1, nbcmpm, &
               nbcmpm)
!
    call intet0(rot2(1), mat1, 3)
    call intet0(rot2(2), mat2, 2)
    call intet0(rot2(3), mat3, 1)
    call r8inir(nbcmpm*nbcmpm, zero, mattmp, 1)
    call pmppr(mat1, nbcmpm, nbcmpm, 1, mat2, &
               nbcmpm, nbcmpm, 1, mattmp, nbcmpm, &
               nbcmpm)
    call r8inir(nbcmpm*nbcmpm, zero, matro2, 1)
    call pmppr(mattmp, nbcmpm, nbcmpm, 1, mat3, &
               nbcmpm, nbcmpm, 1, matro2, nbcmpm, &
               nbcmpm)
!
!
!-----RECUPERATION MAILLAGE ET INTERFACE AMONT DES SOUS-STRUCTURES
!
    call mgutdm(modgen, sst1, ibid, 'NOM_MAILLAGE', ibid, &
                mail1)
    call mgutdm(modgen, sst1, ibid, 'NOM_LIST_INTERF', ibid, &
                lint1)
!
    call mgutdm(modgen, sst2, ibid, 'NOM_MAILLAGE', ibid, &
                mail2)
    call mgutdm(modgen, sst2, ibid, 'NOM_LIST_INTERF', ibid, &
                lint2)
!
!-----RECUPERATION DU NOMBRE DES NOEUDS DE L'INTERFACE
!
    int1 = lint1//'.IDC_LINO'
    call jenonu(jexnom(int1(1:13)//'NOMS', intf1), ibid)
    call jelira(jexnum(int1(1:17), ibid), 'LONMAX', nbno1)
!
    int2 = lint2//'.IDC_LINO'
    call jenonu(jexnom(int2(1:13)//'NOMS', intf2), ibid)
    call jelira(jexnum(int2(1:17), ibid), 'LONMAX', nbno2)
!
    if (nbno1 .ne. nbno2) then
        valk(1) = sst1
        valk(2) = intf1
        valk(3) = sst2
        valk(4) = intf2
        if (option(1:6) .eq. 'REDUIT') then
            nliais = 0
        else
            call utmess('F', 'ALGORITH16_44', nk=4, valk=valk)
        end if
        goto 999
    end if
    nbno = nbno1
!
!C
!CC---ON VERIFIE LA COINCIDENCE DE CHAQUE COUPLE DE NOEUDS
!C
!
    call jenonu(jexnom(lint1//'.IDC_NOMS', intf1), ibid)
    call jeveuo(jexnum(lint1//'.IDC_LINO', ibid), 'L', llint1)
    call jeveuo(lint1//'.IDC_DEFO', 'L', vi=nldesc1)
    call jeveuo(mail1//'.COORDO    .VALE', 'L', vr=nllcoo1)
!
    call jenonu(jexnom(lint2//'.IDC_NOMS', intf2), ibid)
    call jeveuo(jexnum(lint2//'.IDC_LINO', ibid), 'L', llint2)
    call jeveuo(lint2//'.IDC_DEFO', 'L', vi=nldesc2)
    call jeveuo(mail2//'.COORDO    .VALE', 'L', vr=nllcoo2)
!
!     --- CONSTITUTION DE LISTA ET LISTB :
!         LE IEME NOEUD DE L'INTERFACE DROITE A POUR VIS-A-VIS
!         LE ZI(LISTA-1+I) EME NOEUD DE L'INTERFACE GAUCHE
!         RECIPROQUEMENT LE NOEUD DE POSITION J DE L'INTERFACE GAUCHE
!         EST LE VIS-A-VIS DU NOEUD DE POSITION ZI(LISTB-1+J) DE
!         L'INTERFACE DROITE.
    AS_ALLOCATE(vi=lista, size=nbno)
    call wkvect('&&VECOMO.LISTB', 'V V I', nbno, llistb)
    dxrm = 0.d0
    lcaram = 0.d0
    ordre = .true.
!
    do i = 1, nbno
!     ---RECUPERATION DES COORDONNEES DES NOEUDS DE L'INTERFACE DROITE
!
        inu1 = zi(llint1-1+i)
        nuno1 = nldesc1(inu1)
!
        do k = 1, 3
            x1(k) = nllcoo1(1+(nuno1-1)*3+k-1)
        end do
        do k = 1, 3
            xr1(k) = 0.d0
            do l = 1, 3
                xr1(k) = xr1(k)+matro1(k, l)*x1(l)
            end do
            xr1(k) = xr1(k)+tra1(k)
        end do
!
        dxr = 0.d0
        do j = 1, nbno
!       ---RECUPERATION DES COORDONNEES DES NOEUDS DE L'INTERFACE GAUCHE
!
            inu2 = zi(llint2-1+j)
            nuno2 = nldesc2(inu2)
!
            saut = .false.
            do k = 1, 3
                x2(k) = nllcoo2(1+(nuno2-1)*3+k-1)
            end do
            do k = 1, 3
                xr2(k) = 0.d0
                do l = 1, 3
                    xr2(k) = xr2(k)+matro2(k, l)*x2(l)
                end do
                xr2(k) = xr2(k)+tra2(k)
                if (j .ne. 1 .and. abs(xr2(k)-xr1(k)) .gt. dxr) then
!               --- COMPARAISON COMPOSANTE AVEC DISTANCE --
!                   (SI COMPOSANTE > DISTANCE MIN ALORS
!                     TEST SUR DISTANCE INUTILE ET SAUT=.TRUE.)
                    if (j .eq. i .and. icrit .eq. 1) then
                        saut = .true.
                    else
                        goto 120
                    end if
                end if
            end do
!
!          ---CALCUL DE LA DIFFERENCE DES DISTANCES NOEUD A NOEUD
!
            if (.not. saut) then
                dxrij = 0.d0
                do k = 1, 3
                    dxrij = dxrij+(xr1(k)-xr2(k))**2
                end do
                dxrij = sqrt(dxrij)
                if (j .eq. 1 .or. dxrij .lt. dxr) then
!             --- CRITERE SUR DISTANCE (RECHERCHE DU MINIMUM)
                    dxr = dxrij
                    jnode = j
                end if
            end if
!
!          ---CALCUL D'UNE LONGUEUR CARACTERISTIQUE SI CRITERE RELATIF
!
            if (icrit .eq. 1 .and. j .eq. i) then
                lcara1 = 0.d0
                lcara2 = 0.d0
                do k = 1, 3
                    lcara1 = lcara1+xr1(k)**2
                    lcara2 = lcara2+xr2(k)**2
                end do
                lcara1 = sqrt(lcara1)
                lcara2 = sqrt(lcara2)
                if (lcaram .lt. lcara1) lcaram = lcara1
                if (lcaram .lt. lcara2) lcaram = lcara2
            end if
!
120         continue
        end do
!
!
        if (dxrm .lt. dxr) dxrm = dxr
        lista(i) = jnode
        if (zi(llistb-1+jnode) .ne. 0) then
!        --- CAS OU JNODE EST DEJA UN VIS-A-VIS ---
            ip = zi(llistb-1+jnode)
            inu = zi(llint1-1+i)
            nuno = nldesc1(inu)
            nomnoi = int_to_char8(nuno)
            inu = zi(llint2-1+jnode)
            nuno = nldesc2(inu)
            nomnoj = int_to_char8(nuno)
            inu = zi(llint1-1+ip)
            nuno = nldesc1(inu)
            nomnop = int_to_char8(nuno)
            valk(1) = nomnoj
            valk(2) = nomnop
            valk(3) = nomnoi
!
            call utmess('F', 'ALGORITH16_45', nk=3, valk=valk)
            goto 999
        end if
        zi(llistb-1+jnode) = i
!
!        SI JNODE EST DIFFERENT DE I, C'EST QUE LES NOEUDS D'INTERFACE
!        ONT ETE DONNES DANS UN ORDRE DE NON CORRESPONDANCE
        if (jnode .ne. i) ordre = .false.
!
    end do
!
!
!-----VERIFICATION FINALE
!
    if (icrit .eq. 1) then
        if (lcaram .eq. 0.d0) then
            valk(1) = sst1
            valk(2) = intf1
            valk(3) = sst2
            valk(4) = intf2
            call utmess('F', 'ALGORITH16_46', nk=4, valk=valk)
            goto 999
        end if
        dxrm = dxrm/lcaram
    end if
    if (dxrm .gt. difmax) then
        valk(1) = sst1
        valk(2) = intf1
        valk(3) = sst2
        valk(4) = intf2
        if (option(1:6) .eq. 'REDUIT') then
            nliais = 0
        else
            call utmess('F', 'ALGORITH16_47', nk=4, valk=valk)
        end if
        goto 999
    end if
!
    if (.not. ordre) then
!
!       --- LES NOEUDS NE SONT PAS EN VIS-A-VIS ---
!           ON REGARDE D'ABORD SI LE TRI EST PLAUSIBLE
        do i = 1, nbno
            if (zi(llistb-1+lista(i)) .ne. i) then
                valk(1) = sst1
                valk(2) = intf1
                valk(3) = sst2
                valk(4) = intf2
                if (option(1:6) .eq. 'REDUIT') then
                    nliais = 0
                else
                    call utmess('F', 'ALGORITH16_48', nk=4, valk=valk)
                end if
                goto 999
            end if
        end do
!
!        ON RECUPERE LE DESCRIPTEUR DE LA LIAISON COURANTE
        famli = modgen//'      .MODG.LIDF'
        call jeveuo(jexnum(famli, nliais), 'E', ldlid)
!        ON DIT OUI AU REORDONANCEMENT DES NOEUDS
        zk8(ldlid+4) = 'OUI'
        call codent(nliais, 'D', nliai)
!
        temp = '&&OP0126'
        ordol = temp//'      .LINO.'//nliai
        ordod = temp//'      .LDAC.'//nliai
        call jeexin(ordol, iret)
        if (iret .eq. 0) then
            call jecreo(ordol, 'V V I')
            call jeecra(ordol, 'LONMAX', nbno)
            call dismoi('NB_EC', lint2, 'INTERF_DYNA', repi=nbec)
            call jecreo(ordod, 'V V I')
            call jeecra(ordod, 'LONMAX', nbno*nbec)
        end if
!
        call jeveuo(ordol, 'E', llint3)
!
!    ---  ON ORDONNE LES NOEUDS DE LLINT2 SUIVANT LLISTA
        do i = 1, nbno
!         --- RECOPIE DE LLINT2 DANS LLISTB
            zi(llistb-1+i) = zi(llint2-1+i)
        end do
        do i = 1, nbno
            zi(llint3-1+i) = zi(llistb-1+lista(i))
        end do
!    ---  ON REORDONNE LES CODES DE CONDITIONS AUX LIMITES
!         AFIN D'AVOIR UNE VERIFICATION CORRECTE DANS VERILI
        call jenonu(jexnom(lint2//'.IDC_NOMS', intf2), ibid)
        call jeveuo(jexnum(lint2//'.IDC_DDAC', ibid), 'L', ldac2)
!
        call jeveuo(ordod, 'E', llint4)
        call jedetr('&&VECOMO.LISTB')
        call wkvect('&&VECOMO.LISTB', 'V V I', nbno*nbec, llistb)
        do i = 1, nbno
!         --- RECOPIE DE LDAC2 DANS LLISTB
            zi(llistb+(i-1)*nbec) = zi(ldac2+(i-1)*nbec)
        end do
        do i = 1, nbno
            zi(llint4+(i-1)*nbec) = zi(llistb+(lista(i)-1)*nbec)
        end do
!
    end if
    nliais = 1
!
999 continue
!
!       --- DESTRUCTION OBJETS SUR VOLATILE
    AS_DEALLOCATE(vi=lista)
    call jedetr('&&VECOMO.LISTB')
!
!
    call jedema()
end subroutine
