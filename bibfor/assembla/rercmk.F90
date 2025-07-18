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
subroutine rercmk(nu, mo, ma, nlili, nm, &
                  nl, nbntt)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/indiis.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/renuu1.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: mo, ma
    character(len=14) :: nu
    integer(kind=8) :: nlili, nm, nl, nbntt
! ----------------------------------------------------------------------
!     BUT:  CETTE ROUTINE SERT A RENUMEROTER LES NOEUDS D'UN NUME_DDL
!           SUIVANT L'ALGORITHME DE REVERSE-CUTHILL-MAC-KEE.
!     IN:
!     ---
!       NU : NOM DU NUME_DDL QUE L'ON RENUMEROTE
!       MO : NOM DU MODELE SOUS-JACENT AU NUME_DDL
!       MA : NOM DU MAILLAGE SOUS-JACENT AU NUME_DDL
!       NLILI: NOMBRE DE LIGREL DE L'OBJET .LILI
!       NM   : NOMBRE DE NOEUDS PHYSIQUES DU MAILLAGE
!       NL   : NOMBRE DE NOEUDS DE LAGRANGE DU MAILLAGE
!       NBNTT: NOMBRE DE NOEUDS MAXI (NUME_DDL)
!
!     OUT:
!     ----
!
!      ON REMPLIT LES VECTEURS NU//'.NEWN' ET NU//'.OLDN'
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=8) :: exiele
    character(len=24) :: nomli2
    character(len=19) :: nomlig
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacoin, iaconx, iagrel, ialiel
    integer(kind=8) :: iamail, ianema
    integer(kind=8) :: ico, icol, icumul
    integer(kind=8) :: iel, ifm, igrel, iinew, iino, iio1
    integer(kind=8) :: iio2, ilconx, ili, illiel, ilnema, ima, ino
    integer(kind=8) :: irempl, iret, j, jjno, jno, jrang, k
    integer(kind=8) :: l1, l2, ll1, ll2, longi, longo, n1i
    integer(kind=8) :: n1j, n2i, n2j, nbco, nbcomp, nbel, nbgrel
    integer(kind=8) :: nbi, nbma, nbnm, nbnmre, nbnoma, nbnot, nbntre
    integer(kind=8) :: nbsma, nbssa, newnno, niv
    integer(kind=8), pointer :: lcoi(:) => null()
    integer(kind=8), pointer :: vnbco(:) => null()
    integer(kind=8), pointer :: new1(:) => null()
    integer(kind=8), pointer :: old1(:) => null()
    integer(kind=8), pointer :: ordo(:) => null()
    integer(kind=8), pointer :: oldn(:) => null()
    integer(kind=8), pointer :: newn(:) => null()
    integer(kind=8), pointer :: nbno(:) => null()
    integer(kind=8), pointer :: sssa(:) => null()
    integer(kind=8), pointer :: exi1(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
!-----RECUPERATION DU NIVEAU D'IMPRESSION
!
    call infniv(ifm, niv)
!----------------------------------------------------------------------
    nbnoma = nm+nl
!
!     -- ALLOCATION DES OBJETS .NEW1 ET .OLD1  (PROVISOIRES) :
!        CES OBJETS REPRESENTENT LA RENUMEROTATION DE TOUS LES NOEUDS
!     ---------------------------------------------------------------
    AS_ALLOCATE(vi=new1, size=nbntt)
    AS_ALLOCATE(vi=old1, size=nbntt)
!
!     ---------------------------------------------------------------
!        ON CALCULE LA DIMENSION DE LA TABLE DE CONNECTIVITE INVERSE:
!             (EN FAIT, ON SUR-DIMENSIONNE)
!     ---------------------------------------------------------------
!
!     -- ORDO EST UNE TABLE DE TRAVAIL QUI DOIT POUVOIR CONTENIR
!        UNE LIGNE DE CONNECTIVITE INVERSE.
    AS_ALLOCATE(vi=ordo, size=nbntt)
!
    AS_ALLOCATE(vi=lcoi, size=nbntt)
!     -- .LCOI(INO) CONTIENDRA L'ADRESSE DANS .COIN DE LA LISTE
!        DES NOEUDS CONNECTES A INO (C'EST LE VECTEUR CUMULE DE .EXI1)
!        C'EST EN QUELQUE SORTE LE POINTEUR DE LONGUEUR CUMULEE DE .COIN
!
    AS_ALLOCATE(vi=vnbco, size=nbntt)
!     -- .NBCO(INO) CONTIENDRA AU FUR ET A MESURE DE LA CONSTRUCTION
!         DE LA TABLE DE CONNECTIVITE INVERSE, LE NOMBRE REEL DE NOEUDS
!         CONNECTES A INO.
!
!
!
!     -----------------------------------------------------------------
!        RECUPERATION DE .EXI1
!        ALLOCATION DE LA TABLE DE CONNECTIVITE INVERSE: .COIN
!        REMPLISSAGE DU "POINTEUR DE LONGUEUR" .LCOI
!     -----------------------------------------------------------------
!
    call jeveuo(nu//'.EXI1', 'L', vi=exi1)
!
    icumul = 0
!     -- NBNTRE EST LE NOMBRE TOTAL DE NOEUDS A RENUMEROTER
    nbntre = 0
!
    lcoi(1) = 1
    do 5, ino = 1, nbntt-1
        icumul = icumul+exi1(ino+1)
        lcoi(ino+1) = lcoi(ino)+exi1(ino+1)
        if (exi1(ino+1) .gt. 0) nbntre = nbntre+1
5   end do
!
    icumul = icumul+exi1(nbntt+1)
    if (exi1(nbntt+1) .gt. 0) nbntre = nbntre+1
!
    call wkvect('&&RERCMK.COIN', 'V V I', icumul, iacoin)
!
!
!     -----------------------
!     --REMPLISSAGE DE .COIN:
!     -----------------------
!
    call dismoi('NB_MA_MAILLA', mo, 'MODELE', repi=nbma)
    if (nbma .gt. 0) then
        call jeveuo(ma//'.CONNEX', 'L', iaconx)
        call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', ilconx)
    end if
!
!
!     -- 1ERE ETAPE : (SUPER)MAILLES DU MAILLAGE:
!     -------------------------------------------
    call dismoi('NB_SS_ACTI', mo, 'MODELE', repi=nbssa)
    call dismoi('NB_SM_MAILLA', mo, 'MODELE', repi=nbsma)
    if (nbssa .gt. 0) then
        call jeveuo(mo//'.MODELE    .SSSA', 'L', vi=sssa)
    else
        goto 12
    end if
!
    do ima = 1, nbsma
        if (sssa(ima) .eq. 1) then
            call jeveuo(jexnum(ma//'.SUPMAIL', ima), 'L', iamail)
            call jelira(jexnum(ma//'.SUPMAIL', ima), 'LONMAX', nbnm)
            do i = 1, nbnm
                ino = zi(iamail-1+i)
                iino = ino
                if (ino .le. 0) then
                    call utmess('F', 'ASSEMBLA_36')
                end if
                do j = i+1, nbnm
                    jno = zi(iamail-1+j)
                    jjno = jno
                    jrang = indiis(zi(iacoin+lcoi(iino)-1) &
                                   , jjno, 1, vnbco(iino))
!
                    if (jrang .eq. 0) then
                        irempl = vnbco(iino)+1
                        vnbco(iino) = irempl
                        zi(iacoin+lcoi(iino)-1+irempl-1) = &
                            jjno
!
                        irempl = vnbco(jjno)+1
                        vnbco(jjno) = irempl
                        zi(iacoin+lcoi(jjno)-1+irempl-1) = &
                            iino
                    end if
                end do
            end do
        end if
    end do
!
12  continue
!
!
!     -- 2EME ETAPE : MAILLES TARDIVES (OU NON) DES LIGRELS
!                     (MODELE + LISTE DE CHARGES)
!     -----------------------------------------------------
!
    nbnot = 0
    do 30, ili = 2, nlili
        call jenuno(jexnum(nu//'.NUME.LILI', ili), nomli2)
        nomlig = nomli2(1:19)
        call dismoi('EXI_ELEM', nomlig, 'LIGREL', repk=exiele)
        if (exiele(1:3) .eq. 'NON') goto 30
!
        call jeveuo(nomlig//'.LIEL', 'L', ialiel)
        call jeveuo(jexatr(nomlig//'.LIEL', 'LONCUM'), 'L', illiel)
        call jelira(nomlig//'.LIEL', 'NMAXOC', nbgrel)
!
        call jeexin(nomlig//'.NEMA', iret)
        if (iret .gt. 0) then
            call jeveuo(nomlig//'.NEMA', 'L', ianema)
            call jeveuo(jexatr(nomlig//'.NEMA', 'LONCUM'), 'L', ilnema)
        end if
!
        do igrel = 1, nbgrel
            nbel = zi(illiel-1+igrel+1)-zi(illiel-1+igrel)-1
            iagrel = ialiel+zi(illiel-1+igrel)-1
            do iel = 1, nbel
                ima = zi(iagrel-1+iel)
                if (ima .gt. 0) then
                    nbnm = zi(ilconx-1+ima+1)-zi(ilconx-1+ima)
                    iamail = iaconx+zi(ilconx-1+ima)-1
                else
                    nbnm = zi(ilnema-1-ima+1)-zi(ilnema-1-ima)-1
                    iamail = ianema+zi(ilnema-1-ima)-1
                end if
!
                do i = 1, nbnm
                    ino = zi(iamail-1+i)
                    iino = ino
                    if (ino .lt. 0) iino = nbnoma+nbnot-ino
!
                    do j = i+1, nbnm
                        jno = zi(iamail-1+j)
                        jjno = jno
                        if (jno .lt. 0) jjno = nbnoma+nbnot-jno
!
                        jrang = indiis(zi(iacoin+lcoi(iino)-1) &
                                       , jjno, 1, vnbco(iino))
!
                        if (jrang .eq. 0) then
                            irempl = vnbco(iino)+1
                            vnbco(iino) = irempl
                            zi(iacoin+lcoi(iino)-1+irempl-1) = &
                                jjno
!
                            irempl = vnbco(jjno)+1
                            vnbco(jjno) = irempl
                            zi(iacoin+lcoi(jjno)-1+irempl-1) = &
                                iino
                        end if
                    end do
                end do
            end do
        end do
!
        call jeveuo(nomlig//'.NBNO', 'L', vi=nbno)
        nbnot = nbnot+nbno(1)
30  end do
!
!
!
!     --CALCUL DES OBJETS .NEW1 ET .OLD1 :
!     ------------------------------------
!
    iinew = 0
!
!     -- NBCOMP COMPTE LE NOMBRE DE COMPOSANTES CONNEXES DU MODELE
    nbcomp = 0
50  continue
    nbcomp = nbcomp+1
!
!     --ON INITIALISE L'ALGORITHME PAR LE NOEUD I QUI A LA CONNECTIVITE
!     -- LA PLUS FAIBLE (PARMI CEUX RESTANT A RENUMEROTER):
!     "I= MIN(NBCO)"
!     ----------------------------------------------------------------
    i = 0
    do k = 1, nbntt
        if (exi1(1+k) .eq. 0) goto 51
        if (new1(k) .ne. 0) goto 51
        if (i .eq. 0) then
            i = k
        else
            if (vnbco(k) .lt. vnbco(i)) i = k
        end if
51      continue
    end do
    ASSERT(i .ne. 0)
!
    iinew = iinew+1
    new1(i) = iinew
    old1(iinew) = i
!     -- SI ON A RENUMEROTE TOUS LES NOEUDS ATTENDUS, ON SORT :
    if (iinew .eq. nbntre) goto 200
    ico = iinew
!
100 continue
    longi = vnbco(i)
    call renuu1(zi(iacoin-1+lcoi(i)), longi, ordo, longo, vnbco, &
                new1)
    do j = 1, longo
        iinew = iinew+1
        new1(ordo(j)) = iinew
        old1(iinew) = ordo(j)
!        -- SI ON A RENUMEROTE TOUS LES NOEUDS ATTENDUS, ON SORT :
        if (iinew .eq. nbntre) goto 200
    end do
    ico = ico+1
    i = old1(ico)
    if (i .eq. 0) then
        goto 50
    else
        goto 100
    end if
!
200 continue
!
!     -- ON COMPACTE .OLD1 DANS .NEWN ET .OLDN
!     POUR NE CONSERVER QUE LES NOEUDS PHYSIQUES :
!     --------------------------------------------
    call jeveuo(nu//'.OLDN', 'E', vi=oldn)
    call jeveuo(nu//'.NEWN', 'E', vi=newn)
!
    icol = 0
    do i = 1, nbntt
        iio1 = old1(i)
        if (iio1 .eq. 0) goto 3
        if (iio1 .gt. nm) then
            icol = icol+1
        else
            iio2 = i-icol
            if ((iio1 .lt. 1) .or. (iio1 .gt. nm)) then
                call utmess('F', 'ASSEMBLA_38')
            end if
            if ((iio2 .lt. 1) .or. (iio2 .gt. nm)) then
                call utmess('F', 'ASSEMBLA_38')
            end if
            newn(iio1) = i-icol
        end if
    end do
3   continue
!     -- NBNMRE EST LE NOMBRE DE NOEUDS PHYSIQUES A RENUMEROTER
    nbnmre = iio2
!
!     -- ON FINIT EN "REVERSANT" LE TOUT :
!     ------------------------------------
    do i = 1, nm
        if (newn(i) .eq. 0) goto 300
        newnno = nbnmre+1-newn(i)
        newn(i) = newnno
        oldn(newnno) = i
300     continue
    end do
!
!
!     -- ON ECRIT LES LARGEURS DE BANDE MOYENNES AVANT ET APRES:
!     ----------------------------------------------------------
    if (niv .ge. 1) then
        write (ifm, *) '--- RENUMEROTATION DES NOEUDS DU MODELE (RCMK) :'
        write (ifm, *) '   --- NOMBRE DE COMPOSANTES CONNEXES DU MODELE :'&
     &  , nbcomp
    end if
!
    nbi = 0
    ll1 = 0
    ll2 = 0
    do i = 1, nm
        if (exi1(1+i) .eq. 0) goto 600
        nbi = nbi+1
        nbco = vnbco(i)
        l1 = 1
        l2 = 1
        do j = 1, nbco
            n1i = i
            n1j = zi(iacoin-2+lcoi(i)+j)
            if (n1j .gt. nm) goto 601
            l1 = max(l1, (n1i-n1j)+1)
!
            n2i = newn(n1i)
            n2j = newn(n1j)
            l2 = max(l2, (n2i-n2j)+1)
601         continue
        end do
        ll1 = ll1+l1
        ll2 = ll2+l2
600     continue
    end do
    if (niv .ge. 1) then
        write (ifm, *) '   --- HAUTEUR DE COLONNE MOYENNE (EN NOEUDS)'
        write (ifm, *) '        (EN NE TENANT COMPTE QUE DES NOEUDS '&
     &                       //'PHYSIQUES)'
        write (ifm, fmt='(A30,1PD10.3)') '        AVANT RENUMEROTATION:',&
     &             dble(ll1)/nbi
        write (ifm, fmt='(A30,1PD10.3)') '        APRES RENUMEROTATION:',&
     &             dble(ll2)/nbi
!
        if (ll1 .le. ll2) then
            write (ifm, *) '   --- LA NOUVELLE NUMEROTATION OBTENUE PAR '&
     &   //'L ALGORITHME "RCMK" NE SEMBLE PAS'
            write (ifm, *) '       MEILLEURE QUE L ORIGINALE. ELLE L''EST'&
     &   //' PEUT ETRE QUAND MEME DU FAIT DE LA '
            write (ifm, *) '       PRISE EN COMPTE DES RELATIONS LINEAIRES'&
     &  //' ENTRE NOEUDS.'
        end if
    end if
!
    AS_DEALLOCATE(vi=new1)
    AS_DEALLOCATE(vi=old1)
    AS_DEALLOCATE(vi=ordo)
    AS_DEALLOCATE(vi=lcoi)
    AS_DEALLOCATE(vi=vnbco)
    call jedetr('&&RERCMK.COIN')
!
    call jedema()
end subroutine
