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

subroutine xddlim(modele, motcle, nomn, ino, valimr, &
                  valimc, valimf, fonree, icompt, lisrel, &
                  ndim, direct, jnoxfv, ch1, ch2, &
                  ch3, cnxinv, mesh, hea_no)
    implicit none
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8pi.h"
#include "asterfort/assert.h"
#include "asterfort/afrela.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xddlimf.h"
!
    integer(kind=8) :: ino, icompt, ndim, jnoxfv
    real(kind=8) :: valimr, direct(3)
    character(len=4) :: fonree
    character(len=8) :: modele, nomn, valimf, motcle
    character(len=8), intent(in) :: mesh
    character(len=19) :: lisrel, ch1, ch2, ch3, cnxinv, hea_no
! person_in_charge: samuel.geniaut at edf.fr
!
!      TRAITEMENT DE DDL_IMPO SUR UN NOEUD X-FEM
!             (POUR MOTCLE = DX, DY ,DZ)
!      TRAITEMENT DE DDL_IMPO SUR UN NOEUD HM-XFEM
!             (POUR MOTCLE = DX, DY, DZ ET/OU PRE1)
!      TRAITEMENT DE FACE_IMPO SUR UN NOEUD X-FEM
!             (POUR DNOR OU DTAN : MOTCLE = DEPL )
!
! IN  MODELE : NOM DE L'OBJET MODELE ASSOCIE AU LIGREL DE CHARGE
! IN  MOTCLE : NOM DE LA COMPOSANTE DU DEPLACEMENT/PRESSION A IMPOSER
! IN  NOMN   : NOM DU NOEUD INO OU EST EFFECTUE LE BLOCAGE
! IN  INO    : NUMERO DU NOEUD OU EST EFFECTUE LE BLOCAGE
! IN  VALIMR : VALEUR DE BLOCAGE SUR CE DDL (FONREE = 'REEL')
! IN  VALIMC : VALEUR DE BLOCAGE SUR CE DDL (FONREE = 'COMP')
! IN  VALIMF : VALEUR DE BLOCAGE SUR CE DDL (FONREE = 'FONC')
! IN  FONREE : AFFE_CHAR_XXXX OU AFFE_CHAR_XXXX_F
! IN  NDIM
! IN  MESH   : NOM DU MAILLAGE
!
! IN/OUT
!     ICOMPT : "COMPTEUR" DES DDLS AFFECTES REELLEMENT
!     LISREL : LISTE DE RELATIONS AFFECTEE PAR LA ROUTINE
!
!
!
!
    integer(kind=8) :: nbxcmp, nfimax
    parameter(nbxcmp=60, nfimax=10)
    integer(kind=8) :: ier, stano(4), jstnol, jstnod, nrel, fisco(2*nfimax)
    integer(kind=8) ::  jlsnl, jlsnd, jlstl, jlstd, jfiscl, jfiscd, fisc(2*nfimax)
    integer(kind=8) ::  jfisnl, jfisnd, nfh, ifh, nfisc
    integer(kind=8) ::  i, j, nterm, irel, dimens(nbxcmp), ifiss, nfiss, ifisc
    integer(kind=8) ::  nbno, nbmano, adrma, ima, numa, nbnoma, nuno, nuno2
    integer(kind=8) ::  jconx2, iad, fisno(4)
    integer(kind=8) ::  jheavnl, jheavnd, ncompn, heavn(5), hea_se
    real(kind=8) :: r, theta(2), he(2, 4), t, coef(nbxcmp), sign
    real(kind=8) :: lsn(4), lst(4), minlsn, maxlsn, lsn2, ljonc(nfimax)
    character(len=8) :: ddl(nbxcmp), noeud(nbxcmp), noma
    character(len=1) :: axes(3)
    character(len=19) :: ch4, ch5
    complex(kind=8) :: cbid, valimc
    aster_logical :: class
    integer(kind=8), pointer :: nunotmp(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    integer(kind=8), pointer :: connex(:) => null(), ihea_no(:) => null()
    integer(kind=8), pointer :: fisnv(:) => null()
    integer(kind=8), pointer :: fiscv(:) => null()
    real(kind=8), pointer :: lsnv(:) => null()
    real(kind=8), pointer :: lstv(:) => null()
    integer(kind=8), pointer :: stnov(:) => null()
    cbid = dcmplx(0.d0, 0.d0)
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- RECUP DU NUMERO LOCAL NUMO DU NOEUD INO DANS LA MAILLE X-FEM NUMA
    numa = zi(jnoxfv-1+2*(ino-1)+1)
    nuno = zi(jnoxfv-1+2*(ino-1)+2)
    axes(1) = 'X'
    axes(2) = 'Y'
    axes(3) = 'Z'
!
!       REMARQUE : FAIRE DES CALL JEVEUO AU COEUR DE LA BOUCLE SUR
!      LES NOEUDS ET SUR LES DDLS BLOQUES N'EST PAS OPTIMAL DU POINT
!      DE VUE DES PERFORMANCES, MAIS A PRIORI, CA NE DEVRAIT PAS ETRE
!      POUR BEAUCOUP DE NOEUDS
    call jeveuo(ch1//'.CESV', 'L', vi=stnov)
    call jeveuo(ch1//'.CESL', 'L', jstnol)
    call jeveuo(ch1//'.CESD', 'L', jstnod)
    call jeveuo(ch2//'.CESV', 'L', vr=lsnv)
    call jeveuo(ch2//'.CESL', 'L', jlsnl)
    call jeveuo(ch2//'.CESD', 'L', jlsnd)
    call jeveuo(ch3//'.CESV', 'L', vr=lstv)
    call jeveuo(ch3//'.CESL', 'L', jlstl)
    call jeveuo(ch3//'.CESD', 'L', jlstd)
!
! --- NOMBRE DE FISSURES VUES PAR LA MAILLE
    nfiss = zi(jstnod-1+5+4*(numa-1)+2)
!
    do i = 1, 2*nfimax
        fisco(i) = 0
        fisc(i) = 0
    end do
    if (nfiss .gt. 1) then
        ch4 = '&&XDDLIM.CHS4'
        call celces(modele//'.FISSNO', 'V', ch4)
        call jeveuo(ch4//'.CESV', 'L', vi=fisnv)
        call jeveuo(ch4//'.CESL', 'L', jfisnl)
        call jeveuo(ch4//'.CESD', 'L', jfisnd)
! --- NOMBRE DE DDLS HEAVISIDES DANS LA MAILLE
        nfh = zi(jfisnd-1+5+4*(numa-1)+2)
        do i = 1, nfh
            call cesexi('S', jfisnd, jfisnl, numa, nuno, &
                        i, 1, iad)
            fisno(i) = fisnv(iad)
        end do
        ch5 = '&&XDDLIM.CHS5'
        call celces(modele//'.FISSCO', 'V', ch5)
        call jeveuo(ch5//'.CESV', 'L', vi=fiscv)
        call jeveuo(ch5//'.CESL', 'L', jfiscl)
        call jeveuo(ch5//'.CESD', 'L', jfiscd)
        do i = 1, nfiss
            call cesexi('S', jfiscd, jfiscl, numa, 1, &
                        i, 1, iad)
            fisco(2*i-1) = fiscv(iad)
            call cesexi('S', jfiscd, jfiscl, numa, 1, &
                        i, 2, iad)
            fisco(2*i) = fiscv(iad)
        end do
    end if
    do ifiss = 1, nfiss
        call cesexi('S', jstnod, jstnol, numa, nuno, &
                    ifiss, 1, iad)
        stano(ifiss) = stnov(iad)
        call cesexi('S', jlsnd, jlsnl, numa, nuno, &
                    ifiss, 1, iad)
        lsn(ifiss) = lsnv(iad)
        call cesexi('S', jlstd, jlstl, numa, nuno, &
                    ifiss, 1, iad)
        lst(ifiss) = lstv(iad)
    end do
!
! --- RECUPERATION DE LA DEFINITION DES DDLS HEAVISIDES
    hea_se = -99
    heavn(1:5) = 0
    if (stano(1) .eq. 1 .or. stano(1) .eq. 3) then
        call jeveuo(hea_no//'.CESV', 'L', vi=ihea_no)
        call jeveuo(hea_no//'.CESL', 'L', jheavnl)
        call jeveuo(hea_no//'.CESD', 'L', jheavnd)
        ncompn = zi(jheavnd-1+5+4*(numa-1)+3)
        ASSERT(ncompn .eq. 5)
        do i = 1, ncompn
            call cesexi('S', jheavnd, jheavnl, numa, nuno, &
                        1, i, iad)
            heavn(i) = ihea_no(iad)
        end do
    end if
!
! --- IDENTIFICATIOND DES CAS A TRAITER :
! --- SI LA RELATION CINEMATIQUE EST IMPOSEE PAR DES VALEURS REELLES ET
! --- SI NOEUD SUR LES LEVRES ET CONNECTÉ À DES NOEUDS (APPARTENANT AU
! --- GROUPE AFFECTÉ) DE PART ET D'AUTRE DE LA LEVRE : 2 RELATIONS
! --- SINON IL FAUT IMPOSER QUE D'UN SEUL COTÉ        : 1 RELATION
    if (nfiss .eq. 1) then
        if (lsn(1) .eq. 0.d0 .and. lst(1) .lt. 0.d0) then
            minlsn = r8maem()
            maxlsn = -1*r8maem()
! ---     RECUPERATION DE LA LISTE DES NOEUDS AFFECTÉS PAR LA CONDITION
            call jeexin('&&CADDLI.LIST_NODE', ier)
            if (ier .ne. 0) then
                call jeveuo('&&CADDLI.LIST_NODE', 'L', vi=nunotmp)
                call jelira('&&CADDLI.LIST_NODE', 'LONMAX', nbno)
            end if
            call jeexin('&&CAFACI.LIST_NODE', ier)
            if (ier .ne. 0) then
                call jeveuo('&&CAFACI.LIST_NODE', 'L', vi=nunotmp)
                call jelira('&&CAFACI.LIST_NODE', 'LONMAX', nbno)
            end if
            call jeexin('&&CAAREI.LIST_NODE', ier)
            if (ier .ne. 0) then
                call jeveuo('&&CAAREI.LIST_NODE', 'L', vi=nunotmp)
                call jelira('&&CAAREI.LIST_NODE', 'LONMAX', nbno)
            end if
! ---     RECUPERATION DU NOM DU MAILLAGE :
            call jeveuo(modele//'.MODELE    .LGRF', 'L', vk8=lgrf)
            noma = lgrf(1)
! ---     RECUPERATION DES MAILLES CONTENANT LE NOEUD
            call jeveuo(noma//'.CONNEX', 'L', vi=connex)
            call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
            call jelira(jexnum(cnxinv, ino), 'LONMAX', nbmano)
            call jeveuo(jexnum(cnxinv, ino), 'L', adrma)
! ---     BOUCLE SUR LES MAILLES CONTENANT LE NOEUD
            do ima = 1, nbmano
                numa = zi(adrma-1+ima)
                nbnoma = zi(jconx2+numa)-zi(jconx2+numa-1)
! ---       BOUCLE SUR LES NOEUDS DE LA MAILLE
! ---       ATTENTION ON NE PREND EN COMPTE QUE LES MAILLES LINEAIRES !
                do i = 1, nbnoma
                    nuno = connex(zi(jconx2+numa-1)+i-1)
! ---         ON REGARDE SI LE NOEUD APPARTIENT AU GRP DE NOEUD AFFECTÉ
                    do j = 1, nbno
                        nuno2 = nunotmp(j)
                        if (nuno2 .eq. nuno) then
                            call cesexi('C', jlsnd, jlsnl, numa, i, &
                                        1, 1, iad)
                            if (iad .le. 0) goto 110
                            lsn2 = lsnv(iad)
!                  LSN2 = ZR(JLSN-1+NUNO)
                            if (lsn2 .lt. minlsn) minlsn = lsn2
                            if (lsn2 .gt. maxlsn) maxlsn = lsn2
                            goto 110
                        end if
                    end do
110                 continue
                end do
            end do
!
            if ((minlsn .eq. 0.d0) .and. (maxlsn .gt. 0.d0)) then
! ---       ON AFFECTE LA RELATION UNIQUEMENT SUR LA PARTIE MAITRE
                nrel = 1
                theta(1) = r8pi()
                he(1, 1) = 1.d0
            else if ((minlsn .lt. 0.d0) .and. (maxlsn .eq. 0.d0)) then
! ---       ON AFFECTE LA RELATION UNIQUEMENT SUR LA PARTIE ESCLAVE
                nrel = 1
                theta(1) = r8pi()
                he(1, 1) = -1.d0
            elseif (((minlsn .lt. 0.d0) .and. (maxlsn .gt. 0.d0)) .or. ( &
                    nbno .eq. 0)) then
! ---       ON AFFECTE LA RELATION SUR LES DEUX PARTIES
                nrel = 2
                theta(1) = r8pi()
                theta(2) = -r8pi()
                he(1, 1) = 1.d0
                he(2, 1) = -1.d0
            else
! ---       SI NOEUD ISOLE, ON AFFECTE RIEN POUR L'INSTANT
                goto 888
            end if
        else
            nrel = 1
            he(1, 1) = sign(1.d0, lsn(1))
            theta(1) = he(1, 1)*abs(atan2(lsn(1), lst(1)))
        end if
    else if (nfiss .gt. 1) then
        do ifiss = 1, nfiss
            do i = 1, ifiss
                fisc(i) = 0
            end do
            ifisc = ifiss
            nfisc = 0
80          continue
            if (fisco(2*ifisc-1) .gt. 0) then
                nfisc = nfisc+1
                fisc(2*(nfisc-1)+2) = fisco(2*ifisc)
                ifisc = fisco(2*ifisc-1)
                fisc(2*(nfisc-1)+1) = ifisc
                goto 80
            end if
            do i = 1, nfisc
                ljonc(i) = lsn(fisc(2*i-1))
            end do
!   MISE A ZERO POUR LA FONCTION JONCTION AU NIVEAU DU BRANCHEMENT
            do i = 1, nfisc
                if (fisc(2*i)*ljonc(i) .gt. 0.d0) then
                    lsn(ifiss) = 0.d0
                end if
            end do
        end do
!
        nrel = 1
        do ifh = 1, nfh
! --- ON NE PREND PAS ENCORE EN COMPTE LE CAS OU ON PASSE PAR UN NOEUD POUR
! --- LES ELEMENTS MULTI-HEAVISIDE
            if (lsn(fisno(ifh)) .eq. 0 .and. nfisc .eq. 0) goto 888
        end do
    end if
!
    do i = 1, nbxcmp
        dimens(i) = 0
        noeud(i) = nomn
    end do
!
    if (nrel .eq. 2 .and. fonree .eq. 'REEL') then
        call utmess('A', 'XFEM_22', sk=nomn)
    end if
!
    if (fonree .eq. 'FONC' .and. nfiss .eq. 1 .and. (stano(1) .eq. 1 &
                                                     .or. stano(1) .eq. 3)) then
! --- SI LA RELATION CINEMATIQUE EST IMPOSEE PAR UNE FONCTION DE L'ESPACE
        call xddlimf(modele, ino, cnxinv, jnoxfv, motcle, &
                     ch2, ndim, lsn, lst, valimr, valimf, valimc, &
                     fonree, lisrel, nomn, direct, class, mesh, &
                     hea_no)
    end if
!     IMPOSITION DES CONDITIONS CINEMATIQUE "TOTALES" (DDL_CLASS +/- DDL_ENR)
    if ((fonree .eq. 'REEL') .or. (nfiss .gt. 1) .or. class .or. (stano(1) .ne. 1 &
                                                                  .and. stano(1) .ne. 3)) then
        do i = 1, nbxcmp
            dimens(i) = 0
            noeud(i) = nomn
        end do
! --- BOUCLE SUR LES RELATIONS
        do irel = 1, nrel
!
! --- CALCUL DU SOUS DOMAINE CORRESPONDANT A CHAQUE RELATION LINEAIRE
            if (nfiss .eq. 1) hea_se = xcalc_code(1, he_real=[he(irel, 1)])
!
!       CALCUL DES COORDONNÉES POLAIRES DU NOEUD (R,T)
            r = sqrt(lsn(1)**2+lst(1)**2)
            t = theta(irel)
!
!       CAS FACE_IMPO DNOR OU DTAN
            if (motcle(1:8) .eq. 'DEPL    ') then
!
                i = 0
                do j = 1, ndim
!
!           COEFFICIENTS ET DDLS DE LA RELATION
                    i = i+1
                    ddl(i) = 'D'//axes(j)
                    coef(i) = direct(j)
!
                    if (nfiss .eq. 1) then
                        if (stano(1) .eq. 1 .or. stano(1) .eq. 3) then
                            i = i+1
                            ddl(i) = 'H1'//axes(j)
                            coef(i) = xcalc_heav(heavn(1), hea_se, heavn(5))*direct(j)
                        end if
!   AVEC L ENRICHISSEMENT SHIFTED [FRIES] LES FONCTIONS D ENRICHISSIMENT SINGULIERES
!      S ANNULENT AU NOEUD => RIEN A FAIRE
                    end if
                end do
!       CAS DDL_IMPO DX DY DZ (ET/OU PRE1 => POUR HM-XFEM ONLY)
            elseif (motcle .eq. 'DX' .or. motcle .eq. 'DY' .or. motcle .eq. 'DZ') then
!         COEFFICIENTS ET DDLS DE LA RELATION
                ddl(1) = 'D'//motcle(2:2)
                coef(1) = 1.d0
                i = 1
                if (nfiss .eq. 1) then
                    if (stano(1) .eq. 1 .or. stano(1) .eq. 3) then
                        i = i+1
                        ddl(i) = 'H1'//motcle(2:2)
                        coef(i) = xcalc_heav(heavn(1), hea_se, heavn(5))
                    end if
!   AVEC L ENRICHISSEMENT SHIFTED [FRIES] LES FONCTIONS D ENRICHISSIMENT SINGULIERES
!      S ANNULENT AU NOEUD => RIEN A FAIRE
                end if
            elseif (motcle .eq. 'PRE1') then
!           COEFFICIENTS ET DDLS DE LA RELATION
                ddl(1) = 'PRE1'
                coef(1) = 1.d0
                i = 1
                if (nfiss .eq. 1) then
                    if (stano(1) .eq. 1 .or. stano(1) .eq. 3) then
                        i = i+1
                        ddl(i) = 'H1'//motcle(1:4)
                        coef(i) = xcalc_heav(heavn(1), hea_se, heavn(5))
                    end if
                end if
            end if
            nterm = i
            call afrela(coef, [cbid], ddl, noeud, dimens, &
                        [0.d0], nterm, valimr, valimc, valimf, &
                        'REEL', fonree, 0.d0, lisrel)
        end do
!
    end if
!
    icompt = icompt+1
!
888 continue
!
    if (nfiss .gt. 1) call detrsd('CHAM_ELEM_S', ch4)
!
    call jedema()
end subroutine
