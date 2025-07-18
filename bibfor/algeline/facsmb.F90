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

subroutine facsmb(nbnd, nbsn, supnd, invsup, parent, &
                  xadj, adjncy, anc, nouv, fils, &
                  frere, local, global, adress, lfront, &
                  nblign, lgsn, debfac, debfsn, chaine, &
                  place, nbass, delg, lgind, ier)
! person_in_charge: olivier.boiteau at edf.fr
! aslint: disable=W1504
    implicit none
#include "jeveux.h"
#include "asterfort/infniv.h"
#include "asterfort/inschn.h"
#include "asterfort/mltalc.h"
    integer(kind=8) :: nbnd, nbsn, lgind, nbnd1
    integer(kind=8) :: supnd(nbsn+1), invsup(nbnd), parent(nbsn)
    integer(kind=8) :: xadj(nbnd+1), adjncy(*)
    integer(kind=8) :: anc(nbnd), nouv(nbnd), fils(nbsn), frere(nbsn), delg(nbnd)
    integer(kind=4) :: global(lgind), local(lgind)
    integer(kind=8) :: adress(nbsn+1), debfsn(nbsn+1)
    integer(kind=8) :: lfront(nbsn), nblign(nbsn), lgsn(nbsn), debfac(nbnd+1), ier
    integer(kind=8) :: debut
!
!================================================================
!     FACTORISATION SYMBOLIQUE POUR LA MULTIFRONTALE
!     PRISE EN COMPTE DES DDL LAGRANGE  D'ASTER
!     UTILISATION D'UNE LISTE CHAINEE ORDONNEE
!     POUR RANGER LA LISTE  DES NOEUDS ET DES VOISINS D'UN  SUPERNOEUD
!----------------------------------------------
!
!     3 CAS
!     I) LE SND A SON PREMIER DDL NON LAGRANGE : CAS STANDARD ,
!     ON MET SES VOISINS DANS LA CHAINE. SI UN  DES NDS DU SND
!     EST DANS UNE RELATION LINEAIRE,ON MET LAMBDA2
!     DANS LA CHAINE (79)
!     II) LE PREMIER EST UN DDL LAGRANGE DE BLOQUAGE
!     ON MET DANS LA CHAINE
!     TOUS LES LAGRANGES DU SUPERNOEUD,
!     (LBD1 DE BLOCS ET TOUS LES LBD2)
!     ON CHERCHE PRMNDI 1ER ND NON LAGRANGE DU SND
!     ,ET ON MET PRMNDI ET SES VOISINS DANS LA CHAINE
!     III) C'EST UN SND LBD1 DE RELATION LINEAIRE:
!     IL EST MIS DANS LA CHAINE,
!     AINSI QUE LES NDS DE LA RELATION LINEAIRE ET LE LAMBDA2
!-------------------------------------------------------------
!     MODIFICATION DU 15/09/98
!     AMDBAR A L' AIR DE FAIRE DANS CERTAINS CAS DE L'AMALGAMME I.E.
!     CERTAINS SUPERNOEUDS SONT REGROUPES EN UN SEUL AU PRIX DE CERTAINS
!     ZEROS. IL FAUT ALORS METTRE DANS LA LISTE CHAINEE, LES VOISINS DE
!     TOUS LES NOEUDS DU SN ET NON PLUS CEUX DU PREMIER NOEUD COMME
!     AUPPARAVANT. (ON POUVAIT EN OUBLIER)
!     LE TABLEAU PLACE SERT DE FLAG DANS CETTE INSERTION
!     AFIN DE GAGNER UN PEU DE TEMPS, CAR LA PLUPART DES INSERTIONS
!     SONT REDONDANTES.
!---------------------------------------------------------------
!     SOUS-PROGRAMME APPELLE : INSCHN
!==================================================================
    integer(kind=8) :: chaine(nbnd), place(nbnd), nbass(nbsn)
    integer(kind=8) :: i, k, j, nd, p, sni, andi, sn, suiv, cour
    integer(kind=8) :: ind, ndk, ndi, dli
    integer(kind=8) :: ifm, niv, long, decal
!
!-----------------------------------------------------------------------
    call infniv(ifm, niv)
    ier = 0
!     CALCUL DE INVSUP FILS ET FRERE
    do i = 1, nbsn
        place(i) = 0
        fils(i) = 0
        frere(i) = 0
        lgsn(i) = supnd(i+1)-supnd(i)
        do j = supnd(i), supnd(i+1)-1
            invsup(j) = i
        end do
    end do
    do nd = 1, nbsn
        p = parent(nd)
        if (p .ne. 0) then
            if (fils(p) .eq. 0) then
                fils(p) = nd
                place(p) = nd
            else
                frere(place(p)) = nd
                place(p) = nd
            end if
        end if
    end do
!
!
    decal = 0
    adress(1) = 1
    debfac(1) = 1
    do sni = 1, nbsn
        lgsn(sni) = supnd(sni+1)-supnd(sni)
    end do
!
    adress(1) = 1
    debfac(1) = 1
    nbnd1 = nbnd+1
!     CORRECTION DE NOV 2006 CETTE BOUCLE 311 REMPLACE LA PRECEDENTE
!     INTERNE A LA BOUCLE 310 CELA ENTRAINAIT BEAUCOUP DE TEMPS CPU
    do i = 1, nbnd
        place(i) = 0
    end do
    do ndi = 1, nbsn
        chaine(ndi) = nbnd1
    end do
    do sni = 1, nbsn
        ndi = supnd(sni)
        andi = anc(ndi)
        chaine(ndi) = nbnd1
        dli = ndi
        debut = ndi
!        LES TERMES INITIAUX DE LA MATRICE SONT MIS DANS LA CHAINE
!        QUE LE DDL ORDINAIRE OU LAGRANGE
        call inschn(andi, dli, xadj, adjncy, chaine, &
                    nouv, place, debut)
        if (delg(andi) .eq. 0) then
!--------------------------------------------------------------
!     1 .....................................   NDI EST UN DDL ORDINAIRE
!
!     ON INSERE AUSSI DANS LA CHAINE LES VOISINS INITIAUX DE TOUTES
!     INCONNNUES DU SUPERNOEUD, AU CAS OU LA RENUMEROTATION
!     FASSE DE L'AMALGAME.
!--------------------------------------------------------------
!
            do dli = ndi+1, supnd(sni+1)-1
                andi = anc(dli)
                if (delg(andi) .ne. 0) goto 151
                call inschn(andi, dli, xadj, adjncy, chaine, &
                            nouv, place, debut)
!
            end do
151         continue
        end if
!-------------------------------------------------------------
!
!     LES NOEUDS VOISINS DES FILS SONT MIS DANS LA CHAINE
!-------------------------------------------------------------
        sn = fils(sni)
230     continue
!        DO WHILE (SN.NE.0)
        if (sn .ne. 0) then
!           K = ADRESS(SN) + LGSN(SN) + 1 CORRECTION DU 15/03/02
            k = adress(sn)+lgsn(sn)
            ind = 1
            nd = ndi
240         continue
!     DO WHILE (K.LT.ADRESS(SN+1))
            if (k .lt. adress(sn+1)) then
                ndk = global(k)
                if (ndk .gt. ndi) then
                    suiv = nd
235                 continue
                    if (suiv .lt. ndk) then
!     DO WHILE(SUIV.LT.NDK)
                        cour = suiv
                        suiv = chaine(cour)
                        goto 235
                    end if
                    if (suiv .gt. ndk) then
                        chaine(cour) = ndk
                        chaine(ndk) = suiv
                        place(ndk) = 1
                    end if
                    nd = ndk
                end if
                k = k+1
                goto 240
!     FIN DO WHILE
            end if
            sn = frere(sn)
            goto 230
!     FIN DO WHILE
        end if
        k = 0
        ind = ndi
!     DO WHILE (IND.NE.NBND1) ( FIN DE LA CHAINE)
280     continue
        if (ind .ne. nbnd1) then
!-------------------------------------------------------------
!     VERIFICATION DE LA LONGUEUR DE GLOBAL
!-------------------------------------------------------------
            if ((adress(sni)+k) .gt. lgind) then
                ier = lgind*2
                if (niv .ge. 2) then
                    write (ifm, *)&
     &             'LONGUEUR DE GLOBAL PEUT ETRE INSUFFISANTE '
                    write (ifm, *) 'LONGUEUR ALLOUEE :', lgind
                    write (ifm, *) 'ON REITERE AVEC :', ier
                end if
                goto 999
            end if
            global(k+adress(sni)) = int(ind, 4)
            place(global(k+adress(sni))) = k+1
            k = k+1
            ind = chaine(ind)
            goto 280
!     FIN DO WHILE
        end if
        adress(sni+1) = k+adress(sni)
!...........................................
        sn = fils(sni)
!     DO WHILE (SN.NE.0)
290     continue
        if (sn .ne. 0) then
            call mltalc(local, global, adress, sn, lgsn, &
                        place, sni, supnd, nbass(sn))
!
            sn = frere(sn)
            goto 290
!     FIN DO WHILE
        end if
        nblign(sni) = adress(sni+1)-adress(sni)
        lfront(sni) = nblign(sni)-lgsn(sni)
        long = nblign(sni)
!     ANCIENNE VERSION SANS DGEMV
!     DO 300 K = SUPND(SNI),SUPND(SNI+1) - 1
!     DEBFAC(K+1) = DEBFAC(K) + L
!     L = L - 1
!     300     CONTINUE
!     MODIFS POUR DGEMV
        do k = 1, lgsn(sni)
            nd = supnd(sni)+k-1
            debfac(nd) = decal+k
            decal = decal+long
        end do
        debfsn(sni) = debfac(supnd(sni))
!   ON REMET LE TABLEAU PLACE A ZERO ICI AU LIEU DE 311
        do k = adress(sni), (adress(sni+1)-1)
            place(global(k)) = 0
        end do
!
    end do
!
    debfac(nbnd+1) = decal+1
    debfsn(nbsn+1) = debfac(nbnd+1)
    if (niv .ge. 2) then
        write (ifm, *) '   --- LONGUEUR DE LA MATRICE FACTORISEE ', decal
    end if
!
!
999 continue
!
end subroutine
