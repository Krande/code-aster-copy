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
subroutine moinsr(j, n, idil, idiich, idsuiv, &
                  nosuiv, idip, noip, iilib, iimax)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
    integer(kind=8) :: j, n, idil, idiich, idsuiv, idip, iilib, iimax
    integer(kind=8) :: k, idebch, ii1, kil, kip
    character(len=*) :: nosuiv, noip
!
!     INSERTION DU TABLEAU D'ENTIERS ORDONNES IL DANS LA CHAINE ORDONNEE
!     J DE LA STRUCTURE CHAINEE:  IICH,SUIV,IP
!     -----------------------------------------------------------------
!     ATTENTION CE PROGRAMME CONTIENT DES POSSIBILITES D'EXTENSION
!     DES VECTEURS D'ENTIERS DE NOM KSUIV ET KIP.
!     -----------------------------------------------------------------
! IN  J         NUMERO DE LA CHAINE DANS LAQUELLE ON INSERE IL
! IN  N         LONGUEUR DU TABLEAU IL
! IN  IDIL      AD. JEVEUX DU TABLEAU IL QUE L'ON INSERE DANS LA
!               CHAINE J
! VAR IDIICH    AD. JEVEUX DE LA TABLE DES ADRESSES DANS LE VECTEUR DE
!               NOM KIP DU DEBUT DE LA CHAINE K.
!               SI IICH(K) <= 0 LA CHAINE K EST VIDE .
! VAR IDSUIV    AD. JEVEUX DU VECTEUR DE NOM NOSUIV.
! VAR NOSUIV    NOM DU TABLEAU DE CHAINAGE : SUIV(K) DONNE L'ADRESSE
!               DANS IP DE L'ELEMENT QUI SUIT L'ELEMENT K.
!               SI SUIV(K) <= 0 IP(K) EST LE DERNIER ELEMENT DE LA
!               CHAINE -SUIV(K)
! VAR IDIP      AD. JEVEUX DU VECTEUR DE NOM NOIP
! VAR NOIP      NOM DE LA TABLE DES ELEMENTS CHAINES
! VAR IILIB     NUMERO DE LA PREMIERE ADRESSE LIBRE DE IP
! VAR IIMAX     DIMENSION (IN: INITIALE, OUT :  FINALE) DES VECTEURS
!               SUIV ET IP QUI SONT AGRANDIS SI C'EST NECESSAIRE
!------------------------------------------------------------------
!- PRECAUTIONS D'EMPLOI:  N               > 0
!                         IL(I+1)         > IL(I)
!                         IP(IISUIV(K)) > IP(K)
!------------------------------------------------------------------
!     -------------------------------------------------------------
!     FONCTIONS JEVEUX
!     -------------------------------------------------------------
!     -------------------------------------------------------------
!     -------------------------------------------------------------
!----------------------------------------------------------------------
!                DEBUT DES INSTRUCTIONS
!----------------------------------------------------------------------
    if (zi(idiich-1+j) .le. 0) then
!
!
!        --- LA CHAINE J EST VIDE. ON L'INITIALISE PAR IL(1:N) ---
        zi(idiich-1+j) = iilib
        if ((iilib+n) .ge. iimax) then
            iimax = nint(1.5d0*iimax)
            ASSERT((iilib+n) .lt. iimax)
            call juveca(noip, iimax)
            call jeveuo(noip, 'E', idip)
            call juveca(nosuiv, iimax)
            call jeveuo(nosuiv, 'E', idsuiv)
        end if
        do k = 1, n
            zi(idsuiv-1+iilib) = iilib+1
            zi(idip-1+iilib) = zi(idil-1+k)
            iilib = iilib+1
        end do
!        MARQUAGE FIN DE CHAINE
        zi(idsuiv-1+iilib-1) = -j
!
    else
!---
!        LA CHAINE J EST NON VIDE:  DETERMINATION DES PREMIERS INDICES
!        KIL ET II1 TELS QUE IP(II1) < IL(KIL)
        idebch = zi(idiich-1+j)
        if (zi(idil) .lt. zi(idip-1+idebch)) then
!
!           INSERTION DE IL(1) EN DEBUT DE CHAINE
            zi(idiich-1+j) = iilib
            zi(idip-1+iilib) = zi(idil)
            zi(idsuiv-1+iilib) = idebch
            ii1 = iilib
            kil = 2
            if ((iilib+1) .ge. iimax) then
                iimax = nint(1.5d0*iimax)
                call juveca(noip, iimax)
                ASSERT((iilib+1) .lt. iimax)
                call jeveuo(noip, 'E', idip)
                call juveca(nosuiv, iimax)
                call jeveuo(nosuiv, 'E', idsuiv)
            end if
            iilib = iilib+1
!
        else if (zi(idil) .eq. zi(idip-1+idebch)) then
!
!           IL(1) EXISTE DEJA DANS LA CHAINE J . PAS D'INSERTION
            ii1 = idebch
            kil = 2
!
        else
!
!           IP(IDEBCH) < IL(1)
            ii1 = idebch
            kil = 1
        end if
!
!        INSERTION DU RESTE DE IL
        kip = zi(idsuiv-1+ii1)
!
20      continue
        if (kil .le. n) then
!
!           TOUS LES ELEMENTS DE IL N'ONT PAS ETE TRAITES
30          continue
            if (kip .gt. 0) then
!
!              LA CHAINE J N A PAS ETE ENTIEREMENT PARCOURUE.
!              INSERTION EVENTUELLE DE IL(KIL) EN MILIEU DE CHAINE.
!
!              ASSERTION : IL(KIL) > IP(II1)
!
                if (zi(idil-1+kil) .eq. zi(idip-1+kip)) then
!
!                  L'ELEMENT IL(KIL) EXISTE DEJA DANS LA CHAINE J
!                  PAS D'INSERTION
                    kil = kil+1
                    ii1 = kip
                    kip = zi(idsuiv-1+ii1)
                    goto 20
!
                else if (zi(idil-1+kil) .gt. zi(idip-1+kip)) then
!
!                  L'ELEMENT IL(KIL) NE S'INSERE PAS AVANT IP(KIP)
                    ii1 = kip
                    kip = zi(idsuiv-1+ii1)
                    goto 30
!
                else
!
!                  IP(II1) <IL(KIL) <IP(KIP) INSERTION DE IL(KIL)
!                  ENTRE CES 2 ELEMENTS
                    zi(idsuiv-1+ii1) = iilib
                    ii1 = iilib
                    zi(idsuiv-1+ii1) = kip
                    zi(idip-1+iilib) = zi(idil-1+kil)
                    if ((iilib+1) .ge. iimax) then
                        iimax = nint(1.5d0*iimax)
                        call juveca(noip, iimax)
                        ASSERT((iilib+1) .lt. iimax)
                        call jeveuo(noip, 'E', idip)
                        call juveca(nosuiv, iimax)
                        call jeveuo(nosuiv, 'E', idsuiv)
                    end if
                    iilib = iilib+1
                    kil = kil+1
                    goto 20
!
                end if
!
            else
!
!              LA CHAINE J A ETE ENTIEREMENT PARCOURUE.
!              INSERTION DES ELEMENT RESTANT DE IL EN FIN DE CHAINE  .
                zi(idsuiv-1+ii1) = iilib
                if ((iilib+1+n-kil) .ge. iimax) then
                    iimax = nint(1.5d0*iimax)
                    call juveca(noip, iimax)
                    ASSERT((iilib+1+n-kil) .lt. iimax)
                    call jeveuo(noip, 'E', idip)
                    call juveca(nosuiv, iimax)
                    call jeveuo(nosuiv, 'E', idsuiv)
                end if
                do k = kil, n
                    zi(idsuiv-1+iilib) = iilib+1
                    zi(idip-1+iilib) = zi(idil-1+k)
                    iilib = iilib+1
                end do
!              MARQUAGE FIN DE CHAINE
                zi(idsuiv-1+iilib-1) = -j
!
            end if
        end if
    end if
end subroutine
