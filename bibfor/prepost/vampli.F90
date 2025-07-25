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
subroutine vampli(vwork, tdisp, liste, nbt, nbordr, &
                  numini, nbp, tspaq, nomopt, cxsr)
! person_in_charge: jean-michel.proix at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rvinvt.h"
!
    integer(kind=8) :: tdisp, nbp, liste(nbp), nbt, nbordr, numini
    integer(kind=8) :: tspaq
    real(kind=8) :: vwork(tdisp)
    character(len=16) :: nomopt
    character(len=19) :: cxsr
! ---------------------------------------------------------------------
! BUT: CALCULER LA VARIATION D'AMPLITUDE MAXIMALE
! ---------------------------------------------------------------------
! ARGUMENTS:
! VWORK     IN    R  : VECTEUR DE TRAVAIL CONTENANT
!                      L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                      ATTACHES A CHAQUE POINT DE GAUSS DES MAILLES
!                      DU <<PAQUET>> DE MAILLES.
! TDISP     IN    I  : DIMENSION DU VECTEUR VWORK
! LISTE     IN    I  : LISTE COMPLETE DES NOEUDS OU DES POINTS DE GAUSS
!                      PAR MAILLE A TRAITER.
! NBT       IN    I  : NOMBRE TOTAL DE POINT DE GAUSS OU DE NOEUDS
!                      A TRAITER.
! NBORDR    IN    I  : NOMBRE DE NUMERO D'ORDRE STOCKE DANS LA
!                      STRUCTURE DE DONNEES RESULTAT.
! NUMINI    IN    I  : NUMERO DE LA 1ERE MAILLE DU <<PAQUET>> DE
!                      MAILLES COURANT OU DU 1ER NOEUD DU <<PAQUET>> DE
!                      NOEUDS COURANT.
! NBP       IN    I  : NOMBRE DE MAILLES DANS LE <<PAQUET>> DE
!                      MAILLES COURANT.OU NOMBRE DE NOEUDS DANS LE
!                      <<PAQUET>> DE NOEUDS COURANT.
! TSPAQ     IN    I  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                      COURANT.
! NOMOPT    IN    K16: POST-TRAITEMENT AUX NOEUDS OU AUX POINTS DE GAUSS
! CXSR      IN    K19: NOM DU CHAMP SIMPLE DESTINE A RECEVOIR LES
!                      RESULTATS :
!                           X = N ==> CNSR = RESULTATS AUX NOEUDS
!                           X = E ==> CESR = RESULTATS AUX ELEMENTS
!
! REMARQUE :
!  - LA TAILLE DU SOUS-PAQUET EST EGALE A LA TAILLE DU <<PAQUET>> DE
!    MAILLES DIVISEE PAR LE NOMBRE DE NUMERO D'ORDRE (NBORDR).
!-----------------------------------------------------------------------
    integer(kind=8) :: nnoini, nbnop, nbnot, jcnrd, jcnrl
    integer(kind=8) :: l, cnbno, kwork, somnow, ibidno, nunoe, inop
    integer(kind=8) :: decal, i, j, adrsi, adrsj, k, icmp
    integer(kind=8) :: jad, nmaini, nbmap, nbpgt, jcerd, jcerl
    integer(kind=8) :: nbpg, nbpgp, sompgw, imap, ipg
!
    real(kind=8) :: vresu(24)
    real(kind=8) :: vavmis, vatres, vmis, tres, trac, detr
    real(kind=8) :: tensi(6), tensj(6), dtens(6)
!
    character(len=19) :: cnsr, cesr
    real(kind=8), pointer :: cesv(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
!
!-----------------------------------------------------------------------
!
    call jemarq()
!
    vavmis = 0.0d0
    vatres = 0.0d0
!
    if (nomopt .eq. 'DOMA_NOEUD') then
!
        cnsr = cxsr
        nnoini = numini
        nbnop = nbp
        nbnot = nbt
!
!  OBTENTION DES ADRESSES '.CNSD', '.CNSL' ET '.CNSV' DU CHAMP SIMPLE
!  DESTINE A RECEVOIR LES RESULTATS : VMIS ET TRESCA
!
        call jeveuo(cnsr//'.CNSD', 'L', jcnrd)
        call jeveuo(cnsr//'.CNSL', 'E', jcnrl)
        call jeveuo(cnsr//'.CNSV', 'E', vr=cnsv)
!
!
!   CONSTRUCTION DU VECTEUR : CONTRAINTE = F(NUMERO D'ORDRE) EN CHAQUE
!   NOEUDS DU PAQUET DE MAILLES.
        l = 1
        cnbno = 0
        kwork = 0
        somnow = 0
        ibidno = 1
!
!  BOUCLE SUR LES NOEUDS
!
        do inop = nnoini, nnoini+(nbnop-1)
!
            nunoe = liste(inop)
!
            if (inop .gt. nnoini) then
                kwork = 1
                somnow = somnow+1
            end if
!
            cnbno = cnbno+1
            if ((l*int(nbnot/10.0d0)) .lt. cnbno) then
                l = l+1
            end if
!
!
!  CALCUL DE LA VARIATION D'AMPLITUDE
!
!
!  IL Y A 6 COMPOSANTES POUR LES CONTRAINTES ==> DECAL=6
            decal = 18
!
!
            do i = 1, (nbordr-1)
!
                do j = (i+1), nbordr
!
                    adrsi = (i-1)*tspaq+kwork*somnow*decal+(ibidno-1)*decal
!
                    adrsj = (j-1)*tspaq+kwork*somnow*decal+(ibidno-1)*decal
!
!   TENSI/J(1) = CPXXI/J   TENSI/J(2) = CPYYI/J   TENSI/J(3) = CPZZI/J
!   TENSI/J(4) = CPXYI/J   TENSI/J(5) = CPXZI/J   TENSI/J(6) = CPYZI/J
!
                    tensi(1) = vwork(adrsi+1)
                    tensi(2) = vwork(adrsi+2)
                    tensi(3) = vwork(adrsi+3)
                    tensi(4) = vwork(adrsi+4)
                    tensi(5) = vwork(adrsi+5)
                    tensi(6) = vwork(adrsi+6)
!
                    tensj(1) = vwork(adrsj+1)
                    tensj(2) = vwork(adrsj+2)
                    tensj(3) = vwork(adrsj+3)
                    tensj(4) = vwork(adrsj+4)
                    tensj(5) = vwork(adrsj+5)
                    tensj(6) = vwork(adrsj+6)
!
!
                    do k = 1, 6
                        dtens(k) = tensi(k)-tensj(k)
                    end do
!
                    call rvinvt(dtens, vmis, tres, trac, detr)
!
!
!
                    if (vmis .gt. vavmis) then
                        vavmis = vmis
                    end if
!
                    if (tres .gt. vatres) then
                        vatres = tres
                    end if
!
                end do
!
            end do
!
!
            do icmp = 1, 24
                vresu(icmp) = 0.0d0
            end do
            vresu(23) = vavmis
            vresu(24) = vatres
!
!  AFFECTATION DES RESULTATS DANS UN CHAM_ELEM SIMPLE
!
            do icmp = 1, 24
                jad = 24*(nunoe-1)+icmp
                zl(jcnrl-1+jad) = .true.
                cnsv(jad) = vresu(icmp)
            end do
!
        end do
!
!
!  POUR LES GROUPES DE MAILLES
!
    else if (nomopt .eq. 'DOMA_ELGA') then
!
        cesr = cxsr
        nmaini = numini
        nbmap = nbp
        nbpgt = nbt
!
!
!  OBTENTION DES ADRESSES '.CESD', '.CESL' ET '.CESV' DU CHAMP SIMPLE
!  DESTINE A RECEVOIR LES RESULTATS : DOMMAGE_MAX, COORDONNEES VECTEUR
!  NORMAL CORRESPONDANT
!
        call jeveuo(cesr//'.CESD', 'L', jcerd)
        call jeveuo(cesr//'.CESL', 'E', jcerl)
        call jeveuo(cesr//'.CESV', 'E', vr=cesv)
!
!
!  CONSTRUCTION DU VECTEUR : CISAILLEMENT = F(NUMERO D'ORDRE) EN CHAQUE
!  POINT DE GAUSS DU PAQUET DE MAILLES.
        l = 1
        nbpg = 0
        nbpgp = 0
        kwork = 0
        sompgw = 0
!
! BOUCLE SUR LES MAILLES
!
        do imap = nmaini, nmaini+(nbmap-1)
            if (imap .gt. nmaini) then
                kwork = 1
                sompgw = sompgw+liste(imap-1)
            end if
            nbpg = liste(imap)
!
! SI LA MAILLE COURANTE N'A PAS DE POINTS DE GAUSS, LE PROGRAMME
! PASSE DIRECTEMENT A LA MAILLE SUIVANTE.
            if (nbpg .eq. 0) then
                goto 100
            end if
!
            nbpgp = nbpgp+nbpg
            if ((l*int(nbpgt/10.0d0)) .lt. nbpgp) then
                l = l+1
            end if
!
!  BOUCLE SUR LES POINTS DE GAUSS
!
            do ipg = 1, nbpg
!
!  CALCUL DE LA VARIATION D'AMPLITUDE
!
!  IL Y A 6 COMPOSANTES POUR LES CONTRAINTES ==> DECAL=6
                decal = 18
!
!
!  BOUCLE SUR LES NUMEROS D'ORDRES
!
                do i = 1, (nbordr-1)
!
                    do j = (i+1), nbordr
!
                        adrsi = (i-1)*tspaq+kwork*sompgw*decal+(ipg-1)*decal
!
                        adrsj = (j-1)*tspaq+kwork*sompgw*decal+(ipg-1)*decal
!
!
!   TENSI/J(1) = CPXXI/J   TENSI/J(2) = CPYYI/J   TENSI/J(3) = CPZZI/J
!   TENSI/J(4) = CPXYI/J   TENSI/J(5) = CPXZI/J   TENSI/J(6) = CPYZI/J
!
                        tensi(1) = vwork(adrsi+1)
                        tensi(2) = vwork(adrsi+2)
                        tensi(3) = vwork(adrsi+3)
                        tensi(4) = vwork(adrsi+4)
                        tensi(5) = vwork(adrsi+5)
                        tensi(6) = vwork(adrsi+6)
!
                        tensj(1) = vwork(adrsj+1)
                        tensj(2) = vwork(adrsj+2)
                        tensj(3) = vwork(adrsj+3)
                        tensj(4) = vwork(adrsj+4)
                        tensj(5) = vwork(adrsj+5)
                        tensj(6) = vwork(adrsj+6)
!
!
                        do k = 1, 6
                            dtens(k) = tensi(k)-tensj(k)
                        end do
!
                        call rvinvt(dtens, vmis, tres, trac, detr)
!
!
                        if (vmis .gt. vavmis) then
                            vavmis = vmis
                        end if
!
                        if (tres .gt. vatres) then
                            vatres = tres
                        end if
!
                    end do
!
                end do
!
! 11. CONSTRUCTION D'UN CHAM_ELEM SIMPLE PUIS D'UN CHAM_ELEM CONTENANT
!     POUR CHAQUE POINT DE GAUSS DE CHAQUE MAILLE LE DOMMAGE_MAX ET LE
!     VECTEUR NORMAL ASSOCIE.
!
                do icmp = 1, 24
                    vresu(icmp) = 0.0d0
                end do
                vresu(23) = vavmis
                vresu(24) = vatres
!
! 12. AFFECTATION DES RESULTATS DANS UN CHAM_ELEM SIMPLE
!
                do icmp = 1, 24
                    call cesexi('C', jcerd, jcerl, imap, ipg, &
                                1, icmp, jad)
!
                    ASSERT(jad .ne. 0)
                    jad = abs(jad)
                    zl(jcerl-1+jad) = .true.
                    cesv(jad) = vresu(icmp)
!
                end do
!
            end do
!
100         continue
        end do
!
    end if
!
!
! MENAGE
!
! PAS DE MENAGE DANS CETTE ROUTINE
!
    call jedema()
end subroutine
