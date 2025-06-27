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
subroutine rcvale(nommaz, phenom, nbpar, nompar, valpar, &
                  nbres, nomres, valres, icodre, iarret)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rccome.h"
#include "asterfort/rcvals.h"
#include "asterfort/utmess.h"
    integer(kind=8), intent(in) :: nbpar, nbres
    character(len=*), intent(in) :: phenom
    integer(kind=8), intent(in) :: iarret
    character(len=*), intent(in) :: nommaz
    integer(kind=8), intent(out) :: icodre(nbres)
    character(len=8), intent(in) :: nompar(nbpar)
    character(len=*), intent(in) :: nomres(nbres)
    real(kind=8), intent(in) :: valpar(nbpar)
    real(kind=8), intent(out) :: valres(nbres)
! ----------------------------------------------------------------------
!     OBTENTION DE LA VALEUR VALRES D'UN "ELEMENT" D'UNE RELATION DE
!     COMPORTEMENT D'UN MATERIAU DONNE (NOUVELLE FORMULE RAPIDE)
!
!     ARGUMENTS D'ENTREE:
!        NOMMAT : NOM UTILISATEUR DU MATERIAU
!        PHENOM : NOM DU PHENOMENE (I.E. MOT CLE FACTEUR)
!        NBPAR  : NOMBRE DE PARAMETRES DANS NOMPAR ET VALPAR
!        NOMPAR : NOMS DES PARAMETRES(EX: 'TEMP' )
!        VALPAR : VALEURS DES PARAMETRES
!        NBRES  : NOMBRE DE RESULTATS
!        NOMRES : NOM DES RESULTATS (EX: E,NU,... )
!                 TELS QU'IL FIGURENT DANS LA COMMANDE DEFI_MATERIAU
!        IARRET : 0/1/2  COMPORTEMENT VOULU EN CAS DE PROBLEME
!           /0 : PAS DE MESSAGE D'ERREUR
!           /1 : ERREUR FATALE AVEC IDENTIFICATION DE LA MAILLE
!           /2 : ERREUR FATALE SANS IDENTIFICATION DE LA MAILLE
!           PLUS PRECISEMENT :
!           * SI LE PARAMETRE (REEL, COMPLEXE OU FONCTION)
!             N'A PAS ETE FOURNI PAR L'UTILISATEUR :
!               SI IARRET > 0  => ERREUR FATALE
!               SI IARRET = 0  => CODE RETOUR > 0
!           * SI LE PARAMETRE EST UNE FONCTION FOURNIE PAR
!             L'UTILISATEUR MAIS QUE SON EVALUATION ECHOUE :
!             => ERREUR FATALE DANS TOUS LES CAS
!
!     ARGUMENTS DE SORTIE:
!       VALRES : VALEURS DES RESULTATS APRES RECUPERATION
!                OU EVALUATION DE LA FONCTION
!       ICODRE : POUR CHAQUE RESULTAT, 0 SI ON A TROUVE, 1 SINON
!
!
!
!
    integer(kind=8) :: nbmx, nbresp, ires, ier, nbr, nbc, nbk, iret
    integer(kind=8) :: nbobj, nbf, ir, ik
    parameter(nbmx=30)
    integer(kind=8) :: nbfp
    real(kind=8) :: valrep(nbmx)
    aster_logical :: change
    integer(kind=8) :: icodr2(nbmx)
    character(len=2) :: kstop
    character(len=32) :: nomphe, phen, phepre
    character(len=8) :: matpre
    character(len=16) :: nomrep(nbmx), nomfop(nbmx)
    character(len=11) :: k11
    character(len=8) :: nommat
    real(kind=8), pointer :: valr(:) => null()
    character(len=16), pointer :: valk(:) => null()
    save matpre, phepre, nbfp, nbresp, nomrep, valrep, icodr2, nomfop
!
    call jemarq()
    nommat = nommaz
    phen = phenom
    kstop = 'F '
    k11 = ''
!
    ASSERT(iarret .ge. 0 .and. iarret .le. 2)
!
!
! --- TESTS: CELA A-T-IL CHANGE ?
    change = .false.
    if (nbres .gt. nbmx) then
        call utmess('F', 'MODELISA6_94', sk=nommat)
    end if
    if (nommat .ne. matpre) change = .true.
    if (phen .ne. phepre) change = .true.
    if (nbres .ne. nbresp) change = .true.
    do ires = 1, nbres
        if (nomres(ires) .ne. nomrep(ires)) change = .true.
    end do
!
!
    if (.not. change) then
        do ires = 1, nbres
            valres(ires) = valrep(ires)
            icodre(ires) = icodr2(ires)
        end do
        if (nbfp .eq. 0) goto 9999
!
        do ires = 1, nbres
            if (nomfop(ires) .ne. ' ') then
                call fointe(kstop, nomfop(ires), nbpar, nompar, valpar, &
                            valres(ires), ier)
                ASSERT(ier .eq. 0)
                icodre(ires) = 0
            end if
        end do
!
!
    else
        nomphe = phen
        call rccome(nommat, nomphe, iret, k11_ind_nomrc=k11)
        call jeexin(nommat//k11//'.VALR', iret)
        if (iret .eq. 0) then
            do ires = 1, nbres
                icodre(ires) = 1
            end do
            goto 999
        end if
!
        call jeveuo(nommat//k11//'.VALR', 'L', vr=valr)
        call jelira(nommat//k11//'.VALR', 'LONUTI', nbr)
        call jelira(nommat//k11//'.VALC', 'LONUTI', nbc)
        call jeveuo(nommat//k11//'.VALK', 'L', vk16=valk)
        call jelira(nommat//k11//'.VALK', 'LONUTI', nbk)
        do ires = 1, nbres
            icodre(ires) = 1
            nomfop(ires) = ' '
        end do
        nbobj = 0
        do ir = 1, nbr
            do ires = 1, nbres
                if (nomres(ires) .eq. valk(ir)) then
                    valres(ires) = valr(ir)
                    icodre(ires) = 0
                    nbobj = nbobj+1
                end if
            end do
        end do
        if (nbobj .ne. nbres) then
            nbf = (nbk-nbr-nbc)/2
            do ires = 1, nbres
                do ik = 1, nbf
                    if (nomres(ires) .eq. valk(nbr+nbc+ik)) then
                        nomfop(ires) = valk(nbr+nbc+nbf+ik)
                        call fointe(kstop, nomfop(ires), nbpar, nompar, valpar, &
                                    valres(ires), ier)
                        ASSERT(ier .eq. 0)
                        icodre(ires) = 0
                    end if
                end do
            end do
        else
            nbf = 0
        end if
!
!       -- SAUVEGARDE DES VALEURS POUR LE PROCHAIN APPEL :
        matpre = nommat
        phepre = phen
        nbfp = nbf
        nbresp = nbres
        do ires = 1, nbresp
            nomrep(ires) = nomres(ires)
            valrep(ires) = valres(ires)
            icodr2(ires) = icodre(ires)
        end do
!
    end if
999 continue
9999 continue
!
    call rcvals(iarret, icodre, nbres, nomres)
!
    call jedema()
end subroutine
