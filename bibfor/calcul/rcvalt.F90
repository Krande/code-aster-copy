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

subroutine rcvalt(fami, kpg, ksp, poum, jmat, nomat, mfact, &
                  nbpar, nompar, valpar, &
                  nbres, valres, icodre, iarret)

    use calcul_module, only: ca_jvcnom_, ca_nbcvrc_

    implicit none
! person_in_charge: jacques.pellet at edf.fr

#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/fointa.h"
#include "asterc/r8nnem.h"

    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in)          :: kpg
    integer(kind=8), intent(in)          :: ksp
    character(len=1), intent(in) :: poum
    integer(kind=8), intent(in)          :: jmat
    character(len=*), intent(in) :: nomat, mfact
    integer(kind=8), intent(in)          :: nbpar
    character(len=*), intent(in) :: nompar(nbpar)
    real(kind=8), intent(in)     :: valpar(nbpar)
    integer(kind=8), intent(in)          :: nbres
    integer(kind=8), intent(out)         :: icodre(*)
    real(kind=8), intent(out)    :: valres(*)
    integer(kind=8), intent(in)          :: iarret

! ------------------------------------------------------------------------
!
!  But : Recuperation de l'ENSEMBLE des parametres materiaux sous un mot cle facteur
!        La liste des parametres retournes (et la liste des code_retour) est
!        ordonnee selon l'ordre du catalogue (mot cle ORDRE_PARAM)
!        Les nbres premiers parametres de ORDRE_PARAM sont retournes.
!        Si un parametre n'est pas renseigne dans DEFI_MATERIAU, on retourne :
!          * icodre(i)=1
!          * valres(i)= NaN

!     arguments d'entree:
!       in   fami    : famille de point de gauss (rigi,mass,...)
!       in   kpg,ksp : numero du (sous)point de gauss
!       in   poum    : '+' /'-'
!       in   jmat    : adresse de la liste des materiaux codes (zi(imate))
!       in   nomat   : nom du materiau dans le cas d'une liste de materiaux
!                      si = ' ', on exploite le premier de la liste
!       in   mfact   : nom du mot cle facteur (ex 'VISCOCHAB')
!       in   nbpar   : nombre de parametres dans nompar et valpar
!       in   nompar  : noms des parametres(ex: 'TEMP' )
!       in   valpar  : valeurs des parametres
!       in   nbres   : dimension de valres et icodre
!       in   iarret  : comportement souhaite si mfact n'est pas renseigne :
!                = 0 : on remplit icodre(*)=1 et on sort sans message.
!                = 1 : on arrete <F> en indiquant le nom de la maille.
!                = 2 : idem que 1 mais on n'indique pas la maille.

!     arguments de sortie:
!       out  valres(*)  : valeur reelle du parametre (ou NaN)
!       out  icodre(*)  : 0 si on a trouve, 1 sinon
! ----------------------------------------------------------------------
    integer(kind=8) :: lmat, icomp, ipi, ipif, iadzi, iazk24, nbk, ivalk, nbr, nbc
    integer(kind=8) :: lfct, imat, nbmat, n1, k, posi, nmcs
    integer(kind=8) :: ier, jordr, jkord, ivalr, kr, kc, kf, ivalc
    integer(kind=8) :: nbpamx, nbpar2, nbpart, ipar, kmat, inom
    parameter(nbpamx=10)
    real(kind=8) ::  valvrc, rundf, valeur, valpa2(nbpamx)
    character(len=8) :: nompa2(nbpamx), novrc, nomi
    parameter(lmat=9, lfct=10)
    character(len=32) :: valk
    character(len=32) :: nomphe
!  ---------------------------------------------------------------------

    nomphe = mfact
    rundf = r8nnem()

!   -- on ne gere pas encore lsup
    ASSERT(nomphe .ne. 'TRACTION')
    ASSERT(nomphe .ne. 'META_TRAC')

!   -- initialisation de icodre(*) et valres(*) :
!   ---------------------------------------------
    do k = 1, nbres
        icodre(k) = 1
        valres(k) = rundf
    end do

!   -- Calcul de imat
!      Si nomat est fourni , on explore l'entete de la sd mater_code pour
!      trouver le "bon" materiau de la la liste
!   ----------------------------------------------------------------------
    nbmat = zi(jmat)
    if (nomat(1:1) .ne. ' ') then
        do kmat = 1, nbmat
            inom = zi(jmat+kmat)
            nomi = zk8(inom)
            if (nomi .eq. nomat) then
                imat = jmat+zi(jmat+nbmat+kmat)
                goto 9
            end if
        end do
        call utmess('F', 'CALCUL_45', sk=nomat)
    else
        ASSERT(nbmat .eq. 1)
        imat = jmat+zi(jmat+nbmat+1)
    end if
9   continue

!   -- calcul de ipi (pour nomphe):
!   -------------------------------
    do icomp = 1, zi(imat+1)
        if (nomphe .eq. zk32(zi(imat)+icomp-1)) then
            ipi = zi(imat+2+icomp-1)
            goto 11
        end if
    end do

!   -- selon la valeur de iarret on arrete ou non :
!   -----------------------------------------------
    if (iarret .ge. 1) then
        valk = nomphe
        call utmess('F+', 'CALCUL_46', sk=valk)
        if (iarret .eq. 1) then
            call tecael(iadzi, iazk24)
            call utmess('F+', 'CALCUL_47', si=zi(iadzi))
        end if
        call utmess('F', 'VIDE_1')
    end if
    goto 999

!   -- calcul de valres(*) :
!   -------------------------
11  continue
    nbr = zi(ipi)
    nbc = zi(ipi+1)
!   -- la routine n'a pas d'argument pour rendre des complexes :
    ASSERT(nbc .eq. 0)
    nbk = zi(ipi+2)

    ivalk = zi(ipi+3)
    ivalr = zi(ipi+4)
    ivalc = zi(ipi+5)

!   -- la routine n'a de sens que pour un mot cle facteur ayant le mot cle ORDRE_PARAM :
    jordr = zi(ipi+6)
    jkord = zi(ipi+7)
    ASSERT(jordr .ne. 1)
    ASSERT(jkord .ne. 1)

    n1 = zi(jkord-1+1)
    nmcs = zi(jkord-1+2)

!   -- si nbres est insuffisant :
    ASSERT(nbres .le. n1)

!   -- Recuperation des variables de commande :
!   -------------------------------------------
    nbpar2 = 0
    do ipar = 1, ca_nbcvrc_
        novrc = zk8(ca_jvcnom_-1+ipar)
        call rcvarc(' ', novrc, poum, fami, kpg, &
                    ksp, valvrc, ier)
        if (ier .eq. 0) then
            nbpar2 = nbpar2+1
            nompa2(nbpar2) = novrc
            valpa2(nbpar2) = valvrc
        end if
    end do

!   -- On ajoute les varc au debut de la liste des parametres
!      car fointa donne priorite aux derniers :
!   ----------------------------------------------------------
    nbpart = nbpar+nbpar2
    ASSERT(nbpart .le. nbpamx)
    do ipar = 1, nbpar
        nompa2(nbpar2+ipar) = nompar(ipar)
        valpa2(nbpar2+ipar) = valpar(ipar)
    end do

!   -- Boucle sur les mots cles simples et calcul de leurs valeurs :
!   ----------------------------------------------------------------
    do k = 1, nmcs
!       -- posi : numero dans .ORDR :
        posi = zi(jkord-1+2+k)

        kr = zi(jkord-1+2+nmcs+k)
        kc = zi(jkord-1+2+2*nmcs+k)
        kf = zi(jkord-1+2+3*nmcs+k)

        if (kr .gt. 0) then
            valeur = zr(ivalr-1+kr)
        elseif (kc .gt. 0) then
            ASSERT(.false.)
        elseif (kf .gt. 0) then
!           -- c'est un concept : une fonction ou une liste
            ipif = ipi+lmat+(kf-1)*lfct-1

!           -- cas d'une fonction :
            if (zi(ipif+9) .eq. 1) then
                call fointa(ipif, nbpart, nompa2, valpa2, valeur)

!           -- cas d'une table TRC :
            elseif (zi(ipif+9) .eq. 2) then
                ASSERT(.false.)

!           -- cas d'une liste de reels :
            elseif (zi(ipif+9) .eq. 3) then
                ASSERT(.false.)

!           -- cas d'une liste de fonctions :
            elseif (zi(ipif+9) .eq. 4) then
                ASSERT(.false.)
!               ipif2=zi(ipif)+3+lfct*(kv-1) ???
!               call fointa(ipif2, nbpart, nompa2, valpa2, valeur)
            else
                ASSERT(.false.)
            end if
        else
            ASSERT(.false.)
        end if
        valres(posi) = valeur
        icodre(posi) = 0
    end do

999 continue

end subroutine
