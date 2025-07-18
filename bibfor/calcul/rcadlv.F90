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

subroutine rcadlv(fami, kpg, ksp, poum, jmat, nomat, mfact, msimp, &
                  nbpar, nompar, valpar, jadr, nbres, icodre, iarret)

    use calcul_module, only: ca_jvcnom_, ca_nbcvrc_

    implicit none

#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcvals.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/fointa.h"
! -----------------------------------------------------------------
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in)          :: kpg
    integer(kind=8), intent(in)          :: ksp
    character(len=1), intent(in) :: poum
    integer(kind=8), intent(in)          :: jmat
    character(len=*), intent(in) :: nomat, mfact, msimp
    integer(kind=8), intent(in)          :: nbpar
    character(len=*), intent(in) :: nompar(nbpar)
    real(kind=8), intent(in)     :: valpar(nbpar)
    integer(kind=8), intent(out)         :: icodre
    integer(kind=8), intent(out)         :: jadr
    integer(kind=8), intent(out)         :: nbres
    integer(kind=8), intent(in)          :: iarret

! -----------------------------------------------------------------
!  Recuperation de l'adresse jeveux (dans zr) des coefficients materiau
!  correspondant a DEFI_MATERIAU/MFACT/MSIMP
!  Cette routine doit etre utilisee quand le mot cle MFACT/MSIMP
!  correspond a une liste de reels (UMAT/LISTE_COEF par exemple)
!
!     Arguments d'entree:
!       in   fami    : famille de point de gauss ('RIGI','MASS',...)
!       in   kpg,ksp : numero du (sous)point de gauss
!       in   poum    : '+' /'-'
!       in   jmat    : adresse de la liste des materiaux codes (zi(imate))
!       in   nomat   : nom du materiau dans le cas d'une liste de materiaux
!                      si = ' ', on exploite le premier de la liste
!       in   mfact   : nom du mot cle facteur (ex 'UMAT')
!       in   msimp   : nom du mot cle simple (ex 'LISTE_COEF')
!                      tels qu'il figurent dans la commande DEFI_MATERIAU
!       in   nbpar   : nombre de parametres dans nompar et valpar
!       in   nompar  : noms des parametres(ex: 'TEMP' )
!       in   valpar  : valeurs des parametres
!
!     Arguments de sortie:
!       out  jadr    : adresse dans zr de la liste de reels
!       out  nbres   : nombre de valeurs dans zr(jadr)
!       out  icodre  : 0 si on a trouve, 1 sinon
! ----------------------------------------------------------------------

    integer(kind=8) :: lmat, icomp, ipi, ipif, iadzi, iazk24, nbk, ivalk, ik, nbr, nbc
    integer(kind=8) :: lfct, imat, nbmat, code, kv, nbv, ipif2, kmat, inom
    real(kind=8) :: valres
    integer(kind=8) :: nbpamx, nbpar2, nbpart, ipar, ier
    parameter(nbpamx=15)
    real(kind=8) :: valpa2(nbpamx), valvrc
    character(len=8) :: nompa2(nbpamx), novrc
    parameter(lmat=9, lfct=10)
    character(len=32) :: valk
    character(len=8) :: nomi
    character(len=32) :: nomphe
!   --------------------------------------------------------------------------

    icodre = 1
    nbres = 0
    nomphe = mfact

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

!   -- calcul de jadr et de zr(jadr:)
!   -----------------------------------------------
11  continue
    nbr = zi(ipi)
    nbc = zi(ipi+1)
    nbk = zi(ipi+2)
    ivalk = zi(ipi+3)
    do ik = 1, nbk
        if (msimp .eq. zk16(ivalk+nbr+nbc+ik-1)) then
            icodre = 0
            ipif = ipi+lmat+(ik-1)*lfct-1
            ASSERT(zi(ipif+9) .eq. 3 .or. zi(ipif+9) .eq. 4)
            code = zi(zi(ipif))
            ASSERT(code .eq. -1 .or. code .eq. -2)

!           -- 1. Cas d'une liste de reels  :
!           ----------------------------------
            if (code .eq. -1) then
                jadr = zi(zi(ipif)+1)

!           -- 2. Cas d'une liste de fonctions
!           ----------------------------------
            else
                jadr = zi(zi(ipif)+1)
                nbv = zi(zi(ipif)+2)

!               -- 2.1 Recuperation des variables de commande :
!               -----------------------------------------------
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

!               -- 2.2 On ajoute les varc au debut de la liste des parametres
!                  car fointa donne priorite aux derniers :
!               --------------------------------------------------------------
                nbpart = nbpar+nbpar2
                ASSERT(nbpart .le. nbpamx)
                do ipar = 1, nbpar
                    nompa2(nbpar2+ipar) = nompar(ipar)
                    valpa2(nbpar2+ipar) = valpar(ipar)
                end do
!               -- 2.3 On evalue les fonctions de la liste :
!               ----------------------------------------------------------
                do kv = 1, nbv
                    ipif2 = zi(ipif)+3+lfct*(kv-1)
                    call fointa(ipif2, nbpart, nompa2, valpa2, valres)
                    zr(jadr-1+1+kv) = valres
                end do
            end if
            goto 999
        end if
    end do

    call rcvals(iarret, [icodre], 1, msimp)

999 continue
    if (icodre .eq. 0) then
        nbres = nint(zr(jadr))
        jadr = jadr+1
    end if
!
end subroutine
