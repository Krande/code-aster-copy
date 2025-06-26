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

subroutine rcvala(jmat, nomat, phenom, nbpar, nompar, &
                  valpar, nbres, nomres, valres, icodre, &
                  iarret, nan)
    use calcul_module, only: ca_iactif_
    implicit none
#include "jeveux.h"
#include "asterfort/fointa.h"
#include "asterfort/rcvals.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
#include "asterc/r8nnem.h"
! ----------------------------------------------------------------------
! person_in_charge: jacques.pellet at edf.fr
    integer(kind=8), intent(in) :: jmat, nbpar, nbres, iarret
    real(kind=8), intent(in) :: valpar(nbpar)
    real(kind=8), intent(out) :: valres(nbres)
    integer(kind=8), intent(out) :: icodre(nbres)
    character(len=*), intent(in) :: nomat, phenom, nompar(nbpar), nomres(nbres)
    character(len=3), intent(in), optional :: nan
!
! But : Recuperation des valeurs d'une liste de coefficients d'une relation de
!       comportement pour un materiau donne.
!
!     arguments d'entree:
!        jmat      : adresse de la liste des materiaux codes
!        nomat     : nom du materiau dans le cas d'une liste de materiaux
!                    si = ' ', on exploite le premier de la liste
!        phenom    : nom du phenomene (mot cle facteur sans le "_FO")
!        nbpar     : nombre de parametres dans nompar(*) et valpar(*)
!        nompar(*) : noms des parametres(ex: 'TEMP', 'INST' )
!        valpar(*) : valeurs des parametres
!        nbres     : nombre de coefficients recherches
!                    (dimension des tableaux nomres(*), valres(*) et icodre(*)
!        nomres(*) : nom des resultats (ex: 'E','NU',... )
!                    tels qu'il figurent dans la commande DEFI_MATERIAU
!       iarret = 0 : on remplit icodre et on sort sans message.
!              = 1 : si un des parametres n'est pas trouve, on arrete
!                       en fatal en indiquant le nom de la maille.
!              = 2 : idem que 1 mais on n'indique pas la maille.
!              = 3 : si le phenomene n'est pas trouve pour l'element, on fait
!                    idem que iarret = 0, si le phenomene est trouve mais pas
!                    le parametre demande, on fait idem que irarret = 1
!       nan    = 'OUI' (defaut) : pour les parametres non trouves, on retourne valres = NaN
!              = 'NON' : pour les parametres non trouves, on ne modifie pas valres
!
!     arguments de sortie:
!        valres(*) : valeurs des resultats apres recuperation et interpolation
!        icodre(*) : pour chaque resultat, 0 si on a trouve, 1 sinon
!
! ----------------------------------------------------------------------
!   -- parameters associes au materiau code :
    integer(kind=8) :: lmat, lfct, lsup
    parameter(lmat=9, lfct=10, lsup=2)

    integer(kind=8) :: ires, icomp, ipi, iadzi, iazk24, nbobj, nbr, nbc, nbf, ivalk
    integer(kind=8) :: ivalr, ir, ipif, ik, nbmat, imat, kmat, inom
    character(len=8) :: nomail, nomi
    character(len=24) :: nomphe
    character(len=24) :: valk(2)
    real(kind=8) :: rundf
    aster_logical :: lnan
!  ---------------------------------------------------------------------

!   -- On est oblige de recopier phenom car il faut le tronquer
!      parfois a 10 avant de le comparer
    nomphe = phenom

!   -- initialisation de icodre(*) et valres(*) :
!   ---------------------------------------------
    rundf = r8nnem()
    lnan = .true.
    if (present(nan)) then
        ASSERT(nan .eq. 'OUI' .or. nan .eq. 'NON')
        if (nan .eq. 'NON') lnan = .false.
    end if
    do ires = 1, nbres
        icodre(ires) = 1
        if (lnan) valres(ires) = rundf
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
        if (nbmat .gt. 1) then
            if (ca_iactif_ .eq. 2) then
                valk(1) = nomres(1)
                call utmess('A', 'MODELISA9_4', nk=1, valk=valk)
            elseif (ca_iactif_ .eq. 1 .or. ca_iactif_ .eq. 3) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                valk(1) = nomail
                valk(2) = nomres(1)
                call utmess('A', 'MODELISA9_3', nk=2, valk=valk)
            else
                ASSERT(.false.)
            end if
        end if
        imat = jmat+zi(jmat+nbmat+1)
    end if
9   continue

!   -- calcul de ipi (pour nomphe):
!   -------------------------------
    do icomp = 1, zi(imat+1)
        if (nomphe .eq. zk32(zi(imat)+icomp-1)) then
            ipi = zi(imat+2+icomp-1)
            goto 22
        end if
    end do

!   -- selon la valeur de iarret on arrete ou non :
!   ----------------------------------------------
    if (iarret .ge. 1) then
        valk(1) = nomphe
        if (iarret .eq. 1) then
            if (ca_iactif_ .eq. 2) then
                call utmess('F', 'MODELISA9_73', nk=1, valk=valk)
            elseif (ca_iactif_ .eq. 1 .or. ca_iactif_ .eq. 3) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                valk(2) = nomail
                call utmess('F', 'MODELISA9_75', nk=2, valk=valk)
            else
                ASSERT(.false.)
            end if
        elseif (iarret .eq. 2) then
            call utmess('F', 'MODELISA9_74', sk=valk(1))
        elseif (iarret .eq. 3) then
            goto 888
        end if
    end if
    goto 999

!   -- calcul de valres(*) :
!   -------------------------
22  continue
    nbobj = 0
    nbr = zi(ipi)
    nbc = zi(ipi+1)
    nbf = zi(ipi+2)
    ivalk = zi(ipi+3)
    ivalr = zi(ipi+4)
    do ires = 1, nbres
        do ir = 1, nbr
            if (nomres(ires) .eq. zk16(ivalk+ir-1)) then
                valres(ires) = zr(ivalr-1+ir)
                icodre(ires) = 0
                nbobj = nbobj+1
                goto 32
            end if
        end do
32      continue
    end do

    if (nbobj .ne. nbres) then
        do ires = 1, nbres
            ipif = ipi+lmat-1
            do ik = 1, nbf
                if (nomres(ires) .eq. zk16(ivalk+nbr+nbc+ik-1)) then
                    ASSERT(zi(ipif+9) .eq. 1)
                    call fointa(ipif, nbpar, nompar, valpar, valres(ires))
                    icodre(ires) = 0
                end if
                ipif = ipif+lfct
                if (nomphe .eq. 'TRACTION') then
                    ipif = ipif+lsup
                else if (nomphe .eq. 'META_TRACT') then
                    ipif = ipif+lsup
                end if
            end do
        end do
    end if

999 continue

    call rcvals(iarret, icodre, nbres, nomres)

888 continue
end subroutine
