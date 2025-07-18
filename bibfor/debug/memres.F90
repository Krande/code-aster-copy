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

subroutine memres(limpr, ldyn, titre, prec, tmax)
    implicit none
#include "asterc/hpalloc.h"
#include "asterc/hpdeallc.h"
#include "asterc/loisem.h"
#include "asterfort/assert.h"
#include "asterfort/jjldyn.h"
#include "asterfort/utgtme.h"
    character(len=*) :: titre, limpr, ldyn
    real(kind=8) :: tmax, prec
! person_in_charge: jacques.pellet at edf.fr
!   BUT: CALCULER LA QUANTITE DE MEMOIRE ENCORE LIBRE
!-------------------------------------------------------------
! IN  LIMPR: /'NON' : ON SE CONTENTE DE CALCULER TMAX
!            /'OUI' : ON IMPRIME TMAX (+ 2 AUTRES INFOS)
!               IMPRIME SUR .MESS 3 TAILLES (EN MO) :
!                 * MEM JEVEUX DYNAMIQUE ALLOUEE
!                 * MEM JEVEUX INDISPENSABLE ("U")
!                 * MEM ENCORE DISPONIBLE (TMAX)
! IN  TITRE : CHAINE DE CARACTERES IMPRIMEE AU DEBUT DES LIGNES
!            SI LIMPR='OUI'
! IN  LDYN : /'OUI' : ON FORCE LA LIBERATION DE LA MEMOIRE JEVEUX
!               DYNAMIQUE QUI N'EST PAS "U". (APPEL A JJLDYN)
!            /'NON' : ON NE FORCE PAS LA LIBERATION DE LA MEMOIRE JEVEUX
!               DYNAMIQUE
! IN  PREC : PRECISION SOUHAITEE POUR TMAX (EN MO)
! OUT TMAX : TAILLE (EN MO) DE LA MEMOIRE ENCORE DISPONIBLE
!
! REMARQUES IMPORTANTES :
! ---------------------
! 1) CETTE ROUTINE NE PEUT FONCTIONNER QUE SI LE "SYSTEME" LIMITE LA
! MEMOIRE VIRTUELLE QUE L'ON PEUT ALLOUER.
! EN BATCH SUR BULL, C'EST RMS QUI LIMITE LA MEMOIRE
! EN INTERACTIF, IL FAUT SE FIXER SOI-MEME LA LIMITE EN FAISANT
! PAR EXEMPLE : ULIMIT -V 800000 (800 MO)
!
! 2) SUR LES MACHINES AVEC DES INTEGER*4, SI LA MEMOIRE DISPONIBLE
! EST TROP GRANDE, IL PEUT ARRIVER QUE L'ON DEPASSE LA CAPACITE DES I4.
! NORMALEMENT, ON EST ALORS ARRETE PAR UN ASSERT DANS MEMRES.F
!-------------------------------------------------------------
!
    real(kind=8) :: rval(2)
    character(len=8) :: k8tab(2)
    integer(kind=8) :: lon1, idep, incr, iret
    integer(kind=8) :: k, jad, ierr, nbfree, ltot
!
    ASSERT(limpr .eq. 'OUI' .or. limpr .eq. 'NON')
    ASSERT(ldyn .eq. 'OUI' .or. ldyn .eq. 'NON')
    ASSERT(prec .gt. 0.d0)
!
    tmax = -99.d0
    if (ldyn .eq. 'OUI') call jjldyn(0, -1, ltot)
!
    nbfree = 0
    idep = 0
    incr = int(prec*1024*1024/loisem()/2)
!
!
10  continue
    lon1 = incr
    do k = 1, 40
!
!        -- ON VERIFIE LA CAPACITE DES ENTIERS :
        lon1 = lon1*2
        ASSERT(lon1 .gt. 0)
        ASSERT(idep+lon1 .gt. 0)
!
!       -- TENTATIVE D'ALLOCATION :
        call hpalloc(jad, idep+lon1, ierr, 0)
        if (ierr .eq. 0) then
            call hpdeallc(jad, nbfree)
            goto 30
        else
            ASSERT(ierr .eq. -2)
            idep = idep+lon1/2
        end if
!
        if (k .gt. 1) then
            goto 10
        else
            tmax = dble(idep)
            goto 40
        end if
30      continue
    end do
!
40  continue
!
!      -- ON VERIFIE QUE TMAX EST ASSEZ BIEN EVALUE (+/- PREC):
    if (.true._1) then
        call hpalloc(jad, idep-2*incr, ierr, 0)
        ASSERT(ierr .eq. 0)
        call hpdeallc(jad, nbfree)
        call hpalloc(jad, idep+2*incr, ierr, 0)
        ASSERT(ierr .ne. 0)
    end if
!
!     -- CONVERSION EN MO :
    tmax = tmax*loisem()/1024/1024
!
!
    if (limpr .eq. 'OUI') then
        k8tab(1) = 'COUR_JV'
        k8tab(2) = 'CUSE_JV'
        call utgtme(2, k8tab, rval, iret)
        write (6, 9000) '<MEMRES>', titre, rval(1), rval(2), tmax
    end if
!
9000 format(a8, 1x, a, 3(1x, f12.2))
end subroutine
