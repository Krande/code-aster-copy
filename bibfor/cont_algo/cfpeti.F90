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

subroutine cfpeti(sdcont_solv, neq, nbliai, nbliac, &
                  rho, llliai, llliac)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/caladu.h"
#include "asterfort/cfelpv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
!

    character(len=24) :: sdcont_solv
    integer(kind=8) :: neq, nbliai
    integer(kind=8) :: nbliac
    real(kind=8) :: rho
    integer(kind=8) :: llliai, llliac
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - RESOLUTION)
!
! VERIFICATION SI L'ENSEMBLE DES LIAISONS SUPPOSEES EST TROP PETIT
!
! ----------------------------------------------------------------------
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NEQ    : NOMBRE D'EQUATIONS
! IN  NBLIAC : NOMBRE DE LIAISONS ACTIVES
! IN  NBLIAI : NOMBRE DE LIAISONS
! OUT RHO    : COEFFICIENT DE MISE A JOUR
!               1   - TOUTES LES LIAISONS SONT ACTIVES
!               VAL - VALEUR A CORRIGER SUR LE JEU POUR LA LIAISON LLMIN
! OUT LLLIAI : NUMERO DE LA LIAISON LA PLUS VIOLEE
! OUT LLLIAC : NUMERO DE LA LIAISON _ACTIVE_ LA PLUS VIOLEE
!
!
!
!
    real(kind=8), parameter :: un = 1.d0
    real(kind=8) :: rhorho
    real(kind=8) :: aadelt, jeuold, jeunew, jeuinc
    aster_logical :: liaiac, delpos, lelpiv
    integer(kind=8) :: btotal, iliai, iliac
    character(len=19) :: liac
    integer(kind=8) :: jliac
    character(len=24) :: apcoef, apddl, appoin
    integer(kind=8) :: japcoe, japddl, japptr
    character(len=24) :: jeuite
    integer(kind=8) :: jjeuit
    character(len=19) :: ddeplc, ddelt
    integer(kind=8) :: nbddl, jdecal
    real(kind=8), pointer :: vddelt(:) => null()
    real(kind=8), pointer :: ddepc(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- ACCES STRUCTURES DE DONNEES DE CONTACT
!
    liac = sdcont_solv(1:14)//'.LIAC'
    appoin = sdcont_solv(1:14)//'.APPOIN'
    apcoef = sdcont_solv(1:14)//'.APCOEF'
    apddl = sdcont_solv(1:14)//'.APDDL'
    jeuite = sdcont_solv(1:14)//'.JEUITE'
    call jeveuo(liac, 'L', jliac)
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apcoef, 'L', japcoe)
    call jeveuo(apddl, 'L', japddl)
    call jeveuo(jeuite, 'L', jjeuit)
!
! --- ACCES AUX CHAMPS DE TRAVAIL
! --- DDEPLC: INCREMENT DE SOLUTION APRES CORRECTION DU CONTACT
! --- DDELT : INCREMENT DE SOLUTION ITERATION DE CONTACT
!
    ddeplc = sdcont_solv(1:14)//'.DELC'
    ddelt = sdcont_solv(1:14)//'.DDEL'
    call jeveuo(ddeplc(1:19)//'.VALE', 'L', vr=ddepc)
    call jeveuo(ddelt(1:19)//'.VALE', 'L', vr=vddelt)
!
! --- INITIALISATIONS
!
    rhorho = r8maem()
    delpos = .false.
    llliai = 0
    btotal = nbliac
!
! --- VERIFICATION : ENSEMBLE DES LIAISONS SUPPOSEES TROP PETIT ?
!
    if (nbliac .eq. nbliai) then
!
! ----- TOUTES LES LIAISONS SONT ACTIVES
!
        rhorho = un
    else if (nbliac .lt. nbliai) then
!
! ----- RECHERCHE DES LIAISONS NON ACTIVES
!
        do iliai = 1, nbliai
!
            liaiac = .false.
!
! ------- LA LIAISON ILIAI EST-ELLE ACTIVE ? (-> LIAIAC)
!
            do iliac = 1, btotal
                if (zi(jliac-1+iliac) .eq. iliai) liaiac = .true.
            end do
!
! ------- CALCUL DE RHOMIN
!
            if (.not. liaiac) then
!
! --------- CALCUL DE [A].{DDELT} SI LA LIAISON N'EST PAS ACTIVE
!
                jdecal = zi(japptr+iliai-1)
                nbddl = zi(japptr+iliai)-zi(japptr+iliai-1)
                call caladu(neq, nbddl, zr(japcoe+jdecal), zi(japddl+jdecal), vddelt, &
                            aadelt)
!
! --------- SI [A].{DDELT} EST POSITIF: LIAISON ACTIVEE
!
                if (aadelt .gt. r8prem()) then
!
! ----------- ON NE PREND PAS EN COMPTE UNE LIAISON A PIVOT NUL
!
                    call cfelpv(iliai, sdcont_solv, nbliai, lelpiv)
                    if (lelpiv) then
                        exit
                    else
                        delpos = .true.
!
! ------------- VALEUR DE {JEU(DEPTOT)} : JEU AVANT ITERATION DE NEWTON
!
                        jeuold = zr(jjeuit+3*(iliai-1)+1-1)
!
! ------------- CALCUL DE [A].{DDEPLC} - CORRECTION DU JEU
!
                        call caladu(neq, nbddl, zr(japcoe+jdecal), zi(japddl+jdecal), ddepc, &
                                    jeuinc)
!
! ------------- CALCUL DE {JEU(DEPTOT) - [A].{DDEPLC}}/[A].{DDELT}
!
                        jeunew = jeuold-jeuinc
                        jeunew = jeunew/aadelt
!
! ------------- RHOMIN = MIN({JEU(DEPTOT) - A.DDEPLC}/{{A.DDELT})
!
                        if (jeunew .lt. rhorho) then
                            rhorho = jeunew
!
! --------------- LLLIAI: LIAISON LA PLUS VIOLEE
!
                            llliai = iliai
                        end if
                    end if
                end if
            end if
        end do
!
! ----- TOUS LES {A.DELTA} SONT NEGATIFS
!
        if (.not. delpos) then
            rhorho = un
        end if
    end if
!
! --- ON FAIT EN SORTE QUE RHORHO <= 1.D0
!
    rho = min(rhorho, un)
!
! --- LIAISON ACTIVE
!
    if (llliai .ne. 0) then
        llliac = nbliac+1
    end if
!
    call jedema()
!
end subroutine
