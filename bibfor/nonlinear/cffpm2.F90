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
subroutine cffpm2(resoco, resigr, nbliai, nbliac, ndim)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/calapr.h"
#include "asterfort/cfcglt.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/r8inir.h"
!
    character(len=24) :: resoco
    real(kind=8) :: resigr
    integer(kind=8) :: nbliai, nbliac, ndim
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (RESOLUTION - PENALISATION)
!
! CALCUL DE LA MATRICE FRO2 = E_T*AT
!
! ----------------------------------------------------------------------
!
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  RESIGR : RESI_GLOB_RELA
! IN  NBLIAI : NOMBRE DE LIAISONS DE CONTACT POSSIBLES
! IN  NBLIAC : NOMBRE DE LIAISONS ACTIVES
! IN  NDIM   : DIMENSION DU PROBLEME
!
!
!
!
    integer(kind=8) :: ndlmax
    parameter(ndlmax=30)
    integer(kind=8) :: jdecal, nbddl
    real(kind=8) :: jeuini, glis
    real(kind=8) :: coefpt, coefff, coefte, beta
    real(kind=8) :: lambdc, lambdf
    integer(kind=8) :: iliac, iliai, iliai2
    character(len=19) :: mu, liac, afmu
    integer(kind=8) :: jmu, jliac, jafmu
    character(len=24) :: appoin, apddl
    integer(kind=8) :: japptr, japddl
    character(len=24) :: jeux
    integer(kind=8) :: jjeux
    character(len=24) :: tacfin
    integer(kind=8) :: jtacf
    integer(kind=8) :: ztacf
    character(len=19) :: fro2
    integer(kind=8) :: jfro2
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    appoin = resoco(1:14)//'.APPOIN'
    apddl = resoco(1:14)//'.APDDL'
    liac = resoco(1:14)//'.LIAC'
    tacfin = resoco(1:14)//'.TACFIN'
    afmu = resoco(1:14)//'.AFMU'
    mu = resoco(1:14)//'.MU'
    fro2 = resoco(1:14)//'.FRO2'
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apddl, 'L', japddl)
    call jeveuo(liac, 'L', jliac)
    call jeveuo(tacfin, 'L', jtacf)
    call jeveuo(afmu, 'L', jafmu)
    call jeveuo(mu, 'E', jmu)
    ztacf = cfmmvd('ZTACF')
!
    jeux = resoco(1:14)//'.JEUX'
    call jeveuo(jeux, 'L', jjeux)
!
! --- CALCUL DE LA MATRICE E_T*AaT
!
    do iliai = 1, nbliai
!
! ----- INITIALISATION DE LA COLONNE
!
        call jeveuo(jexnum(fro2, iliai), 'E', jfro2)
        call r8inir(ndlmax, 0.d0, zr(jfro2), 1)
!
! ----- LA LIAISON EST-ELLE ACTIVE ?
!
        jeuini = zr(jjeux+3*(iliai-1)+1-1)
!
! ----- CALCUL
!
        if (jeuini .lt. 0.d0) then
!
! ------- REPERAGE DE LA LIAISON
!
            jdecal = zi(japptr+iliai-1)
            nbddl = zi(japptr+iliai)-zi(japptr+iliai-1)
!
! ------- CALCUL DU JEU TANGENT
!
            call cfcglt(resoco, iliai, glis)
!
! ------- PARAMETRES
!
            coefff = zr(jtacf+ztacf*(iliai-1)+1-1)
            coefpt = zr(jtacf+ztacf*(iliai-1)+3-1)
            coefte = zr(jtacf+ztacf*(iliai-1)+4-1)
!
! ------- LAMBDA DE FROTTEMENT
!
            do iliac = 1, nbliac
                iliai2 = zi(jliac+iliac-1)
                lambdc = zr(jmu+iliac-1)
                if (iliai2 .eq. iliai) then
                    if (lambdc .gt. 0.d0) then
                        lambdf = coefff*lambdc
                    else
                        lambdf = 0.d0
                    end if
                end if
            end do
!
! ------- ACTIVATION GLISSEMENT/ADHERENCE
!
            if (lambdf .eq. 0.d0) then
                beta = 0.d0
            else
                if (zr(jmu+2*nbliai+iliai-1) .ne. 0.d0) then
                    if (glis .le. lambdf/coefpt) then
                        beta = 0.d0
                    else
                        beta = sqrt(1.d0/(lambdf*glis))
                    end if
                else
                    beta = 0.d0
                end if
                if (resigr .ge. 1.d-03) then
                    beta = sqrt(coefte)*beta
                end if
                call calapr(nbddl, beta, zr(jafmu), zi(japddl+jdecal), zr(jfro2))
                zr(jmu+2*nbliai+iliai-1) = 1.d0
            end if
        else
            zr(jmu+2*nbliai+iliai-1) = 0.d0
        end if
        call jelibe(jexnum(fro2, iliai))
    end do
!
    call jedema()
!
end subroutine
