! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine nmceai(numedd, depdel, deppr1, deppr2, depold, &
                  sdpilo, rho, eta, isxfe, f, &
                  indic)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer :: indic
    character(len=24) :: numedd
    character(len=19) :: sdpilo, depdel, depold, deppr1, deppr2
    real(kind=8) :: eta, rho, f
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PILOTAGE - SELECTION PARAMETRE)
!
! CALCUL DU PARAMETRE DE SELECTION DE TYPE ANGL_INCR_DEPL
!
! ----------------------------------------------------------------------
!
!
! IN  NUMEDD : NUME_DDL
! IN  SDPILO : SD PILOTAGE
! IN  DEPDEL : INCREMENT DE DEPLACEMENT DEPUIS DEBUT PAS DE TEMPS
! IN  DEPOLD : INCREMENT DE DEPLACEMENT PAS DE TEMPS PRECEDENT
! IN  DEPPR1 : INCREMENT DE DEPLACEMENT K-1.F_DONNE
! IN  DEPPR2 : INCREMENT DE DEPLACEMENT K-1.F_PILO
! IN  RHO    : PARAMETRE DE RECHERCHE LINEAIRE
! IN  ETA    : PARAMETRE DE PILOTAGE
! IN  ISXFE  : INDIQUE SI LE MODELE EST UN MODELE XFEM
! OUT F      : VALEUR DU CRITERE
! OUT INDIC  : 0 CRITERE NON UTILISABLE
!              1 CRITERE UTILISABLE
!
!
!
!
    real(kind=8) :: sca, nodup, coef, nodup1, nodup2
    real(kind=8) :: norm_depold
    integer :: jdu1
    integer :: neq, i, j
    character(len=19) :: profch, chapil, chapic, selpil
    real(kind=8) :: dn, dc, dp, da
    aster_logical :: isxfe
    real(kind=8), pointer :: coee(:) => null()
    real(kind=8), pointer :: vcoef(:) => null()
    real(kind=8), pointer :: depde(:) => null()
    real(kind=8), pointer :: depol(:) => null()
    real(kind=8), pointer :: du0(:) => null()
    real(kind=8), pointer :: plsl(:) => null()
    integer, pointer :: deeq(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    if (isxfe) then
        chapil = sdpilo(1:14)//'.PLCR'
        call jeveuo(chapil(1:19)//'.VALE', 'L', vr=vcoef)
        chapic = sdpilo(1:14)//'.PLCI'
        call jeveuo(chapic(1:19)//'.VALE', 'L', vr=coee)
    else
!
! --- ACCES VECTEUR DE SELCTION CMP DX/DY/DZ
!
        selpil = sdpilo(1:14)//'.PLSL'
        call jeveuo(selpil(1:19)//'.VALE', 'L', vr=plsl)
    end if
!
!
! --- INITIALISATIONS--------------------------------
!
    sca = 0.d0
    nodup = 0.d0
    nodup1 = 0.d0
    nodup2 = 0.d0
    norm_depold = 0.d0
    f = 0.d0
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
    call dismoi('NUME_EQUA', depdel, 'CHAM_NO', repk=profch)
    call jeveuo(profch(1:19)//'.DEEQ', 'L', vi=deeq)
!
! --- ACCES AUX VECTEURS SOLUTIONS
!
    call jeveuo(depdel(1:19)//'.VALE', 'L', vr=depde)
    call jeveuo(deppr1(1:19)//'.VALE', 'L', vr=du0)
    call jeveuo(deppr2(1:19)//'.VALE', 'L', jdu1)
    call jeveuo(depold(1:19)//'.VALE', 'L', vr=depol)
!
!
! --- CALCUL DE L'ANGLE
!
    if (isxfe) then
        do i = 1, neq
            if (deeq(2*i) .gt. 0) then
                norm_depold = norm_depold+depol(i)**2
                if (coee(i) .eq. 0.d0) then
                    sca = sca+depol(i)*vcoef(i)**2*(depde(i)+rho*du0(1+i-1)+eta*zr(jdu1+i-&
                          &1))
                    nodup1 = nodup1+vcoef(i)**2*(depde(i)+rho*du0(i)+eta*zr(jdu1+i-1))**2
                    nodup2 = nodup2+vcoef(i)**2*depol(i)**2
                else
                    da = 0.d0
                    dn = 0.d0
                    dc = 0.d0
                    dp = 0.d0
                    do j = i+1, neq
                        if (coee(i) .eq. coee(j)) then
                            da = da+vcoef(i)*depol(i)+vcoef(j)*depol(j)
                            dn = dn+vcoef(i)*depde(i)+vcoef(j)*depde(j)
                            dc = dc+vcoef(i)*du0(i)+vcoef(j)*du0(j)
                            dp = dp+vcoef(i)*zr(jdu1-1+i)+vcoef(j)*zr(jdu1-1+j)
                        end if
                    end do
                    sca = sca+da*(dn+rho*dc+eta*dp)
                    nodup1 = nodup1+(dn+rho*dc+eta*dp)**2
                    nodup2 = nodup2+da**2
                end if
            end if
        end do
        nodup = nodup1*nodup2
    else
        do i = 1, neq
            coef = plsl(i)
            sca = sca+(depol(i)*(depde(i)+rho*du0(i)+eta*zr(jdu1+i-1)))*coef
            nodup = nodup+(depde(i)+rho*du0(i)+eta*zr(jdu1+i-1))**2
            norm_depold = norm_depold+depol(i)**2
        end do
    end if
!
    if (nodup .eq. 0.d0 .or. norm_depold .eq. 0.d0) then
        indic = 0
        f = 0.d0
    else
        indic = 1
        f = sca/sqrt(nodup)
    end if
    f = -f
!
    call jedema()
end subroutine
