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

subroutine ldc_disjvp(ppr, ppi, ppc, yy0, dy0, dyy, decoup)
    implicit none
#include "asterf_types.h"
    real(kind=8)     :: ppr(*)
    integer(kind=8)          :: ppi(*)
    character(len=*) :: ppc(*)
    real(kind=8)     :: yy0(*)
    real(kind=8)     :: dy0(*)
    real(kind=8)     :: dyy(*)
    aster_logical    :: decoup
!
! person_in_charge: jean-luc.flejou at edf.fr
! ----------------------------------------------------------------------
!
#include "asterc/r8prem.h"
!
! ----------------------------------------------------------------------
!        MODELE ELASTO PLASTIQUE ENDOMMAGEANT EN FLEXION (1DDL=RZ)
!
!  IN
!     ppr      : paramètres réels
!     ppi      : paramètres entiers
!     ppc      : paramètres character
!     nbeq     : nombre d'équations
!     yy0      : valeurs initiales
!     dy0      : dérivées initiales
!     pf       : adresse des fonctions
!
!  OUT
!     dyy      : dérivées calculées
!     decoup   : pour forcer l'adaptation du pas de temps
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: mm, xm, sim, my, fp, re, dre, hrem, hrep, hdrp, hdrm, kplus, kmoins, lamb
    real(kind=8) :: factp, factm, fact, dm1, dm2, dd1, dd2, eps
!
!   système d'équations :
    integer(kind=8) :: imoment, itheta, ithetap, idp, idm, ixm, idiss
    parameter(imoment=1, itheta=2, ithetap=3, idp=4, idm=5, ixm=6, idiss=7)
!   paramètres du modèle : KE, KP, KDP, KDM, RDP, RDM, MYP, MYM
    integer(kind=8) :: ike, ikp, ikdp, ikdm, irdp, irdm, imyp, imym
    parameter(ike=1, ikp=2, ikdp=3, ikdm=4, irdp=5, irdm=6, imyp=7, imym=8)
!
    decoup = ASTER_FALSE
    eps = r8prem()
!   initialisation
    dyy(itheta) = dy0(itheta)
!   par defaut rien n evolue
    dyy(ithetap) = 0.0
    dyy(idp) = 0.0
    dyy(idm) = 0.0
    dyy(ixm) = 0.0
    dyy(imoment) = 0.0
    dyy(idiss) = 0.0
!   raideurs secantes
    if (abs(yy0(idp)) .lt. eps .or. abs(yy0(idm)) .lt. eps) then
        decoup = ASTER_TRUE
        goto 999
    end if
    kplus = ppr(ikdp)+(ppr(ike)-ppr(ikdp))*ppr(irdp)/yy0(idp)
    kmoins = ppr(ikdm)+(ppr(ike)-ppr(ikdm))*abs(ppr(irdm))/yy0(idm)
!   rotation elastique
    re = yy0(itheta)-yy0(ithetap)
!   fonctions de heaviside
    if (re < 0.) then
        hrem = 1.
        hrep = 0.
    else
        hrem = 0.
        hrep = 1.
    end if
!   calcul du moment - projection elastique
    mm = kplus*re*hrep+kmoins*re*hrem
!   calcul du seuil
    xm = yy0(ixm)
!   signe de M-Xm : sim
    if (mm < xm) then
        sim = -1.0
        my = abs(ppr(imym))
    else
        sim = 1.0
        my = ppr(imyp)
    end if
!   fonction seuil : fp
    fp = sim*(mm-xm)-my
    if ((dyy(itheta)-dyy(ithetap)) < 0.) then
        hdrm = 1.
        hdrp = 0.
    else
        hdrm = 0.
        hdrp = 1.
    end if
!   cas plastique
    if (fp > 0.) then
        factp = ppr(ikdp)+(ppr(ike)-ppr(ikdp))*(ppr(irdp)/yy0(idp))*(1.+re*hdrp/yy0(idp))
        factm = ppr(ikdm)+(ppr(ike)-ppr(ikdm))*(abs(ppr(irdm))/yy0(idm))*(1.-re*hdrm/yy0(idm))
        fact = hrep*factp+hrem*factm
        if ((ppr(ikp)+fact) > eps) then
            lamb = sim*dyy(itheta)*fact/(ppr(ikp)+fact)
            if (lamb > 0.) then
                dyy(ithetap) = sim*lamb
            end if
        end if
        dyy(ixm) = ppr(ikp)*dyy(ithetap)
    end if
!   derivee de la rotation elastique
    dre = dyy(itheta)-dyy(ithetap)
!   mise a jour de l endommagement
    if (re > yy0(idp)) then
        dyy(idp) = max(dre, 0.)
    end if
    if ((-1.0*(re)) > yy0(idm)) then
        dyy(idm) = max(-1.0*dre, 0.)
    end if
    dm1 = kplus*dre*hrep-(ppr(ike)-ppr(ikdp))*(ppr(irdp)/(yy0(idp)*yy0(idp)))*dyy(idp)*re*hrep
    dm2 = kmoins*dre*hrem-(ppr(ike)-ppr(ikdm))*(abs(ppr(irdm))/(yy0(idm)*yy0(idm)))*dyy(idm)*re*hrem
    dyy(imoment) = dm1+dm2
    dd1 = (1./2.)*(ppr(ike)-ppr(ikdp))*ppr(irdp)*dyy(idp)
    dd2 = (1./2.)*(ppr(ike)-ppr(ikdm))*abs(ppr(irdm))*dyy(idm)+my*abs(dyy(ithetap))
    dyy(idiss) = dd1+dd2
!
999 continue
end subroutine
