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
subroutine rsingu(ndim, nelem, nbr, nalpha, degre, &
                  prec, erreur, alpha, types, re)
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndim, nelem, nbr(nelem), nalpha, degre
    real(kind=8) :: prec, erreur(nelem), alpha(nelem), re(nelem)
    character(len=16) :: types
!
!     BUT:
!         CALCUL DU RAPPORT ENTRE L ANCIENNE ET LA NOUVELLE TAILLE
!         OPTION : 'SING_ELEM'
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   NELEM         : NOMBRE D ELEMENTS FINIS
! IN   NBR(NELEM)    : NOMBRE DE COMPOSANTES A STOCKER PAR EF
!      3 SI EF SURFACIQUES EN 2D OU VOLUMIQUES EN 3D
!      0 SINON
! IN   NALPHA        : NOMBRE DE CPE PAR ELEMENT DIFFERENTS
!                      1 PAR DEFAUT SI PAS DE SINGULARITE
! IN   DEGRE         : DEGRE DES EF 1 EF P1 2 POUR EF P2
! IN   PREC          : % DE L ERREUR TOTALE SOUHAITE POUR CALCULER
!                      LA NOUVELLE CARTE DE TAILLE DES EF
! IN   ERREUR(NELEM) : ERREUR DE CHAQUE EF
! IN   ALPHA(NELEM)  : DEGRE DE LA SINGULARITE PAR ELEMENT
!
!      SORTIE :
!-------------
! OUT  RE(NELEM)     : RAPPORT ENTRE ANCIENNE ET NOUVELLE TAILLE
!
! ......................................................................
!
    integer(kind=8) :: inel, iter
    real(kind=8) :: errtot, prec0, cumm, ordre, ae(nelem), be(nelem)
    real(kind=8) :: d, mu, fonc, dfonc
    aster_logical :: lqi
!
! 0 - ERREUR EN NORME DE L ENERGIE OU EN QUANTITE D INTERET
!
    lqi = .false.
    if (types(1:2) .eq. 'QI') lqi = .true.
!
! 1 - CALCUL DE L ERREUR TOTALE
!
    errtot = 0.d0
    do inel = 1, nelem
        errtot = errtot+erreur(inel)**2
! ------POUR L ERREUR EN QUANTITE D INTERET QUI PEUT ETRE NEGATIVE
        erreur(inel) = abs(erreur(inel))
    end do
    errtot = sqrt(errtot)
!
! 2 - CALCUL DE RE
!
    ordre = degre
    d = ndim
    prec0 = prec*errtot
!
! 2.1 - CAS OU AUCUN ELEMENT N EST SINGULIER
!
    if (nalpha .eq. 1) then
!
        cumm = 0.d0
!
        do inel = 1, nelem
            if (nbr(inel) .eq. 3) then
                if (lqi) then
                    cumm = cumm+(erreur(inel)**(d/(2.d0*ordre+d)))
                else
                    cumm = cumm+(erreur(inel)**(2.d0*d/(2.d0*ordre+d)))
                end if
            end if
        end do
!
        cumm = cumm**(1.d0/(2.d0*ordre))
!
        do inel = 1, nelem
            if (nbr(inel) .eq. 3) then
                if (lqi) then
                    re(inel) = erreur(inel)**(1.d0/(2.d0*ordre+d))
                    re(inel) = 1.d0/(re(inel)*cumm)
                    re(inel) = (prec0**(1.d0/ordre))*re(inel)
                else
                    re(inel) = erreur(inel)**(2.d0/(2.d0*ordre+d))
                    re(inel) = 1.d0/(re(inel)*cumm)
                    re(inel) = (prec0**(1.d0/2.d0*ordre))*re(inel)
                end if
            end if
        end do
!
! 2.2 - CAS OU CERTAINS ELEMENTS SONT SINGULIERS
! 2.2.1 - CALCUL DES COEFFICIENTS AE ET BE POUR SIMPLIFIER EXPRESSION
!
    else
        mu = 0.d0
        do inel = 1, nelem
            if (nbr(inel) .eq. 3) then
                ae(inel) = 2.d0*alpha(inel)/(2.d0*alpha(inel)+d)
                if (lqi) then
                    be(inel) = erreur(inel)**(d/(2.d0*alpha(inel)))
                else
                    be(inel) = erreur(inel)**(d/alpha(inel))
                end if
                be(inel) = d*be(inel)/(2.d0*alpha(inel))
                be(inel) = be(inel)**ae(inel)
                mu = mu+(erreur(inel)**(2.d0*ordre/(2.d0*ordre+d)))
            end if
        end do
!
! 2.2.2 - RECHERCHE DU LAGRANGIEN MU PAR MEHODE DE NEWTON
!         OU DICHOTOMIE
!
        mu = mu/(prec0**2.d0)
        mu = mu**((2.d0*ordre+d)/(2.d0*ordre))
        mu = d*mu/(2.d0*ordre)
        iter = 1
!
70      continue
!
        if (iter .le. 15) then
            if (lqi) then
                fonc = -prec0
            else
                fonc = -(prec0**2.d0)
            end if
            dfonc = 0.d0
            do inel = 1, nelem
                if (nbr(inel) .eq. 3) then
                    fonc = fonc+be(inel)/(mu**ae(inel))
                    dfonc = dfonc-ae(inel)*be(inel)/(mu**(ae(inel)+1.d0) &
                                                     )
                end if
            end do
!
            if (abs(fonc) .le. 1.d-06) goto 60
            if (fonc .ge. 0.d0) then
                mu = mu-fonc/dfonc
                iter = iter+1
                goto 70
            else
                mu = mu/2.d0
                iter = iter+1
                goto 70
            end if
!
60          continue
!
            do inel = 1, nelem
                if (nbr(inel) .eq. 3) then
                    if (lqi) then
                        re(inel) = d/(2.d0*mu*alpha(inel)*erreur(inel))
                    else
                        re(inel) = d/(2.d0*mu*alpha(inel)*(erreur(inel) &
                                                           **2.d0))
                    end if
                    re(inel) = re(inel)**(1.d0/(2.d0*alpha(inel)+d))
                end if
            end do
!
        else
!
            call utmess('F', 'CALCULEL3_99')
!
        end if
!
!
    end if
!
end subroutine
