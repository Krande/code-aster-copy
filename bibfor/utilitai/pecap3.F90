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

subroutine pecap3(chgeoz, tempez, iomega)
!.======================================================================
    implicit none
!
!      PECAP3  -- CALCUL DE LA CONSTANTE DE GAUCHISSEMENT D'UNE POUTRE
!                 DEFINIE PAR SA SECTION MAILLEE EN ELEMENTS
!                 MASSIFS 2D.
!                 CETTE CONSTANTE DE GAUCHISSEMENT EST AUSSI
!                 APPELEE INERTIE DE GAUCHISSEMENT.
!
!          .LE DOMAINE SUR-LEQUEL ON TRAVAILLE REPRESENTE LA
!           SECTION DE LA POUTRE MAILLEE AVEC DES ELEMENTS 2D
!           ISOPARAMETRIQUES THERMIQUES (THERMIQUES CAR ON
!           DOIT RESOUDRE UNE EQUATION DE LAPLACE).
!
!          .LA CONSTANTE DE GAUCHISSEMENT IOMEGA EST DETERMINEE
!           DE LA MANIERE SUIVANTE :
!             EN SE PLACANT DANS LE REPERE PRINCIPAL D'INERTIE
!             AVEC COMME ORIGINE LE CENTRE DE TORSION
!
!             1)SI L'ON ECRIT L'EQUATION D'EQUILIBRE LOCALE :
!                DIV(SIGMA) = 0 SELON L'AXE DE LA POUTRE
!                ON OBTIENT POUR LE PROBLEME DE TORSION L'EQUATION :
!                LAPLACIEN(OMEGA) = 0     DANS LA SECTION
!
!             2)D'AUTRE-PART, LA SECTION ETANT EN EQUILIBRE
!               LA FORCE NORMALE AU CONTOUR EN TOUT POINT DE CE
!               CONTOUR EST NULLE , SOIT (SIGMA).N = 0
!              CE QUI DONNE POUR LE PROBLEME DE TORSION :
!     D(OMEGA)/D(N) = Z*NY-Y*NZ   SUR LE CONTOUR DE LA SECTION
!     NY ET NZ ETANT LES COMPOSANTES DU VECTEUR N NORMAL A CE CONTOUR
!
!             3)ON OBTIENT LA CONDITION DIRICHLET PERMETTANT
!     DE RESOUDRE LE PROBLEME EN ECRIVANT :
!             SOMME/SECTION(OMEGA.DS) = 0
!     (CA VIENT DE L'EQUATION D'EQUILIBRE SELON L'AXE DE LA POUTRE
!      SOIT N = 0 , N ETANT L'EFFORT NORMAL
!      OR N = SOMME/SECTION(SIGMA_XX.DS)
!         N = SOMME/SECTION(E*OMEGA(Y,Z)*THETA_X,XX.DS) )
!
!     ON A ALORS IOMEGA = SOMME_S(OMEGA**2.DS)
!
!     L'OPTION : 'CARA_GAUCHI'   CALCULE :
!       IOMEGA  =  SOMME/SECTION(OMEGA**2.DS)
!
!   ARGUMENT        E/S  TYPE         ROLE
!    CHGEOZ         IN    K*      COORDONNEES DES CONNECTIVITES
!                                 DANS LE REPERE PRINCIPAL D'INERTIE
!    TEMPEZ         IN    K*      RESULTAT DE TYPE EVOL_THER
!                                 REFERENCANT LE CHAMP DE SCALAIRES
!                                 SOLUTION DE L'EQUATION 1
!    IOMEGA         OUT   R       CONSTANTE DE GAUCHISSEMENT
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/mesomm.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/utmess.h"
    character(len=*) :: chgeoz, tempez
    real(kind=8) :: iomega
! -----  VARIABLES LOCALES
    character(len=8) :: lpain(2), lpaout(1)
    character(len=8) :: temper
    character(len=8) :: crit, modele
    character(len=19) :: nume_equa
    character(len=14) :: typres
    character(len=19) :: knum, ligrth
    character(len=24) :: lchin(2), lchout(1), chgeom
    character(len=24) :: chtemp
    real(kind=8) :: work(9)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! ---- INITIALISATIONS
!      ---------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iret, nbordr
    real(kind=8) :: prec, zero
!-----------------------------------------------------------------------
    zero = 0.0d0
    prec = 1.0d-3
    chgeom = chgeoz
    temper = tempez
    knum = '&&PECAP3.NUME_ORD_1'
    crit = 'RELATIF'
!
    do i = 1, 9
        work(i) = zero
    end do
!
! --- ON VERIFIE QUE LE RESULTAT EST DE TYPE EVOL_THER :
!     ------------------------------------------------
    call dismoi('TYPE_RESU', temper, 'RESULTAT', repk=typres)
    if (typres .ne. 'EVOL_THER') then
        call utmess('F', 'UTILITAI3_57')
    end if
!
! --- RECUPERATION DU NOMBRE D'ORDRES DU RESULTAT :
!     -------------------------------------------
    call rsutnu(temper, ' ', 0, knum, nbordr, &
                prec, crit, iret)
    if (nbordr .ne. 1) then
        call utmess('F', 'UTILITAI3_58', sk=temper)
    end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      CALL UTIMSD(IFM,2,.FALSE.,.TRUE.,TEMPER,1,' ')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! --- RECUPERATION DU CHAMP DE TEMPERATURES DU RESULTAT :
!     -------------------------------------------------
    call rsexch('F', temper, 'TEMP', 1, chtemp, &
                iret)
!
! --- RECUPERATION DU NUME_DDL ASSOCIE AU CHAMP DE TEMPERATURES :
!     ---------------------------------------------------------
    call dismoi('NUME_EQUA', chtemp, 'CHAM_NO', repk=nume_equa)
!
! --- RECUPERATION DU MODELE ASSOCIE AU NUME_DDL  :
!     ------------------------------------------
    call dismoi('NOM_MODELE', nume_equa, 'NUME_EQUA', repk=modele)
!
! --- RECUPERATION DU LIGREL DU MODELE  :
!     --------------------------------
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrth)
!
! --- CALCUL POUR CHAQUE ELEMENT DE LA SECTION DE L'INTEGRALE DU
! --- CHAMP DE SCALAIRES SOLUTION DE L'EQUATION DE LAPLACE AU CARRE :
!     -------------------------------------------------------------
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PTEMPER'
    lchin(2) = chtemp
    lpaout(1) = 'PCASECT'
    lchout(1) = '&&PECAP3.INTEG1'
!
    call calcul('S', 'CARA_GAUCHI', ligrth, 2, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
!
! --- SOMMATION DES INTEGRALES PRECEDENTES SUR LA SECTION DE LA POUTRE
! --- (I.E. CALCUL DE SOMME_SECTION_POUTRE(OMEGA**2.DS)) :
!     --------------------------------------------------
    call mesomm(lchout(1), 9, vr=work)
    iomega = work(1)
!
    call detrsd('CHAMP_GD', '&&PECAP3.INTEG')
!.============================ FIN DE LA ROUTINE ======================
end subroutine
