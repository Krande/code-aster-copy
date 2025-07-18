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

subroutine xpraju(noma, fiss, cnslt, cnsvt, cnsvn, &
                  deltat, vmax)
!
    implicit none
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/dismoi.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=8) :: noma, fiss
    character(len=19) :: cnsvt, cnsvn, cnslt
    real(kind=8) :: deltat, vmax
!
! person_in_charge: patrick.massin at edf.fr
!     ------------------------------------------------------------------
!
!       XPRAJU   : X-FEM PROPAGATION : AJUSTEMENT DE VN
!       ------     -     --            ---
!    AJUSTEMENT DU CHAMP DE VITESSE VN :
!          SI  LT <=0 , VN AJUSTEE = 0
!          SINON, VN AJUSTEE = (VN*LST)/(VT*DELTAT)
!
!    ENTREE
!        NOMA    : NOM DU CONCEPT MAILLAGE
!        FISS    : NOM DU CONCEPT FISSURE X-FEM DE LA FISSURE A PROPAGER
!        CNSLT   : CHAM_NO_S LEVEL SET TANGENTIELLE
!        CNSVT   : CHAM_NO_S VITESSE TANGENTIELLE DE PROPAGATION
!        CNSVN   : CHAM_NO_S VITESSE NORMALE DE PROPAGATION
!        DELTAT  : TEMPS TOTAL DE PROPAGATION
!
!    SORTIE
!        CNSVN   : CHAM_NO_S VITESSE NORMALE DE PROPAGATION AJUSTEE
!        VMAX    : VALEUR MAXIMALE DES COMPOSANTES DE VITESSE
!
!     ------------------------------------------------------------------
!
!
    integer(kind=8) :: i, nbno, ifm, niv, cptzo, cptaju
    integer(kind=8) :: jlisno
    real(kind=8) :: modzon, dmin
    real(kind=8), pointer :: ltno(:) => null()
    real(kind=8), pointer :: vnno(:) => null()
    real(kind=8), pointer :: vtno(:) => null()
!
!-----------------------------------------------------------------------
!     DEBUT
!-----------------------------------------------------------------------
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
!
!  RECUPERATION DU NOMBRE DE NOEUDS DU MAILLAGE
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbno)
!
!   RECUPERATION DES ADRESSES DES CHAMPS DE VITESSE AUX NOEUDS
    call jeveuo(cnsvt//'.CNSV', 'L', vr=vtno)
    call jeveuo(cnsvn//'.CNSV', 'E', vr=vnno)
!
!   RECUPERATION DE L'ADRESSE DES VALEURS DE LST
    call jeveuo(cnslt//'.CNSV', 'L', vr=ltno)
!
!   RETRIEVE THE LIST OF THE NODES THAT MUST TO BE USED IN THE
!   CALCULUS
    call jeveuo(fiss//'.PRO.NOEUD_TORE', 'L', jlisno)
!
    cptzo = 0
    cptaju = 0
    dmin = r8miem()
    vmax = 0.d0
!
!   BOUCLE SUR TOUS LES NOEUDS DU MAILLAGE
    do i = 1, nbno
!
        if (zl(jlisno-1+i)) then
!
            if (ltno(i) .le. dmin) then
!
!             THE NODE (OR ITS PROJECTION) IS ON THE EXISTING CRACK
!             SURFACE. ITS NORMAL SPEED MUST BE SET TO ZERO.
                vnno(i) = 0
!
!             CALCULATE THE MAXIMUM VALUE OF THE SPEED COMPONENTS
                if (abs(vtno(i)) .gt. vmax) vmax = abs(vtno(i))
!
                cptzo = cptzo+1
!
            else
!
!             THE NODE (OR ITS PROJECTION) IS AHEAD OF THE CRACK TIP.
!             ITS NORMAL SPEED MUST BE RECALCULATED USING A LINEAR
!             EXTRAPOLATION.
                modzon = vtno(i)*deltat
                vnno(i) = vnno(i)*ltno(i)/modzon
!
!             CALCULATE THE MAXIMUM VALUE OF THE SPEED COMPONENTS
                if (abs(vtno(i)) .gt. vmax) vmax = abs(vtno(i))
                if (abs(vnno(i)) .gt. vmax) vmax = abs(vnno(i))
!
                cptaju = cptaju+1
!
            end if
!
        end if
!
    end do
!
    if (niv .ge. 1) then
        write (ifm, *) '   NOMBRE DE NOEUDS DONT VN EST ANNULEE :', cptzo
        write (ifm, *) '   NOMBRE DE NOEUDS DONT VN EST AJUSTEE :', cptaju
    end if
!
!-----------------------------------------------------------------------
!     FIN
!-----------------------------------------------------------------------
    call jedema()
end subroutine
