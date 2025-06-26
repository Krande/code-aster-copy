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
subroutine lcafyd(compor, materd, materf, nbcomm, cpmono, &
                  nmat, mod, nvi, vind, &
                  sigd, nr, yd)
    implicit none
!
!
!     CHOIX DES VALEURS DE VIND A AFFECTER A YD
!     CAS PARTICULIER DU  MONOCRISTAL  :
!     ON GARDE 1 VARIABLE INTERNE PAR SYSTEME DE GLISSEMENT SUR 3
!     ----------------------------------------------------------------
!     IN
!          COMP   :  NOM MODELE DE COMPORTEMENT
!          MATERD :  COEF MATERIAU A T
!          MATERF :  COEF MATERIAU A T+DT
!          NBCOMM :  INDICES DES COEF MATERIAU
!          NMAT   :  DIMENSION MATER
!          COMP   :  TYPE DE MODELISATION
!          NVI    :  NOMBRE DE VARIABLES INTERNES
!          VIND   :  VARIABLES INTERNES A T
!          VINF   :  VARIABLES INTERNES A T+DT (BASE SUR PRED_ELAS)
!          NR     :  DIMENSION VECTEUR INCOONUES
!          SIGD   :  ETAT DE CONTRAINTES A T
!     OUT  YD     :  VECTEUR INITIAL
!
!     ----------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/assert.h"
#include "asterfort/lcgrla.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcopli.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
    integer(kind=8) :: ndt, nvi, nmat, ndi, ns, i, nbcomm(nmat, 3), nr
    real(kind=8) :: yd(*), materd(nmat, 2), materf(nmat, 2), vind(*)
    real(kind=8) :: id(3, 3), hookf(6, 6), dkooh(6, 6), epsegl(6), fe(3, 3)
    real(kind=8) :: dtot, sigd(6)
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    character(len=16) :: rela_comp
    character(len=24) :: cpmono(5*nmat+1), necoul
    character(len=8) :: mod
    common/tdim/ndt, ndi
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    blas_int :: b_incx, b_incy, b_n
    common/polycr/irr, decirr, nbsyst, decal, gdef
    data id/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!     ----------------------------------------------------------------
!
!     INITIALISATION DE YD EN IMPLICITE
    rela_comp = compor(RELA_NAME)
!
!
!     AFFECTATION DE YD = ( SIGD , VIND , (EPSD(3)) )
!
    yd(1:ndt) = sigd(1:ndt)
!
    if (rela_comp .eq. 'MONOCRISTAL') then
! ATTENTION !         NS=(NVI-8)/3
        ns = nr-ndt
        irr = 0
        decirr = 0
        if (materf(nbcomm(1, 1), 2) .ge. 4) then
!           KOCKS-RAUCH ET DD_CFC : VARIABLE PRINCIPALE=DENSITE DISLOC
!           UNE SEULE FAMILLE
            ASSERT(nbcomm(nmat, 2) .eq. 1)
            do i = 1, ns
                yd(ndt+i) = vind(6+3*(i-1)+1)
            end do
            necoul = cpmono(3)
            if (necoul .eq. 'MONO_DD_CC_IRRA') then
                irr = 1
                decirr = 6+3*ns
            end if
            if (necoul .eq. 'MONO_DD_CFC_IRRA') then
                irr = 1
                decirr = 6+3*ns
            end if
        else
!           AUTRES COMPORTEMENTS MONOCRISTALLINS
            do i = 1, ns
                yd(ndt+i) = vind(6+3*(i-1)+2)
            end do
        end if
!
!
        if (gdef .eq. 1) then
! les 9 variables internes  de 6+3*ns+1 à 6+3*ns+9
! REPRESENTENT FE - ID
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vind(nvi-3-18+10), b_incx, fe, b_incy)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, +1.d0, id, b_incx, fe, &
                       b_incy)
            call lcgrla(fe, epsegl)
            if (materf(nmat, 2) .eq. 0) then
                call lcopli('ISOTROPE', mod, materf(1, 1), hookf)
            else if (materf(nmat, 2) .eq. 1) then
                call lcopli('ORTHOTRO', mod, materf(1, 1), hookf)
            end if
! Y contient H*(FeT.Fe-Id)/2, ce ne sont pas exactement les PK2
! Y contient ensuite les ns alpha_s ou gamma_s suivant la rela_comp
            yd(1:ndt) = matmul(hookf(1:ndt, 1:ndt), epsegl(1:ndt))
        end if
!
!
!
    else if (rela_comp .eq. 'IRRAD3M') then
!        CORRESPONDANCE ENTRE LES VARIABLES INTERNES ET LES EQUATIONS
!        DU SYSTEME DIFFERENTIEL
!        DEFORMATION PLASTIQUE CUMULEE
        yd(ndt+1) = vind(1)
!        FONCTION SEUIL DE FLUAGE
        yd(ndt+2) = vind(2)
!        DEFORMATION EQUIVALENTE DE FLUAGE
        yd(ndt+3) = vind(3)
!        DEFORMATION DE GONFLEMENT
        yd(ndt+4) = vind(4)
!
    else if (rela_comp .eq. 'LETK') then
! --- INITIALISATION A ZERO DU MULTIPLICATEUR PLASTIQUE
        yd(ndt+1) = 0.d0
! --- INITIALISATION A XIP
        yd(ndt+2) = vind(1)
! --- INITIALISATION A XIVP
        yd(ndt+3) = vind(3)
!
    else if (rela_comp .eq. 'LKR') then
! --- INITIALISATION A ZERO DU MULTIPLICATEUR PLASTIQUE
        yd(ndt+1) = 0.d0
! --- INITIALISATION A XIP
        yd(ndt+2) = vind(1)
! --- INITIALISATION A XIVP
        yd(ndt+3) = vind(3)
!
    else if (rela_comp .eq. 'HAYHURST') then
        call lcopil('ISOTROPE', mod, materd(1, 1), dkooh)
!        DEFORMATION ELASTIQUE INSTANT PRECEDENT
        yd(1:ndt) = matmul(dkooh(1:ndt, 1:ndt), sigd(1:ndt))
        dtot = 1.d0/(1.d0-vind(11))
        yd(1:ndt) = dtot*yd(1:ndt)
!
!        CORRESPONDANCE ENTRE LES VARIABLES INTERNES ET LES EQUATIONS
!        DU SYSTEME DIFFERENTIEL
!        DEFORMATION PLASTIQUE CUMULEE
        yd(ndt+1) = vind(7)
!        H1
        yd(ndt+2) = vind(8)
!        H2
        yd(ndt+3) = vind(9)
!        D
        yd(ndt+4) = vind(11)
!
    else
!     CAS GENERAL :
!        TOUTES LES VARIABLES INTERNES SONT RECOPIES
!        LA DERNIERE C'EST TOUJOURS L'INDICATEUR PLASTIQUE
        yd(ndt+1:ndt+nvi-1) = vind(1:nvi-1)
    end if
!
end subroutine
