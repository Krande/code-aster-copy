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
! ==================================================================================================
!
! Module for stabilization of HEXA cells for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Stabilization_Hexa_module
! ==================================================================================================
    use SolidShell_type
    use SolidShell_Utilities_module
    use SolidShell_Kinematic_Hexa_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public  :: compStabMatrMateHexa, compStabModulusHexa, compStabForcHexa, compStabMatrGeomHexa, &
               compStabSigmHexa
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/SolidShell_type.h"
#include "asterfort/assert.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! Compute effective shear modulus for HEXA
!
! In  Ueff             : scalar for stabilization
! In  sigm             : stress tensor
! In  epsi             : strain tensor
! In  dsidep           : material jacobian matrix (dSigm / dEpsi)
! In  kineHexa         : kinematic quantities for HEXA cell
! IO  stabHexa         : stabilization for HEXA cell
! Out Ueff             : effective shear modulus
!
! --------------------------------------------------------------------------------------------------
    subroutine compStabModulusHexa(sigm, epsi, dsidep, Ueff)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        real(kind=8), intent(in) :: sigm(SSH_SIZE_TENS), epsi(SSH_SIZE_TENS)
        real(kind=8), intent(in) :: dsidep(SSH_SIZE_TENS, SSH_SIZE_TENS)
        real(kind=8), intent(out) :: Ueff
! - Local
        real(kind=8) :: sigmVol, epsiVol
        real(kind=8) :: sigmDev(SSH_SIZE_TENS), epsiDev(SSH_SIZE_TENS)
        real(kind=8) :: aux1, aux2
!   ------------------------------------------------------------------------------------------------
!
        Ueff = 0.d0

! - Volumic part
        sigmVol = (sigm(1)+sigm(2)+sigm(3))/3.d0
        epsiVol = (epsi(1)+epsi(2)+epsi(3))/3.d0

! - Deviatoric part

        sigmDev = sigm-sigmVol*tensorIden
        epsiDev = epsi-epsiVol*tensorIden

! - Compute
        aux1 = sigmDev(1)*sigmDev(1)+ &
               sigmDev(2)*sigmDev(2)+ &
               sigmDev(3)*sigmDev(3)+ &
               2.d0*(sigmDev(4)*sigmDev(4)+sigmDev(5)*sigmDev(5)+sigmDev(6)*sigmDev(6))
        aux2 = epsiDev(1)*epsiDev(1)+ &
               epsiDev(2)*epsiDev(2)+ &
               epsiDev(3)*epsiDev(3)+ &
               (epsiDev(4)*epsiDev(4)+epsiDev(5)*epsiDev(5)+epsiDev(6)*epsiDev(6))/2.d0

        if (abs(aux2) .le. hexaStabMini) then
            Ueff = (dsidep(4, 4)+dsidep(5, 5)+dsidep(6, 6))/3.d0
        else
            Ueff = 0.5d0*sqrt(aux1/aux2)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compStabMatrMateHexa
!
! Compute stabilization matrix for material part
! - Compute the deformation gradient B in cartesian base - Deviatoric part
! - D'apres Hughes ces termes sont deja deviatoric dans leur forme initiale.
!
! In  geomHexa         : geometric properties for HEXA cell
! In  kineHexa         : kinematic quantities for HEXA cell
! In  Ueff             : scalar for stabilization
! IO  stabHexa         : stabilization for HEXA cell
!
! --------------------------------------------------------------------------------------------------
    subroutine compStabMatrMateHexa(geomHexa, kineHexa, Ueff, stabHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_GEOM_HEXA), intent(in)    :: geomHexa
        type(SSH_KINE_HEXA), intent(in)    :: kineHexa
        real(kind=8), intent(in)           :: Ueff
        type(SSH_STAB_HEXA), intent(inout) :: stabHexa
! - Local
        real(kind=8) :: coef1, coef2
!   ------------------------------------------------------------------------------------------------
!
        stabHexa%matrStabMate = 0.d0
        stabHexa%SXI = 0.d0
        stabHexa%SETA = 0.d0
        stabHexa%SETAZETA = 0.d0
        stabHexa%SXIZETA = 0.d0

! - Coefficient for stabilization
        coef1 = geomHexa%cellGeom%detJac0*hexaStabCoef1
        coef2 = geomHexa%cellGeom%detJac0*hexaStabCoef2

! - Compute products
        call prodBTDB(Ueff*hexaStabDMatr, SSH_SIZE_TENS, SSH_NBDOFG_HEXA, kineHexa%BCartX, &
                      stabHexa%SXI)
        call prodBTDB(Ueff*hexaStabDMatr, SSH_SIZE_TENS, SSH_NBDOFG_HEXA, kineHexa%BCartY, &
                      stabHexa%SETA)
        call prodBTDB(Ueff*hexaStabDMatr, SSH_SIZE_TENS, SSH_NBDOFG_HEXA, kineHexa%BCartYZ, &
                      stabHexa%SETAZETA)
        call prodBTDB(Ueff*hexaStabDMatr, SSH_SIZE_TENS, SSH_NBDOFG_HEXA, kineHexa%BCartXZ, &
                      stabHexa%SXIZETA)

! - Compute stabilization matrix for material part
        stabHexa%matrStabMate = coef1*(stabHexa%SXI+stabHexa%SETA)+ &
                                coef2*(stabHexa%SXIZETA+stabHexa%SETAZETA)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compStabSigmHexa
!
! Compute stabilized stresses for HEXA cell
!
! In  geomHexa         : geometric properties for HEXA cell
! In  kineHexa         : kinematic quantities for HEXA cell
! In  Ueff             : scalar for stabilization
! In  disp             : displacements of element
! IO  stabHexa         : stabilization for HEXA cell
! In  epsg             : Green-Lagrange strains
!
! --------------------------------------------------------------------------------------------------
    subroutine compStabSigmHexa(geomHexa, kineHexa, Ueff, disp, stabHexa, epsg_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_GEOM_HEXA), intent(in)           :: geomHexa
        type(SSH_KINE_HEXA), intent(in)           :: kineHexa
        real(kind=8), intent(in)                  :: Ueff, disp(SSH_NBDOF_HEXA)
        type(SSH_STAB_HEXA), intent(inout)        :: stabHexa
        type(SSH_EPSG_HEXA), optional, intent(in) :: epsg_
! - Local
    real(kind=8) :: EH1(SSH_SIZE_TENS), EH2(SSH_SIZE_TENS), EH23(SSH_SIZE_TENS), EH13(SSH_SIZE_TENS)
!   ------------------------------------------------------------------------------------------------
!
        stabHexa%sigmStabXI = 0.d0
        stabHexa%sigmStabETA = 0.d0
        stabHexa%sigmStabXIZETA = 0.d0
        stabHexa%sigmStabETAZETA = 0.d0

! - Compute the deformation in cartesian frame
        if (kineHexa%lLarge) then
            EH1 = matmul(geomHexa%T0, epsg_%ECovaXI)+ &
                  matmul(geomHexa%TXI, epsg_%ECova0)
            EH2 = matmul(geomHexa%T0, epsg_%ECovaETA)+ &
                  matmul(geomHexa%TETA, epsg_%ECova0)
            EH23 = matmul(geomHexa%T0, epsg_%ECovaETAZETA)+ &
                   matmul(geomHexa%TETA, epsg_%ECovaZETA)+ &
                   matmul(geomHexa%TZETA, epsg_%ECovaETA)
            EH13 = matmul(geomHexa%T0, epsg_%ECovaXIZETA)+ &
                   matmul(geomHexa%TXI, epsg_%ECovaZETA)+ &
                   matmul(geomHexa%TZETA, epsg_%ECovaXI)
        else
            EH1 = matmul(kineHexa%BCartX, disp(1:SSH_NBDOFG_HEXA))
            EH2 = matmul(kineHexa%BCartY, disp(1:SSH_NBDOFG_HEXA))
            EH23 = matmul(kineHexa%BCartYZ, disp(1:SSH_NBDOFG_HEXA))
            EH13 = matmul(kineHexa%BCartXZ, disp(1:SSH_NBDOFG_HEXA))
        end if

! - Compute stabilisation stresses
        stabHexa%sigmStabXI = matmul(Ueff*hexaStabDMatr, EH1)
        stabHexa%sigmStabETA = matmul(Ueff*hexaStabDMatr, EH2)
        stabHexa%sigmStabXIZETA = matmul(Ueff*hexaStabDMatr, EH13)
        stabHexa%sigmStabETAZETA = matmul(Ueff*hexaStabDMatr, EH23)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compStabForcHexa
!
! Compute stabilized force for HEXA cell
!
! In  kineHexa         : kinematic quantities for HEXA cell
! IO  stabHexa         : stabilization for HEXA cell
!
! --------------------------------------------------------------------------------------------------
    subroutine compStabForcHexa(kineHexa, stabHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_KINE_HEXA), intent(in)           :: kineHexa
        type(SSH_STAB_HEXA), intent(inout)        :: stabHexa
!   ------------------------------------------------------------------------------------------------
!
        stabHexa%forcStab = 0.d0

! - Compute stabilization for internal forcess
        stabHexa%forcStab = &
            hexaStabCoef1*(matmul(transpose(kineHexa%BCartX), stabHexa%sigmStabXI)+ &
                           matmul(transpose(kineHexa%BCartY), stabHexa%sigmStabETA))+ &
            hexaStabCoef2*(matmul(transpose(kineHexa%BCartYZ), stabHexa%sigmStabETAZETA)+ &
                           matmul(transpose(kineHexa%BCartXZ), stabHexa%sigmStabXIZETA))
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compStabMatrGeomHexa
!
! Compute stabilization matrix for geometric part
!
! In  geomHexa         : geometric properties for HEXA cell
! In  kineHexa         : kinematic quantities for HEXA cell
! IO  stabHexa         : stabilization for HEXA cell
!
! --------------------------------------------------------------------------------------------------
    subroutine compStabMatrGeomHexa(geomHexa, kineHexa, stabHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_GEOM_HEXA), intent(in)    :: geomHexa
        type(SSH_KINE_HEXA), intent(in)    :: kineHexa
        type(SSH_STAB_HEXA), intent(inout) :: stabHexa

! - Local
        integer(kind=8), parameter :: nbNodeGeom = SSH_NBNODEG_HEXA
        integer(kind=8) :: iNodeGeom, jNodeGeom
        real(kind=8) :: const, coef1, coef2
        real(kind=8) :: GCova0(SSH_SIZE_TENS), GCovaZETA(SSH_SIZE_TENS)
        real(kind=8) :: GCartXI(SSH_SIZE_TENS), GCartETA(SSH_SIZE_TENS)
        real(kind=8) :: GCartETAZETA(SSH_SIZE_TENS), GCartXIZETA(SSH_SIZE_TENS)
        real(kind=8) :: GCovaXI(SSH_SIZE_TENS), GCovaETA(SSH_SIZE_TENS)
        real(kind=8) :: GCovaETAZETA(SSH_SIZE_TENS), GCovaXIZETA(SSH_SIZE_TENS)
!   ------------------------------------------------------------------------------------------------
!
        stabHexa%matrStabGeom = 0.d0
        ASSERT(kineHexa%lLarge)

! - Coefficient for stabilization
        coef1 = geomHexa%cellGeom%detJac0*hexaStabCoef1
        coef2 = geomHexa%cellGeom%detJac0*hexaStabCoef2

! - For "standard" nodes
        do iNodeGeom = 1, nbNodeGeom
            do jNodeGeom = 1, nbNodeGeom

! --------- Compute gradients (covariant) for geometric matrix
                call compGCovaMatrHexa(iNodeGeom, jNodeGeom, &
                                       GCova0, GCovaZETA)
                call compGCovaSMatrHexa(iNodeGeom, jNodeGeom, &
                                        GCovaXI, GCovaETA, &
                                        GCovaETAZETA, GCovaXIZETA)

! --------- Project gradients in cartesian frame
                GCartXi = matmul(geomHexa%T0, GCovaXI)+ &
                          matmul(geomHexa%TXI, GCova0)
                GCartETA = matmul(geomHexa%T0, GCovaETA)+ &
                           matmul(geomHexa%TETA, GCova0)
                GCartETAZETA = matmul(geomHexa%T0, GCovaETAZETA)+ &
                               matmul(geomHexa%TETA, GCovaZETA)+ &
                               matmul(geomHexa%TZETA, GCovaETA)
                GCartXIZETA = matmul(geomHexa%T0, GCovaXIZETA)+ &
                              matmul(geomHexa%TXI, GCovaZETA)+ &
                              matmul(geomHexa%TZETA, GCovaXI)

! --------- Compute matrix
                const = coef1*sum(stabHexa%sigmStabXI*GCartXi)+ &
                        coef1*sum(stabHexa%sigmStabETA*GCartETA)+ &
                        coef2*(sum(stabHexa%sigmStabETAZETA*GCartETAZETA)+ &
                               sum(stabHexa%sigmStabXIZETA*GCartXIZETA))
                stabHexa%matrStabGeom(3*(iNodeGeom-1)+1:3*(iNodeGeom-1)+3, &
                                      3*(jNodeGeom-1)+1:3*(jNodeGeom-1)+3) = const*matr3Iden
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module SolidShell_Stabilization_Hexa_module
