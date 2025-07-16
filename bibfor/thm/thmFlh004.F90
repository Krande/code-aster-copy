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
! person_in_charge: sylvie.granet at edf.fr
! aslint: disable=W1504
!
subroutine thmFlh004(ds_thm, lMatr, lSigm, ndim, j_mater, &
                     dimdef, dimcon, &
                     addep1, addep2, adcp11, adcp12, adcp21, &
                     addeme, addete, &
                     t, p2, pvp, &
                     grat, grap1, grap2, &
                     rho11, h11, h12, &
                     satur, dsatur, gravity, tperm, &
                     congep, dsde)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/hmderp2.h"
#include "asterfort/thmEvalPermLiquGaz.h"
#include "asterfort/thmEvalFickSteam.h"
!
    type(THM_DS), intent(in) :: ds_thm
    aster_logical, intent(in) :: lMatr, lSigm
    integer(kind=8), intent(in) :: j_mater
    integer(kind=8), intent(in) :: ndim, dimdef, dimcon
    integer(kind=8), intent(in) :: addeme, addep1, addep2, addete, adcp11, adcp12, adcp21
    real(kind=8), intent(in) :: rho11, satur, dsatur
    real(kind=8), intent(in) :: grat(3), grap1(3), grap2(3)
    real(kind=8), intent(in) :: p2, pvp, t
    real(kind=8), intent(in) :: gravity(3), tperm(ndim, ndim)
    real(kind=8), intent(in) :: h11, h12
    real(kind=8), intent(inout) :: congep(1:dimcon)
    real(kind=8), intent(inout) :: dsde(1:dimcon, 1:dimdef)
!
! ----------------------------------------------------------------------------------------------
!
! THM
!
! Compute flux and stress for hydraulic - 'LIQU_VAPE_GAZ'
!
! ------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  ndim             : dimension of space (2 or 3)
! In  j_mater          : coded material address
! In  dimdef           : dimension of generalized strains vector
! In  dimcon           : dimension of generalized stresses vector
! In  addeme           : adress of mechanic dof in vector of generalized strains
! In  addete           : adress of thermic dof in vector of generalized strains
! In  addep1           : adress of first hydraulic dof in vector of generalized strains
! In  addep2           : adress of second hydraulic dof in vector of generalized strains
! In  adcp11           : adress of first hydraulic/first component dof in vector of gene. stresses
! In  adcp12           : adress of first hydraulic/second component dof in vector of gene. stresses
! In  adcp21           : adress of second hydraulic/first component dof in vector of gene. stresses
! In  t                : temperature - At end of current step
! In  p2               : gaz pressure - At end of current step
! In  pvp              : steam pressure
! In  grat             : gradient of temperature
! In  grap1            : gradient of capillary pressure
! In  grap2            : gradient of gaz pressure
! In  rho11            : volumic mass for liquid
! In  h11              : enthalpy of liquid
! In  h12              : enthalpy of steam
! In  satur            : saturation
! In  dsatur           : derivative of saturation (/pc)
! In  gravity          : gravity
! In  tperm            : permeability tensor
! IO  congep           : generalized stresses - At end of current step
! IO  dsde             : derivative matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k
    real(kind=8) :: rgaz
    real(kind=8) :: permli, dperml
    real(kind=8) :: permgz, dperms, dpermp
    real(kind=8) :: dfickt, dfickg, krel2, fick
    real(kind=8) :: krel1, dkrel1, rho12, rho21, masart, masvrt
    real(kind=8) :: cliq, alpliq
    real(kind=8) :: viscl, dviscl, viscg, dviscg
    real(kind=8) :: mamola, mamolv, rhog
    real(kind=8) :: lambd1(5), lambd2(5), fv(5)
    real(kind=8) :: cvp, gpv(3), gc(3)
    real(kind=8) :: dp11p1, dp11p2, dp11t
    real(kind=8) :: dp21p1, dp21p2, dp21t
    real(kind=8) :: dp12p1, dp12p2, dp12t
    real(kind=8) :: dr11p1, dr11p2, dr11t
    real(kind=8) :: dr12p1, dr12p2, dr12t
    real(kind=8) :: dr21p1, dr21p2, dr21t
    real(kind=8) :: dauxp1, dauxp2, dauxt
    real(kind=8) :: dgpvp1(3), dgpvp2(3), dgpvt(3), gr12(3), gr21(3)
    real(kind=8) :: dgcvp1(3), dgcvp2(3), dgcvt(3)
    real(kind=8) :: dgpgt(2)
    real(kind=8) :: dgcgp1(2), dgcgp2(2), dgcgt(2)
    real(kind=8) :: dgpgp1(2), dgpgp2(2)
    real(kind=8) :: dp1pp1, dp2pp1, dtpp1, dp1pp2, dp2pp2
    real(kind=8) :: dtpp2, dp1pt, dp2pt, dtpt
! derivative of grad_rhovp :
    real(kind=8) :: dgrvp1(3), dgrvp2(3), dgrvt(3), dgrvgp1, dgrvgp2, dgrvgt
! derivative of grad_rhoas :
    real(kind=8) :: dgrasp1(3), dgrasp2(3), dgrast(3), dgrasgp1, dgrasgp2, dgrasgt
!
! ------------------------------------------------------------------------------------------------
!
    lambd1 = 0.d0
    lambd2 = 0.d0
    fv = 0.d0
    cvp = 0.d0
    gpv = 0.d0
    gr12 = 0.d0
    gr21 = 0.d0
    gc = 0.d0
    dr11p1 = 0.d0
    dr11p2 = 0.d0
    dr11t = 0.d0
    dr12p1 = 0.d0
    dr12p2 = 0.d0
    dr12t = 0.d0
    dr21p1 = 0.d0
    dr21p2 = 0.d0
    dr21t = 0.d0
    dauxp1 = 0.d0
    dauxp2 = 0.d0
    dauxt = 0.d0
    dgpvp1 = 0.d0
    dgpvp2 = 0.d0
    dgpvt = 0.d0
    dgcvp1 = 0.d0
    dgcvp2 = 0.d0
    dgcvt = 0.d0
    dgpgt = 0.d0
    dgcgp1 = 0.d0
    dgcgp2 = 0.d0
    dgcgt = 0.d0
    dgpgp1 = 0.d0
    dgpgp2 = 0.d0
!
    dgrvp1 = 0.d0
    dgrvp2 = 0.d0
    dgrvt = 0.d0
    dgrvgp1 = 0.d0
    dgrvgp2 = 0.d0
    dgrvgt = 0.d0
    dgrvp1 = 0.d0
    dgrvp2 = 0.d0
    dgrvt = 0.d0
    dgrvgp1 = 0.d0
    dgrvgp2 = 0.d0
    dgrvgt = 0.d0
!
    dgrasp1 = 0.d0
    dgrasp2 = 0.d0
    dgrast = 0.d0
    dgrasgp1 = 0.d0
    dgrasgp2 = 0.d0
    dgrasgt = 0.d0
!
! - Evaluate permeability for liquid and gaz
!
    call thmEvalPermLiquGaz(ds_thm, &
                            j_mater, satur, p2, t, &
                            permli, dperml, &
                            permgz, dperms, dpermp)
!
! - Evaluate Fick coefficients for steam in gaz
!
    call thmEvalFickSteam(j_mater, &
                          satur, p2, pvp, t, &
                          fick, dfickt, dfickg)
!
! - Get parameters
!
    krel1 = permli
    dkrel1 = dperml*dsatur
    krel2 = permgz
    rgaz = ds_thm%ds_material%solid%r_gaz
    cliq = ds_thm%ds_material%liquid%unsurk
    alpliq = ds_thm%ds_material%liquid%alpha
    viscl = ds_thm%ds_material%liquid%visc
    dviscl = ds_thm%ds_material%liquid%dvisc_dtemp
    mamolv = ds_thm%ds_material%steam%mass_mol
    viscg = ds_thm%ds_material%steam%visc
    dviscg = ds_thm%ds_material%steam%dvisc_dtemp
    mamola = ds_thm%ds_material%gaz%mass_mol
    rho12 = mamolv*pvp/rgaz/t
    rho21 = mamola*(p2-pvp)/rgaz/t
    rhog = rho12+rho21
    masart = mamola/rgaz/t
    masvrt = mamolv/rgaz/t
!
!  Concentration massique
    cvp = rho12/rhog
!
! - Fick
!
    fv(1) = fick
    fv(2) = 0.d0
    fv(3) = 0.d0
    fv(4) = dfickg
    fv(5) = dfickt
!
! - Thermic conductivity
!
    lambd2(1) = krel2/viscg
    lambd2(2) = 0.d0
    lambd2(3) = dperms*dsatur/viscg
    lambd2(4) = dpermp/viscg
    lambd2(5) = -krel2/viscg/viscg*dviscg
    lambd1(1) = krel1/viscl
    lambd1(2) = 0.d0
    lambd1(3) = dkrel1/viscl
    lambd1(4) = 0.d0
    lambd1(5) = -krel1/viscl/viscl*dviscl
!

! - Compute some derivatives for LIQU_VAPE_GAZ
!
    call hmderp2(ds_thm, t, pvp, &
                 rho11, rho12, h11, h12, &
                 dp11p1, dp11p2, dp11t, &
                 dp12p1, dp12p2, dp12t, &
                 dp21p1, dp21p2, dp21t, &
                 dp1pp1, dp2pp1, dtpp1, &
                 dp1pp2, dp2pp2, dtpp2, &
                 dp1pt, dp2pt, dtpt)
!
! - Pressure gradient (Eq. 5.5.1-7)
!
    do i = 1, ndim
! gradients de pression de vapeur
        gpv(i) = dp12p2*grap2(i)+dp12p1*grap1(i)
        if (ds_thm%ds_elem%l_dof_ther) then
            gpv(i) = gpv(i)+dp12t*grat(i)
        end if
!gradient of  density
        gr12(i) = masvrt*gpv(i)
        gr21(i) = masart*(grap2(i)-gpv(i))
        if (ds_thm%ds_elem%l_dof_ther) then
            gr12(i) = gr12(i)-rho12/t*grat(i)
            gr21(i) = gr21(i)-rho21/t*grat(i)
        end if
        gc(i) = gr12(i)/rhog-rho12/rhog/rhog*(gr12(i)+gr21(i))
    end do
!
! - Volumic mass - Derivative
!
    if (lMatr) then
        dr11p1 = rho11*dp11p1*cliq
        dr11p2 = rho11*dp11p2*cliq
        dr12p1 = rho12/pvp*dp12p1
        dr12p2 = rho12/pvp*dp12p2
        dr21p1 = masart*dp21p1
        dr21p2 = masart*dp21p2
        if (ds_thm%ds_elem%l_dof_ther) then
            dr11t = -3.d0*alpliq*rho11
            dr12t = rho12*(dp12t/pvp-1.d0/t)
            dr21t = masart*dp12t-rho21/t
        end if
! ----- GRADPVP and GRADCVP - Derivative
        do i = 1, ndim
! --------- GRADPVP - Derivative
            dgpvp1(i) = dp1pp2*grap2(i)+dp1pp1*grap1(i)
            dgpvp2(i) = dp2pp2*grap2(i)+dp2pp1*grap1(i)
            if (ds_thm%ds_elem%l_dof_ther) then
                dgpvp1(i) = dgpvp1(i)+dp1pt*grat(i)
                dgpvp2(i) = dgpvp2(i)+dp2pt*grat(i)
                dgpvt(i) = dtpp2*grap2(i)+dtpp1*grap1(i)+dtpt*grat(i)
            end if
            dgpgp1(1) = dp12p1
            dgpgp2(1) = dp12p2
            if (ds_thm%ds_elem%l_dof_ther) then
                dgpgt(1) = dp12t
            end if
!
! --------- GRADRVP - Derivative of gradient of vapor density and as for gas density
!
            dgrvp1(i) = masvrt*dgpvp1(i)
            dgrvp2(i) = masvrt*dgpvp2(i)
            dgrasp1(i) = -masart*dgpvp1(i)
            dgrasp2(i) = -masart*dgpvp2(i)

            dgrvgp1 = masvrt*dgpgp1(1)
            dgrvgp2 = masvrt*dgpgp2(1)
            dgrasgp1 = -masart*dgpgp1(1)
            dgrasgp2 = -masart*dgpgp2(1)+masart
            if (ds_thm%ds_elem%l_dof_ther) then
                dgrvp1(i) = dgrvp1(i)-grat(i)/t*dr12p1
                dgrvp2(i) = dgrvp2(i)-grat(i)/t*dr12p2
                dgrvt(i) = masvrt*dgpvt(i)-masvrt/t*gpv(i)-grat(i)/t*dr12t+grat(i)/t/t*rho12
                dgrvgt = masvrt*dgpgt(1)-rho12/t
!
                dgrasp1(i) = dgrasp1(i)+grat(i)/t/t*dr21p1
                dgrasp2(i) = dgrasp2(i)+grat(i)/t/t*dr21p2
                dgrast(i) = -masart*dgpvt(i)-masart/t*(grap2(i)-gpv(i))+grat(i)/t/t*dr21t
                dgrasgt = -masart*dgpgt(1)-rho21/t
            end if
!
! --------- GRADCVP - Derivative
            dgcvp1(i) = dgrvp1(i)/rhog+ &
                        (2*rho12*(gr12(i)+gr21(i))/rhog/rhog/rhog-gr12(i)/rhog/rhog) &
                        *(dr12p1+dr21p1) &
                        -(gr12(i)+gr21(i))/rhog/rhog*dr12p1 &
                        -rho12/rhog/rhog*(dgrvp1(i)+dgrasp1(i))
            dgcvp2(i) = dgrvp2(i)/rhog+ &
                        (2*rho12*(gr12(i)+gr21(i))/rhog/rhog/rhog-gr12(i)/rhog/rhog) &
                        *(dr12p2+dr21p2)-(gr12(i)+gr21(i))/rhog/rhog*dr12p2 &
                        -rho12/rhog/rhog*(dgrvp2(i)+dgrasp2(i))
!
            if (ds_thm%ds_elem%l_dof_ther) then
!
                dgcvt(i) = dgrvt(i)/rhog+ &
                           (2*rho12*(gr12(i)+gr21(i))/rhog/rhog/rhog-gr12(i)/rhog/rhog) &
                           *(dr12t+dr21t)-(gr12(i)+gr21(i))/rhog/rhog*dr12t &
                           -rho12/rhog/rhog*(dgrvt(i)+dgrast(i))
            end if
            dgcgp1(1) = dgrvgp1/rhog-rho12/rhog/rhog*(dgrvgp1+dgrasgp1)
            dgcgp2(1) = dgrvgp2/rhog-rho12/rhog/rhog*(dgrvgp2+dgrasgp2)
!
            if (ds_thm%ds_elem%l_dof_ther) then
                dgcgt(1) = dgrvgt/rhog-rho12/rhog/rhog*(dgrvgt+dgrasgt)
            end if
        end do
    end if
!
! - Hydraulic flux
!
    if (lSigm) then
        do i = 1, ndim
            congep(adcp11+i) = 0.d0
            congep(adcp12+i) = -rhog*fv(1)*gc(i)
            congep(adcp21+i) = +rhog*fv(1)*gc(i)
            do j = 1, ndim
                congep(adcp11+i) = congep(adcp11+i)+ &
                                   rho11*lambd1(1)*tperm(i, j)* &
                                   (-grap2(j)+grap1(j)+rho11*gravity(j))
                congep(adcp12+i) = congep(adcp12+i)+ &
                                   rho12*lambd2(1)*tperm(i, j)* &
                                   (-grap2(j)+(rho12+rho21)*gravity(j))
                congep(adcp21+i) = congep(adcp21+i)+ &
                                   rho21*lambd2(1)*tperm(i, j)* &
                                   (-grap2(j)+(rho12+rho21)*gravity(j))
            end do
        end do
    end if
!
! - Update matrix
!
    if (lMatr) then
        do i = 1, ndim
            do j = 1, ndim
                dsde(adcp11+i, addep1) = dsde(adcp11+i, addep1)+ &
                                         dr11p1*lambd1(1)*tperm(i, j)* &
                                         (-grap2(j)+grap1(j)+rho11*gravity(j))
                dsde(adcp11+i, addep1) = dsde(adcp11+i, addep1)+ &
                                         rho11*lambd1(3)*tperm(i, j)* &
                                         (-grap2(j)+grap1(j)+rho11*gravity(j))
                dsde(adcp11+i, addep1) = dsde(adcp11+i, addep1)+ &
                                         rho11*lambd1(1)*tperm(i, j)* &
                                         (dr11p1*gravity(j))
                dsde(adcp11+i, addep2) = dsde(adcp11+i, addep2)+ &
                                         dr11p2*lambd1(1)*tperm(i, j)* &
                                         (-grap2(j)+grap1(j)+rho11*gravity(j))
                dsde(adcp11+i, addep2) = dsde(adcp11+i, addep2)+ &
                                         rho11*lambd1(4)*tperm(i, j)* &
                                         (-grap2(j)+grap1(j)+rho11*gravity(j))
                dsde(adcp11+i, addep2) = dsde(adcp11+i, addep2)+ &
                                         rho11*lambd1(1)*tperm(i, j)*(dr11p2*gravity(j))
                dsde(adcp11+i, addep1+j) = dsde(adcp11+i, addep1+j)+ &
                                           rho11*lambd1(1)*tperm(i, j)
                dsde(adcp11+i, addep2+j) = dsde(adcp11+i, addep2+j)- &
                                           rho11*lambd1(1)*tperm(i, j)
            end do
            do j = 1, ndim
                dsde(adcp12+i, addep1) = dsde(adcp12+i, addep1)+ &
                                         dr12p1*lambd2(1)*tperm(i, j)* &
                                         (-grap2(j)+(rho12+rho21)*gravity(j))
                dsde(adcp12+i, addep1) = dsde(adcp12+i, addep1)+ &
                                         rho12*lambd2(3)*tperm(i, j)* &
                                         (-grap2(j)+(rho12+rho21)*gravity(j))
                dsde(adcp12+i, addep1) = dsde(adcp12+i, addep1)+ &
                                         rho12*lambd2(1)*tperm(i, j)* &
                                         ((dr12p1+dr21p1)*gravity(j))
            end do
            dsde(adcp12+i, addep1) = dsde(adcp12+i, addep1)-(dr12p1+dr21p1)*fv(1)*gc(i)
            dsde(adcp12+i, addep1) = dsde(adcp12+i, addep1)-rhog*fv(3)*gc(i)
            dsde(adcp12+i, addep1) = dsde(adcp12+i, addep1)-rhog*fv(1)*dgcvp1(i)

            do j = 1, ndim
                dsde(adcp12+i, addep2) = dsde(adcp12+i, addep2)+ &
                                         dr12p2*lambd2(1)*tperm(i, j)* &
                                         (-grap2(j)+(rho12+rho21)*gravity(j))
                dsde(adcp12+i, addep2) = dsde(adcp12+i, addep2)+ &
                                         rho12*lambd2(4)*tperm(i, j)* &
                                         (-grap2(j)+(rho12+rho21)*gravity(j))
                dsde(adcp12+i, addep2) = dsde(adcp12+i, addep2)+ &
                                         rho12*lambd2(1)*tperm(i, j)* &
                                         ((dr12p2+dr21p2)*gravity(j))
            end do
            dsde(adcp12+i, addep2) = dsde(adcp12+i, addep2)-(dr12p2+dr21p2)*fv(1)*gc(i)
            dsde(adcp12+i, addep2) = dsde(adcp12+i, addep2)-rhog*fv(4)*gc(i)
            dsde(adcp12+i, addep2) = dsde(adcp12+i, addep2)-rhog*fv(1)*dgcvp2(i)
            dsde(adcp12+i, addep1+i) = dsde(adcp12+i, addep1+i)-rhog*fv(1)*dgcgp1(1)
            do j = 1, ndim
                dsde(adcp12+i, addep2+j) = dsde(adcp12+i, addep2+j)- &
                                           rho12*lambd2(1)*tperm(i, j)
            end do
            dsde(adcp12+i, addep2+i) = dsde(adcp12+i, addep2+i)-rhog*fv(1)*dgcgp2(1)

            do j = 1, ndim
                dsde(adcp21+i, addep1) = dsde(adcp21+i, addep1)+ &
                                         dr21p1*lambd2(1)*tperm(i, j)* &
                                         (-grap2(j)+(rho12+rho21)*gravity(j))
                dsde(adcp21+i, addep1) = dsde(adcp21+i, addep1)+ &
                                         rho21*lambd2(3)*tperm(i, j)* &
                                         (-grap2(j)+(rho12+rho21)*gravity(j))
                dsde(adcp21+i, addep1) = dsde(adcp21+i, addep1)+ &
                                         rho21*lambd2(1)*tperm(i, j)* &
                                         ((dr12p1+dr21p1)*gravity(j))
            end do
            dsde(adcp21+i, addep1) = dsde(adcp21+i, addep1)+(dr12p1+dr21p1)*fv(1)*gc(i)
            dsde(adcp21+i, addep1) = dsde(adcp21+i, addep1)+rhog*fv(3)*gc(i)
            dsde(adcp21+i, addep1) = dsde(adcp21+i, addep1)+rhog*fv(1)*dgcvp1(i)
            do j = 1, ndim
                dsde(adcp21+i, addep2) = dsde(adcp21+i, addep2)+ &
                                         dr21p2*lambd2(1)*tperm(i, j)* &
                                         (-grap2(j)+(rho12+rho21)*gravity(j))
                dsde(adcp21+i, addep2) = dsde(adcp21+i, addep2)+ &
                                         rho21*lambd2(4)*tperm(i, j)* &
                                         (-grap2(j)+(rho12+rho21)*gravity(j))
                dsde(adcp21+i, addep2) = dsde(adcp21+i, addep2)+ &
                                         rho21*lambd2(1)*tperm(i, j)* &
                                         ((dr12p2+dr21p2)*gravity(j))
            end do
            dsde(adcp21+i, addep2) = dsde(adcp21+i, addep2)+(dr12p2+dr21p2)*fv(1)*gc(i)
            dsde(adcp21+i, addep2) = dsde(adcp21+i, addep2)+rhog*fv(4)*gc(i)
            dsde(adcp21+i, addep2) = dsde(adcp21+i, addep2)+rhog*fv(1)*dgcvp2(i)
            dsde(adcp21+i, addep1+i) = dsde(adcp21+i, addep1+i)+rhog*fv(1)*dgcgp1(1)
            do j = 1, ndim
                dsde(adcp21+i, addep2+j) = dsde(adcp21+i, addep2+j)- &
                                           rho21*lambd2(1)*tperm(i, j)
            end do
            dsde(adcp21+i, addep2+i) = dsde(adcp21+i, addep2+i)+rhog*fv(1)*dgcgp2(1)
            if (ds_thm%ds_elem%l_dof_meca) then
                do j = 1, 3
                    do k = 1, ndim
                        dsde(adcp11+i, addeme+ndim-1+j) = &
                            dsde(adcp11+i, addeme+ndim-1+j)+ &
                            rho11*lambd1(2)*tperm(i, k)* &
                            (-grap2(k)+grap1(k)+rho11*gravity(k))
                        dsde(adcp12+i, addeme+ndim-1+j) = &
                            dsde(adcp12+i, addeme+ndim-1+j)+ &
                            rho12*lambd2(2)*tperm(i, k)* &
                            (-grap2(k)+(rho12+rho21)*gravity(k))
                    end do
                    dsde(adcp12+i, addeme+ndim-1+j) = dsde(adcp12+i, addeme+ndim-1+j)- &
                                                      rhog*fv(2)*gc(i)
                    do k = 1, ndim
                        dsde(adcp21+i, addeme+ndim-1+j) = &
                            dsde(adcp21+i, addeme+ndim-1+j)+ &
                            rho21*lambd2(2)*tperm(i, k)* &
                            (-grap2(k)+(rho12+rho21)*gravity(k))
                    end do
                    dsde(adcp21+i, addeme+ndim-1+j) = &
                        dsde(adcp21+i, addeme+ndim-1+j)+ &
                        rhog*fv(2)*gc(i)
                end do
            end if
            if (ds_thm%ds_elem%l_dof_ther) then
                do j = 1, ndim
                    dsde(adcp11+i, addete) = dsde(adcp11+i, addete)+ &
                                             dr11t*lambd1(1)*tperm(i, j)* &
                                             (-grap2(j)+grap1(j)+rho11*gravity(j))
                    dsde(adcp11+i, addete) = dsde(adcp11+i, addete)+ &
                                             rho11*lambd1(5)*tperm(i, j)* &
                                             (-grap2(j)+grap1(j)+rho11*gravity(j))
                    dsde(adcp11+i, addete) = dsde(adcp11+i, addete)+ &
                                             rho11*lambd1(1)*tperm(i, j)* &
                                             ((dr11t)*gravity(j))
                    dsde(adcp12+i, addete) = dsde(adcp12+i, addete)+ &
                                             dr12t*lambd2(1)*tperm(i, j)* &
                                             (-grap2(j)+(rho12+rho21)*gravity(j))
                    dsde(adcp12+i, addete) = dsde(adcp12+i, addete)+ &
                                             rho12*lambd2(5)*tperm(i, j)* &
                                             (-grap2(j)+(rho12+rho21)*gravity(j))
                    dsde(adcp12+i, addete) = dsde(adcp12+i, addete)+ &
                                             rho12*lambd2(1)*tperm(i, j)* &
                                             ((dr12t+dr21t)*gravity(j))
                end do
                dsde(adcp12+i, addete) = dsde(adcp12+i, addete)-(dr12t+dr21t)*fv(1)*gc(i)
                dsde(adcp12+i, addete) = dsde(adcp12+i, addete)-rhog*fv(5)*gc(i)
                dsde(adcp12+i, addete) = dsde(adcp12+i, addete)-rhog*fv(1)*dgcvt(i)
                dsde(adcp12+i, addete+i) = dsde(adcp12+i, addete+i)-rhog*fv(1)*dgcgt(1)

                do j = 1, ndim
                    dsde(adcp21+i, addete) = dsde(adcp21+i, addete)+ &
                                             dr21t*lambd2(1)*tperm(i, j)* &
                                             (-grap2(j)+(rho12+rho21)*gravity(j))
                    dsde(adcp21+i, addete) = dsde(adcp21+i, addete)+ &
                                             rho21*lambd2(5)*tperm(i, j)* &
                                             (-grap2(j)+(rho12+rho21)*gravity(j))
                    dsde(adcp21+i, addete) = dsde(adcp21+i, addete)+ &
                                             rho21*lambd2(1)*tperm(i, j)* &
                                             ((dr12t+dr21t)*gravity(j))
                end do
                dsde(adcp21+i, addete) = dsde(adcp21+i, addete)+(dr12t+dr21t)*fv(1)*gc(i)
                dsde(adcp21+i, addete) = dsde(adcp21+i, addete)+rhog*fv(5)*gc(i)
                dsde(adcp21+i, addete) = dsde(adcp21+i, addete)+rhog*fv(1)*dgcvt(i)
                dsde(adcp21+i, addete+i) = dsde(adcp21+i, addete+i)+rhog*fv(1)*dgcgt(1)

            end if
        end do
    end if
!
end subroutine
