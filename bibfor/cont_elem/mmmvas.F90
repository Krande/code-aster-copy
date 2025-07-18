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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine mmmvas(ndim, nne, nnm, nnl, nbdm, nbcps, &
                  vectce, vectcm, vectfe, vectfm, &
                  vectcc, vectff, &
                  vcont, vfric)
!
    implicit none
!
    integer(kind=8), intent(in) :: nbdm, ndim, nnl, nne, nnm, nbcps
    real(kind=8), intent(in) :: vectce(27), vectcm(27)
    real(kind=8), intent(in) :: vectfe(27), vectfm(27)
    real(kind=8), intent(in) :: vectcc(9), vectff(18)
    real(kind=8), intent(inout) :: vcont(81), vfric(81)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Assembling vectors
!
! --------------------------------------------------------------------------------------------------
!
! In  ndim             : dimension of problem (2 or 3)
! In  nne              : number of slave nodes
! In  nnm              : number of master nodes
! In  nnl              : number of nodes with Lagrange multiplicators (contact and friction)
! In  nbdm             : number of components by node for all dof
! In  nbcps            : number of components by node for Lagrange multiplicators
! In  vectcc           : vector for DOF [contact]
! In  vectff           : vector for DOF [friction]
! In  vectce           : vector for DOF [slave] - For contact part
! In  vectfe           : vector for DOF [slave] - For friction part
! In  vectcm           : vector for DOF [master] - For contact part
! In  vectfm           : vector for DOF [master] - For friction part
! IO  vcont            : resultant vector for contact
! IO  vfric            : resultant vector for friction
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ii, jj
    integer(kind=8) :: inoe, inom, inoc, inof, idim, icmp
    integer(kind=8) :: nbcpf
!
! --------------------------------------------------------------------------------------------------
!
    do inoe = 1, nne
        do idim = 1, ndim
            ii = nbdm*(inoe-1)+idim
            jj = ndim*(inoe-1)+idim
            vcont(ii) = vcont(ii)+vectce(jj)
            vfric(ii) = vfric(ii)+vectfe(jj)
        end do
    end do
!
    do inom = 1, nnm
        do idim = 1, ndim
            ii = nbdm*nne+ndim*(inom-1)+idim
            jj = ndim*(inom-1)+idim
            vcont(ii) = vcont(ii)+vectcm(jj)
            vfric(ii) = vfric(ii)+vectfm(jj)
        end do
    end do
!
    do inoc = 1, nnl
        ii = nbdm*(inoc-1)+ndim+1
        jj = inoc
        vcont(ii) = vcont(ii)+vectcc(jj)
    end do
!
    nbcpf = nbcps-1
    do inof = 1, nnl
        do icmp = 1, nbcpf
            ii = nbdm*(inof-1)+ndim+1+icmp
            jj = nbcpf*(inof-1)+icmp
            vfric(ii) = vfric(ii)+vectff(jj)
        end do
    end do
!
end subroutine
