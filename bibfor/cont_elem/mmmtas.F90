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
! aslint: disable=W1504
!
subroutine mmmtas(nbdm, ndim, nnl, nne, nnm, nbcps, &
                  matrcc, matree, matrmm, matrem, &
                  matrme, matrce, matrcm, matrmc, matrec, &
                  matr, &
                  matrff_, matrfe_, matrfm_, matrmf_, matref_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/mmmtdb.h"
!
    integer(kind=8), intent(in) :: nbdm, ndim, nnl, nne, nnm, nbcps
    real(kind=8), intent(in) :: matrcc(9, 9), matree(27, 27), matrmm(27, 27)
    real(kind=8), intent(in) :: matrem(27, 27), matrme(27, 27), matrce(9, 27), matrcm(9, 27)
    real(kind=8), intent(in) :: matrec(27, 9), matrmc(27, 9)
    real(kind=8), intent(inout) :: matr(81, 81)
    real(kind=8), optional, intent(in) :: matrff_(18, 18)
    real(kind=8), optional, intent(in) :: matrfe_(18, 27), matrfm_(18, 27)
    real(kind=8), optional, intent(in) :: matrmf_(27, 18), matref_(27, 18)
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
! In  matrcc           : matrix for DOF [contact x contact]
! In  matrff           : matrix for DOF [friction x friction]
! In  matrce           : matrix for DOF [contact x slave]
! In  matrcm           : matrix for DOF [contact x master]
! In  matrfe           : matrix for DOF [friction x slave]
! In  matrfm           : matrix for DOF [friction x master]
! In  matree           : matrix for DOF [slave x slave]
! In  matrmm           : matrix for DOF [master x master]
! In  matrem           : matrix for DOF [slave x master]
! In  matrme           : matrix for DOF [master x slave]
! In  matrec           : matrix for DOF [slave x contact]
! In  matrmc           : matrix for DOF [master x contact]
! In  matref           : matrix for DOF [slave x friction]
! In  matrmf           : matrix for DOF [master x friction]
! IO  matr             : resultant matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ii, jj, kk, ll
    integer(kind=8) :: nbcpf
    integer(kind=8) :: inoc, inoe, inom, inof, icmp, idim
    integer(kind=8) :: inoe1, inoe2, inom1, inom2, idim1, idim2
    integer(kind=8) :: inoc1, inoc2, inof1, inof2, icmp1, icmp2
    aster_logical :: debug
!
! --------------------------------------------------------------------------------------------------
!
    nbcpf = nbcps-1
    debug = ASTER_FALSE
!
    do inoc1 = 1, nnl
        do inoc2 = 1, nnl
            ii = nbdm*(inoc1-1)+ndim+1
            jj = nbdm*(inoc2-1)+ndim+1
            kk = inoc1
            ll = inoc2
            matr(ii, jj) = matr(ii, jj)+matrcc(kk, ll)
            if (debug) call mmmtdb(matrcc(kk, ll), 'CC', ii, jj)
        end do
    end do
!
    do inoc = 1, nnl
        do inoe = 1, nne
            do idim = 1, ndim
                ii = nbdm*(inoc-1)+ndim+1
                jj = nbdm*(inoe-1)+idim
                kk = inoc
                ll = ndim*(inoe-1)+idim
                matr(ii, jj) = matr(ii, jj)+matrce(kk, ll)
                if (debug) call mmmtdb(matrce(kk, ll), 'CE', ii, jj)
            end do
        end do
    end do
!
    do inoc = 1, nnl
        do inoe = 1, nne
            do idim = 1, ndim
                ii = nbdm*(inoc-1)+ndim+1
                jj = nbdm*(inoe-1)+idim
                kk = inoc
                ll = ndim*(inoe-1)+idim
                matr(jj, ii) = matr(jj, ii)+matrec(ll, kk)
                if (debug) call mmmtdb(matrec(ll, kk), 'EC', jj, ii)
            end do
        end do
    end do
!
    do inoc = 1, nnl
        do inom = 1, nnm
            do idim = 1, ndim
                ii = nbdm*(inoc-1)+ndim+1
                jj = nbdm*nne+ndim*(inom-1)+idim
                kk = inoc
                ll = ndim*(inom-1)+idim
                matr(ii, jj) = matr(ii, jj)+matrcm(kk, ll)
                if (debug) call mmmtdb(matrcm(kk, ll), 'CM', ii, jj)
            end do
        end do
    end do
!
    do inoc = 1, nnl
        do inom = 1, nnm
            do idim = 1, ndim
                ii = nbdm*(inoc-1)+ndim+1
                jj = nbdm*nne+ndim*(inom-1)+idim
                kk = ndim*(inom-1)+idim
                ll = inoc
                matr(jj, ii) = matr(jj, ii)+matrmc(kk, ll)
                if (debug) call mmmtdb(matrmc(kk, ll), 'MC', jj, ii)
            end do
        end do
    end do
!
    do inoe1 = 1, nne
        do inoe2 = 1, nne
            do idim2 = 1, ndim
                do idim1 = 1, ndim
                    ii = nbdm*(inoe1-1)+idim1
                    jj = nbdm*(inoe2-1)+idim2
                    kk = ndim*(inoe1-1)+idim1
                    ll = ndim*(inoe2-1)+idim2
                    matr(ii, jj) = matr(ii, jj)+matree(kk, ll)
                    if (debug) call mmmtdb(matree(kk, ll), 'EE', ii, jj)
                end do
            end do
        end do
    end do
!
    do inom1 = 1, nnm
        do inom2 = 1, nnm
            do idim2 = 1, ndim
                do idim1 = 1, ndim
                    kk = ndim*(inom1-1)+idim1
                    ll = ndim*(inom2-1)+idim2
                    ii = nbdm*nne+ndim*(inom1-1)+idim1
                    jj = nbdm*nne+ndim*(inom2-1)+idim2
                    matr(ii, jj) = matr(ii, jj)+matrmm(kk, ll)
                    if (debug) call mmmtdb(matrmm(kk, ll), 'MM', ii, jj)
                end do
            end do
        end do
    end do
!
    do inoe = 1, nne
        do inom = 1, nnm
            do idim2 = 1, ndim
                do idim1 = 1, ndim
                    kk = ndim*(inoe-1)+idim1
                    ll = ndim*(inom-1)+idim2
                    ii = nbdm*(inoe-1)+idim1
                    jj = nbdm*nne+ndim*(inom-1)+idim2
                    matr(ii, jj) = matr(ii, jj)+matrem(kk, ll)
                    if (debug) call mmmtdb(matrem(kk, ll), 'EM', ii, jj)
                end do
            end do
        end do
    end do
!
    do inoe = 1, nne
        do inom = 1, nnm
            do idim2 = 1, ndim
                do idim1 = 1, ndim
                    ii = nbdm*nne+ndim*(inom-1)+idim2
                    jj = nbdm*(inoe-1)+idim1
                    kk = ndim*(inom-1)+idim2
                    ll = ndim*(inoe-1)+idim1
                    matr(ii, jj) = matr(ii, jj)+matrme(kk, ll)
                    if (debug) call mmmtdb(matrme(kk, ll), 'ME', ii, jj)
                end do
            end do
        end do
    end do
!
! - For friction
!
    if (nbcpf .gt. 0) then
!
        do inof1 = 1, nnl
            do inof2 = 1, nnl
                do icmp1 = 1, nbcpf
                    do icmp2 = 1, nbcpf
                        ii = nbdm*(inof1-1)+ndim+1+icmp1
                        jj = nbdm*(inof2-1)+ndim+1+icmp2
                        kk = (ndim-1)*(inof1-1)+icmp1
                        ll = (ndim-1)*(inof2-1)+icmp2
                        matr(ii, jj) = matr(ii, jj)+matrff_(kk, ll)
                        if (debug) call mmmtdb(matrff_(kk, ll), 'FF', ii, jj)
                    end do
                end do
            end do
        end do
!
        do inof = 1, nnl
            do inoe = 1, nne
                do icmp = 1, nbcpf
                    do idim = 1, ndim
                        kk = nbcpf*(inof-1)+icmp
                        ll = ndim*(inoe-1)+idim
                        ii = nbdm*(inof-1)+ndim+1+icmp
                        jj = nbdm*(inoe-1)+idim
                        matr(ii, jj) = matr(ii, jj)+matrfe_(kk, ll)
                        if (debug) call mmmtdb(matrfe_(kk, ll), 'FE', ii, jj)
                    end do
                end do
            end do
        end do
!
        do inoe = 1, nne
            do inof = 1, nnl
                do icmp = 1, nbcpf
                    do idim = 1, ndim
                        kk = ndim*(inoe-1)+idim
                        ll = nbcpf*(inof-1)+icmp
                        ii = nbdm*(inoe-1)+idim
                        jj = nbdm*(inof-1)+ndim+1+icmp
                        matr(ii, jj) = matr(ii, jj)+matref_(kk, ll)
                        if (debug) call mmmtdb(matrfe_(kk, ll), 'FE', ii, jj)
                    end do
                end do
            end do
        end do
!
        do inof = 1, nnl
            do inom = 1, nnm
                do icmp = 1, nbcpf
                    do idim = 1, ndim
                        kk = nbcpf*(inof-1)+icmp
                        ll = ndim*(inom-1)+idim
                        ii = nbdm*(inof-1)+ndim+1+icmp
                        jj = nbdm*nne+ndim*(inom-1)+idim
                        matr(ii, jj) = matr(ii, jj)+matrfm_(kk, ll)
                        if (debug) call mmmtdb(matrfm_(kk, ll), 'FM', ii, jj)
                    end do
                end do
            end do
        end do
!
        do inom = 1, nnm
            do inof = 1, nnl
                do icmp = 1, nbcpf
                    do idim = 1, ndim
                        kk = ndim*(inom-1)+idim
                        ll = nbcpf*(inof-1)+icmp
                        ii = nbdm*nne+ndim*(inom-1)+idim
                        jj = nbdm*(inof-1)+ndim+1+icmp
                        matr(ii, jj) = matr(ii, jj)+matrmf_(kk, ll)
                        if (debug) call mmmtdb(matrmf_(kk, ll), 'MF', ii, jj)
                    end do
                end do
            end do
        end do
    end if
!
end subroutine
