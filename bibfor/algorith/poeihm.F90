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
subroutine poeihm(nomte, option, modint, jgao, nno1, &
                  nno2, ncmp, nvim, vpg, vno)
    implicit none
#include "asterfort/assert.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/ppgan2.h"
    integer(kind=8) :: jgao, ncmp, nvim
    real(kind=8) :: vno(*), vpg(*)
    character(len=3) :: modint
    character(len=8) :: lielrf(10)
    character(len=16) :: option, nomte
!
! --- ROUTINE DE POST-TRAITEMENT JOINT HM --------------------------
! =====================================================================
! IN NOMTE  : NOM DE L'ELEMENT
! IN OPTION : OPTION DE CALCUL
! IN MODINT : MODE D'INTEGRATION
! IN JGAO : MATRICE DE PASSAGE NOEUDS -> POINTS D'INTEGRATION
! IN NNO1 : NOMBRE DE NOEUDS BORD INF ET SUP
! IN NNO2 : NOMBRE DE NOEUDS SEGMENT MILIEU
! IN NCMP :
! IN NVIM :
! IN VPG  : CHAMPS AUX POINTS D'INTEGRATION
! =====================================================================
! OUT VNO : CHAMPS AUX NOEUDS DE L'ELEMENT
! =====================================================================
    integer(kind=8) :: i, j, jgapg1, jgaso1, jgapg2, jgaso2
    integer(kind=8) :: ndim, nno, npg, ndim2, nno2, nnos2, npg2, nno3, nnos3
    integer(kind=8) :: nno1, nnos1
    integer(kind=8) :: nvmax, npgmax, nnosma, dimmax, nnomax, ntrou
    parameter(nvmax=60)
    parameter(npgmax=8)
    parameter(nnosma=8)
    parameter(dimmax=31)
    parameter(nnomax=20)
    real(kind=8) :: vpg1(npgmax*nvmax), vpg2(nnosma*nvmax)
    real(kind=8) :: spg1(npgmax*dimmax), spg2(nnosma*dimmax)
    real(kind=8) :: varpg1(nnomax*nvmax), varso1(nnomax*nvmax)
    real(kind=8) :: varpg2(nnomax*nvmax), varso2(nnomax*nvmax)
    real(kind=8) :: sefpg1(nnomax*dimmax), sefso1(nnomax*dimmax)
    real(kind=8) :: sefpg2(nnomax*dimmax), sefso2(nnomax*dimmax)
    real(kind=8) :: vno1(nnomax*dimmax)
    integer(kind=8) :: next(3), next2(3), nmil(2)
!
    data next/1, 2, 5/
    data next2/4, 3, 7/
    data nmil/8, 6/
! =====================================================================
    if (modint .ne. 'RED') then
        call ppgan2(jgao, 1, ncmp, vpg, vno1)
!
        do i = 1, nno1
            do j = 1, nvim
                vno((next(i)-1)*ncmp+j) = vno1((i-1)*ncmp+j)
                vno((next2(i)-1)*ncmp+j) = vno1((i-1)*ncmp+j)
            end do
            do j = nvim+1, ncmp
                vno((next(i)-1)*ncmp+j) = vno1((i-1)*ncmp+j)
                vno((next2(i)-1)*ncmp+j) = vno1((i-1)*ncmp+j)
            end do
        end do
        do i = 1, nno2
            do j = 1, nvim
                vno((nmil(i)-1)*ncmp+j) = vno1((i-1)*ncmp+j)
            end do
            do j = nvim+1, ncmp
                vno((nmil(i)-1)*ncmp+j) = vno1((i-1)*ncmp+j)
            end do
        end do
!
    else
!
! =====================================================================
! --- MATRICE DE PASSAGE POINTS DE GAUSS -> SOMMETS JGAPG ------------
! =====================================================================
        call elref2(nomte, 2, lielrf, ntrou)
!
!
        call elrefe_info(elrefe=lielrf(1), fami='MASS', ndim=ndim, nno=nno1, nnos=nnos1, &
                         npg=npg, jgano=jgapg1)
!
        call elrefe_info(elrefe=lielrf(2), fami='MASS', ndim=ndim, nno=nno2, nnos=nnos2, &
                         npg=npg, jgano=jgapg2)
! =====================================================================
! --- MATRICE DE PASSAGE SOMMETS -> SOMMETS : JGASO ------------------
! =====================================================================
        call elrefe_info(elrefe=lielrf(1), fami='NOEU_S', ndim=ndim2, nno=nno3, nnos=nnos3, &
                         npg=npg2, jgano=jgaso1)
!
        call elrefe_info(elrefe=lielrf(2), fami='NOEU_S', ndim=ndim2, nno=nno3, nnos=nnos3, &
                         npg=npg2, jgano=jgaso2)
!
        nno = 2*nno1+nno2
! =====================================================================
! --- ON VERIFIE QUE LES DIMENSIONNEMENTS SONT A JOUR -----------------
! =====================================================================
        ASSERT(nno .le. nnomax)
        ASSERT(npg .le. npgmax)
        ASSERT(nnos1 .le. nnosma)
        if (option .eq. 'SIEF_ELNO  ') then
! =====================================================================
! --- ON VERIFIE QUE LES DIMENSIONNEMENTS SONT A JOUR -----------------
! =====================================================================
            ASSERT(ncmp .le. dimmax)
            do i = 1, ncmp*npg
                spg1(i) = vpg(i)
            end do
            do i = 1, ncmp*npg2
                spg2(i) = vpg(ncmp*npg+i)
            end do
            call ppgan2(jgapg1, 1, ncmp, spg1, sefpg1)
            call ppgan2(jgaso1, 1, ncmp, spg2, sefso1)
            do i = 1, nno1
                do j = 1, nvim
                    vno((next(i)-1)*ncmp+j) = sefpg1((i-1)*ncmp+j)
                    vno((next2(i)-1)*ncmp+j) = sefpg1((i-1)*ncmp+j)
                end do
                do j = nvim+1, ncmp
                    vno((next(i)-1)*ncmp+j) = sefso1((i-1)*ncmp+j)
                    vno((next2(i)-1)*ncmp+j) = sefso1((i-1)*ncmp+j)
                end do
            end do
            call ppgan2(jgapg2, 1, ncmp, spg1, sefpg2)
            call ppgan2(jgaso2, 1, ncmp, spg2, sefso2)
            do i = 1, nno2
                do j = 1, nvim
                    vno((nmil(i)-1)*ncmp+j) = sefpg1((i-1)*ncmp+j)
                end do
                do j = nvim+1, ncmp
                    vno((nmil(i)-1)*ncmp+j) = sefso1((i-1)*ncmp+j)
                end do
            end do
        end if
        if (option .eq. 'VARI_ELNO  ') then
! =====================================================================
! --- ON VERIFIE QUE LES DIMENSIONNEMENTS SONT A JOUR -----------------
! =====================================================================
            ASSERT(ncmp .le. nvmax)
            do i = 1, ncmp*npg
                vpg1(i) = vpg(i)
            end do
            do i = 1, ncmp*npg2
                vpg2(i) = vpg(ncmp*npg+i)
            end do
            call ppgan2(jgapg1, 1, ncmp, vpg1, varpg1)
            call ppgan2(jgaso1, 1, ncmp, vpg2, varso1)
            do i = 1, nno1
                do j = 1, nvim
                    vno((next(i)-1)*ncmp+j) = varpg1((i-1)*ncmp+j)
                    vno((next2(i)-1)*ncmp+j) = varpg1((i-1)*ncmp+j)
                end do
                do j = nvim+1, ncmp
                    vno((next(i)-1)*ncmp+j) = varso1((i-1)*ncmp+j)
                    vno((next2(i)-1)*ncmp+j) = varso1((i-1)*ncmp+j)
                end do
            end do
            call ppgan2(jgapg2, 1, ncmp, vpg1, varpg2)
            call ppgan2(jgaso2, 1, ncmp, vpg2, varso2)
            do i = 1, nno2
                do j = 1, nvim
                    vno((nmil(i)-1)*ncmp+j) = varpg1((i-1)*ncmp+j)
                end do
                do j = nvim+1, ncmp
                    vno((nmil(i)-1)*ncmp+j) = varso1((i-1)*ncmp+j)
                end do
            end do
        end if
    end if
! =====================================================================
end subroutine
