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
subroutine posthm(option, modint, jgano, ncmp, nvim, &
                  vpg, vno)
    implicit none
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/ppgan2.h"
    integer(kind=8) :: jgano, ncmp, nvim
    real(kind=8) :: vno(*), vpg(*)
    character(len=3) :: modint
    character(len=16) :: option
! --- ROUTINE DE POST-TRAITEMENT POUR LA THM --------------------------
! --- ON DEFINIT LES VECTEURS SUIVANT : -------------------------------
! =====================================================================
! --- * VPG1 DE DIMENSION NPGMAX*NVMAX - OU NVMAX EST LE --------------
! ---   NOMBRE MAXIMAL DE VARIABLE INTERNE (A CE JOUR : 20 : 16 POUR --
! ---   LA LOI CJS ET 4 POUR LES VINT LIEES A LA THM NON MECANIQUE) ---
! ---   NPGMAX ETANT LE NOMBRE DE POINTS DE GAUSS MAXIMAL -------------
! ---   (A CE JOUR : 8 : 8 POINTS DE GAUSS POUR LES HEXA20) -----------
! =====================================================================
! --- * VPG2 DE DIMENSION NNOSMA*NVMAX - OU ---------------------------
! ---   NNOSMA ETANT LE NOMBRE DE NOEUDS SOMMETS MAXIMAL --------------
! ---   (A CE JOUR : 8 : 8 SOMMETS POUR LES HEXA20) -------------------
! =====================================================================
! --- * SPG1 DE DIMENSION NPGMAX*DIMMAX - OU DIMMAX -------------------
! ---   EST LA DIMENSION MAXIMAL DE DIMCON (A CE JOUR : 31+5) ----------
! =====================================================================
! --- * SPG2 DE DIMENSION NNOSMA*DIMMAX -------------------------------
! =====================================================================
! --- * VARIPG ET VARISO DE DIMENSION NNOMAX*NVMAX - OU NNOMAX EST LE -
! ---   NOMBRE MAXIMAL DE NOEUDS (A CE JOUR : 20 : 20 SOMMETS POUR ----
! ---   LES HEXA20) ---------------------------------------------------
! =====================================================================
! --- * SIEFPG ET SIEFSO DE DIMENSION NNOMAX*DIMMAX -------------------
! =====================================================================
    integer(kind=8) :: i, j, jganpg, jganso
    integer(kind=8) :: ndim, nno, nnos, npg, ndim2, nno2, nnos2, npg2
    integer(kind=8) :: nvmax, npgmax, nnosma, dimmax, nnomax
    parameter(nvmax=60)
    parameter(npgmax=8)
    parameter(nnosma=8)
    parameter(dimmax=31+5)
    parameter(nnomax=20)
    real(kind=8) :: vpg1(npgmax*nvmax), vpg2(nnosma*nvmax)
    real(kind=8) :: spg1(npgmax*dimmax), spg2(nnosma*dimmax)
    real(kind=8) :: varipg(nnomax*nvmax), variso(nnomax*nvmax)
    real(kind=8) :: siefpg(nnomax*dimmax), siefso(nnomax*dimmax)
! =====================================================================
    if (modint .ne. 'RED') then
        call ppgan2(jgano, 1, ncmp, vpg, vno)
    else
! =====================================================================
! --- MATRICE DE PASSAGE POINTS DE GAUSS -> SOMMETS JGANPG ------------
! =====================================================================
        call elrefe_info(fami='MASS', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jgano=jganpg)
! =====================================================================
! --- MATRICE DE PASSAGE SOMMETS -> SOMMETS : JGANSO ------------------
! =====================================================================
        call elrefe_info(fami='NOEU_S', ndim=ndim2, nno=nno2, nnos=nnos2, npg=npg2, &
                         jgano=jganso)
! =====================================================================
! --- ON VERIFIE QUE LES DIMENSIONNEMENTS SONT A JOUR -----------------
! =====================================================================
        ASSERT(nno .le. nnomax)
        ASSERT(npg .le. npgmax)
        ASSERT(nnos .le. nnosma)
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
            call ppgan2(jganpg, 1, ncmp, spg1, siefpg)
            call ppgan2(jganso, 1, ncmp, spg2, siefso)
            do i = 1, nno
                do j = 1, nvim
                    vno((i-1)*ncmp+j) = siefpg((i-1)*ncmp+j)
                end do
                do j = nvim+1, ncmp
                    vno((i-1)*ncmp+j) = siefso((i-1)*ncmp+j)
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
            call ppgan2(jganpg, 1, ncmp, vpg1, varipg)
            call ppgan2(jganso, 1, ncmp, vpg2, variso)
            do i = 1, nno
                do j = 1, nvim
                    vno((i-1)*ncmp+j) = varipg((i-1)*ncmp+j)
                end do
                do j = nvim+1, ncmp
                    vno((i-1)*ncmp+j) = variso((i-1)*ncmp+j)
                end do
            end do
        end if
    end if
! =====================================================================
end subroutine
