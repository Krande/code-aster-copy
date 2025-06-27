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
!
interface
subroutine nocart(carte, code, ncmp, groupma, mode, nma,&
                  limano, limanu, ligrel,&
                  rapide,jdesc,jnoma,jncmp,jnoli,jvale,&
                  jvalv,jnocmp,ncmpmx,nec, ctype,&
                  jlima0,jlimac,lontav)

    character(len=*), intent(in) :: carte
    integer(kind=8), intent(in) :: code
    integer(kind=8), intent(in) :: ncmp
    character(len=*), intent(in), optional :: groupma
    character(len=*),intent(in), optional :: mode
    integer(kind=8), intent(in), optional :: nma
    character(len=*), intent(in), optional :: limano(*)
    integer(kind=8), intent(in), optional :: limanu(*)
    character(len=*), intent(in), optional ::  ligrel

!   -- arguments optionnels pour gagner du CPU :
    character(len=3), intent(in), optional ::  rapide
    integer(kind=8), intent(inout), optional ::  jdesc
    integer(kind=8), intent(inout), optional ::  jnoma
    integer(kind=8), intent(inout), optional ::  jncmp
    integer(kind=8), intent(inout), optional ::  jnoli
    integer(kind=8), intent(inout), optional ::  jvale
    integer(kind=8), intent(inout), optional ::  jvalv
    integer(kind=8), intent(in)   , optional ::  jnocmp
    integer(kind=8), intent(in)   , optional ::  ncmpmx
    integer(kind=8), intent(in)   , optional ::  nec
    character(len=8), intent(in), optional ::  ctype
    integer(kind=8), intent(inout), optional ::  jlima0
    integer(kind=8), intent(inout), optional ::  jlimac
    integer(kind=8), intent(inout), optional ::  lontav

end subroutine nocart
end interface
