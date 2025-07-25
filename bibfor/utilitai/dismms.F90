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

subroutine dismms(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(MATR_ASSE) (MARCHE AUSSI PARFOIS SUR MATR_ASSE_GENE)
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/dismnu.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"

    integer(kind=8) :: repi, ierd
    character(len=*) :: questi
    character(len=*) :: nomobz, repkz
! ----------------------------------------------------------------------
!    IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE CONCEPT MATR_ASSE  (K19)
!    OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)

! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=32) :: repk
    character(len=24) ::   k24
    character(len=19) :: nomob, solveu
    character(len=2) :: typmat
    character(len=8) :: nommai
!-----------------------------------------------------------------------
    integer(kind=8) ::  ier
    character(len=16) :: typeco
    character(len=24), pointer :: refa(:) => null()
    character(len=24), pointer :: slvk(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    repk = ' '
    repi = 0
    ierd = 0

    nomob = nomobz
    call jeveuo(nomob//'.REFA', 'L', vk24=refa)

    if (questi(1:9) .eq. 'NUM_GD_SI') then
        call dismnu(questi, refa(2) (1:14), repi, repk, ierd)
    else if (questi(1:9) .eq. 'NOM_GD_SI') then
        call dismnu('NOM_GD', refa(2) (1:14), repi, repk, ierd)

    else if (questi .eq. 'TYPE_MATRICE') then
        typmat = refa(9) (1:2)
        if (typmat .eq. 'MS') then
            repk = 'SYMETRI'
        else if (typmat .eq. 'MR') then
            repk = 'NON_SYM'
        else
            ASSERT(.false.)
        end if

    else if (questi .eq. 'NB_EQUA') then
        call dismnu(questi, refa(2) (1:14), repi, repk, ierd)

    else if (questi .eq. 'MATR_DISTRIBUEE') then
        call dismnu(questi, refa(2) (1:14), repi, repk, ierd)
        if (repk .eq. 'OUI') then
            ASSERT(refa(11) .eq. 'MATR_DISTR')
        end if

    else if (questi .eq. 'MATR_HPC') then
        nommai = refa(1) (1:8)
        call gettco(nommai, typeco)
        if (typeco .eq. 'MAILLAGE_P') then
            repk = 'OUI'
        else
            repk = 'NON'
        end if

    else if (questi .eq. 'SOLVEUR') then
        if (refa(7) .ne. ' ') then
            repk = refa(7)
        end if

    else if (questi .eq. 'NOM_MODELE') then
        call dismnu(questi, refa(2) (1:14), repi, repk, ierd)

    else if (questi .eq. 'NOM_MAILLA') then
        repk = refa(1) (1:8)

    else if (questi .eq. 'NOM_NUME_DDL') then
        repk = refa(2) (1:14)

    else if (questi .eq. 'EXIS_LAGR') then
        call jeexin(nomob//'.CONL', ier)
        if (ier .eq. 0) then
            repk = 'NON'
        else
            repk = 'OUI'
        end if

    else if (questi .eq. 'XFEM') then
        repk = refa(17)

    else if (questi .eq. 'XFEM_PC') then
        repk = refa(18) (1:19)

    else if (questi .eq. 'XFEM_PC_INV') then
        repk = refa(16) (1:19)

    else if (questi .eq. 'METH_RESO' .or. questi .eq. 'RENUM_RESO') then
        if (refa(7) .ne. ' ') then
            solveu = refa(7) (1:19)
            call jeveuo(solveu//'.SLVK', 'L', vk24=slvk)
            if (questi .eq. 'METH_RESO') then
                repk = slvk(1)
            else
                repk = slvk(4)
            end if
        end if
    else if (questi .eq. 'NUME_EQUA') then
        repk = refa(2) (1:14)//'.NUME'

    else if (questi .eq. 'PHENOMENE') then
        call dismnu(questi, refa(2) (1:14), repi, repk, ierd)

    else if (questi .eq. 'SUR_OPTION') then
        repk = refa(4) (1:16)

    else if (questi .eq. 'MPI_COMPLET') then
        k24 = refa(11)
        ASSERT(k24 .eq. 'MPI_COMPLET' .or. k24 .eq. 'MPI_INCOMPLET' .or. k24 .eq. 'MATR_DISTR')
        if (k24 .eq. 'MPI_COMPLET') then
            repk = 'OUI'
        else
            repk = 'NON'
        end if

    else if (questi .eq. 'MATR_DISTR') then
        k24 = refa(11)
        ASSERT(k24 .eq. 'MPI_COMPLET' .or. k24 .eq. 'MPI_INCOMPLET' .or. k24 .eq. 'MATR_DISTR')
        if (k24 .eq. 'MATR_DISTR') then
            repk = 'OUI'
        else
            repk = 'NON'
        end if

    else if (questi .eq. 'EXIS_CINE') then
        call jeexin(nomob//'.CCID', ier)
        ! cas MATR_HPC
        nommai = refa(1) (1:8)
        call gettco(nommai, typeco)
        if (typeco .eq. 'MAILLAGE_P') then
            call asmpi_comm_vect('MPI_SUM', 'I', 1, 0, sci=ier)
        end if
        if (ier .eq. 0) then
            repk = 'NON'
        else
            repk = 'OUI'
        end if

    else
        ierd = 1
    end if

    repkz = repk
    call jedema()
end subroutine
