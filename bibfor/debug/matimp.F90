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
subroutine matimp(matz, ific, typimz)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/isParallelMatrix.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*) :: matz, typimz
    integer(kind=8) :: ific
! ---------------------------------------------------------------------
! BUT: IMPRIMER UNE MATRICE SUR UN LISTING
! ---------------------------------------------------------------------
!     ARGUMENTS:
! MATZ   IN/JXIN  K19 : MATR_ASSE A IMPRIMER
! IFIC   IN       I   : NUMERO DE L'UNITE LOGIQUE D'IMPRESSION
! TYPIMP IN       K8  : FORMAT DE L'IMPRESSION ' ', 'ASTER' OU 'MATLAB'
!                       SI TYPIMP=' ', TYPIMP='ASTER'
! ---------------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
    integer(kind=8) :: iligl, jcoll, kterm, n, nz, nsmdi, jsmhc, nsmhc
    integer(kind=8) :: jdelg, n1, nvale, jvale, nlong, jval2, nuno, nucmp, k, jcmp
    integer(kind=8) :: iligg, jcolg, jnlogl, coltmp, jprddl, rang, nbproc
    character(len=8) :: nomgd, nocmp, noma, nono, typimp
    character(len=14) :: nonu, localOrGhost
    character(len=1) :: ktyp
    character(len=3) :: mathpc
    character(len=19) :: mat19
    aster_logical :: ltypr, lsym, lmd, l_parallel_matrix
    real(kind=8) :: dble, dimag
    integer(kind=8), pointer :: deeq(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    character(len=24), pointer :: refa(:) => null()
    character(len=24), pointer :: refn(:) => null()
    mpi_int :: mrank, msize
!
!     ------------------------------------------------------------------
    call jemarq()
!
!
    mat19 = matz
    typimp = typimz
!
    call jeveuo(mat19//'.REFA', 'L', vk24=refa)
    noma = refa(1) (1:8)
    nonu = refa(2) (1:14)
!
    lmd = (refa(11) .eq. 'MATR_DISTR')
!
    l_parallel_matrix = isParallelMatrix(mat19)
    localOrGhost = ' '
!
    call jeveuo(nonu//'.SMOS.SMDI', 'L', vi=smdi)
    call jelira(nonu//'.SMOS.SMDI', 'LONMAX', nsmdi)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)
    call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
    if (lmd) then
        call jeveuo(nonu//'.NUML.DELG', 'L', jdelg)
        call jelira(nonu//'.NUML.DELG', 'LONMAX', n1)
        call jeveuo(nonu//'.NUML.NULG', 'L', jnlogl)
    else
        call jeveuo(nonu//'.NUME.DELG', 'L', jdelg)
        call jelira(nonu//'.NUME.DELG', 'LONMAX', n1)
        jnlogl = 0
        if (l_parallel_matrix) then
            call jeveuo(nonu//'.NUME.PDDL', 'L', jprddl)
            call asmpi_info(rank=mrank, size=msize)
            rang = to_aster_int(mrank)
            nbproc = to_aster_int(msize)
        end if
    end if
    ASSERT(n1 .eq. nsmdi)
!     --- CALCUL DE N
    n = nsmdi
!     --- CALCUL DE NZ
    nz = smdi(n)
!
    ASSERT(nz .le. nsmhc)
    call jelira(mat19//'.VALM', 'NMAXOC', nvale)
    if (nvale .eq. 1) then
        lsym = .true.
    else if (nvale .eq. 2) then
        lsym = .false.
    else
        ASSERT(.false.)
    end if
!
    call jeveuo(jexnum(mat19//'.VALM', 1), 'L', jvale)
    call jelira(jexnum(mat19//'.VALM', 1), 'LONMAX', nlong)
    ASSERT(nlong .eq. nz)
    if (.not. lsym) then
        call jeveuo(jexnum(mat19//'.VALM', 2), 'L', jval2)
        call jelira(jexnum(mat19//'.VALM', 2), 'LONMAX', nlong)
        ASSERT(nlong .eq. nz)
    end if
!
    call jelira(jexnum(mat19//'.VALM', 1), 'TYPE', cval=ktyp)
    ltypr = (ktyp .eq. 'R')
!
!     --- ENTETES
    write (ific, *) ' '
    write (ific, *) '% --------------------------------------------'//&
     &                '----------------------------------------------'
!     --- ENTETE FORMAT ASTER
    if ((typimp .eq. ' ') .or. (typimp .eq. 'ASTER')) then
        write (ific, *) 'DIMENSION DE LA MATRICE :', n
        write (ific, *) 'NOMBRE DE TERMES NON NULS (MATRICE SYMETRIQUE) :'&
     &                , nz
        write (ific, *) 'MATRICE A COEEFICIENTS REELS :', ltypr
        write (ific, *) 'MATRICE SYMETRIQUE :', lsym
        write (ific, *) 'MATRICE DISTRIBUEE :', lmd
        write (ific, *) 'MATRICE HPC :', l_parallel_matrix
        if (l_parallel_matrix) write (ific, *) 'NUMEROTATION LOCALE AU PROCESSUS :', &
         &rang, ' / ', nbproc
        write (ific, *) ' '
!     --- ENTETE FORMAT MATLAB
    else if (typimp .eq. 'MATLAB') then
        write (ific, *) '% IMPRESSION DE LA MATRICE '//mat19//' AU FORMAT'&
     &                //' MATLAB.'
        if (lmd) then
            write (ific, *) '% 0- FUSIONNER LES FICHIERS PROVENANT DES'//&
     &                  ' DIFFERENTS PROCESSEURS (MATR_DISTRIBUEE).'
        end if
        write (ific, *) '% 1- COPIER DANS UN FICHIER mat.dat.'
        write (ific, *) '% 2- CHARGER DANS MATLAB OU OCTAVE PAR :'
        write (ific, *) '% 2-1 >> load -ascii mat.dat;'
        write (ific, *) '% 2-2 >> A=spconvert(mat);'
        write (ific, *) ' '
    else
        ASSERT(.false.)
    end if
!
!
!     ------------------------------------------------
!     IMPRESSION DES TERMES DE LA MATRICE
!     ------------------------------------------------
    if ((typimp .eq. ' ') .or. (typimp .eq. 'ASTER')) write (ific, 1003) 'ILIGL', 'JCOLL', 'VALEUR'
    jcoll = 1
    do kterm = 1, nz
!
!       --- PARTIE TRIANGULAIRE SUPERIEURE
        if (smdi(jcoll) .lt. kterm) jcoll = jcoll+1
        iligl = zi4(jsmhc-1+kterm)
        if (lmd) then
            iligg = zi(jnlogl+iligl-1)
            jcolg = zi(jnlogl+jcoll-1)
        else
            iligg = iligl
            jcolg = jcoll
        end if
        if ((.not. lsym) .and. (iligg .ge. jcolg)) then
            coltmp = jcolg
            jcolg = iligg
            iligg = coltmp
        end if
        if (ltypr) then
            write (ific, 1001) iligg, jcolg, zr(jvale-1+kterm)
        else
            write (ific, 1002) iligg, jcolg, dble(zc(jvale-1+kterm)), &
                dimag(zc(jvale-1+kterm))
        end if
!
!        --- PARTIE TRIANGULAIRE INFERIEURE
        if ((.not. lsym) .and. (iligg .ne. jcolg)) then
            if (ltypr) then
                write (ific, 1001) jcolg, iligg, zr(jval2-1+kterm)
            else
                write (ific, 1002) jcolg, iligg, dble(zc(jval2-1+kterm)), &
                    dimag(zc(jval2-1+kterm))
            end if
        end if
!
!       --- SI 'MATLAB' ET SYMETRIQUE , PSEUDO PARTIE INFERIEURE
        if (lsym .and. (typimp .eq. 'MATLAB') .and. (iligg .ne. jcolg)) then
            if (ltypr) then
                write (ific, 1001) jcolg, iligg, zr(jvale-1+kterm)
            else
                write (ific, 1002) jcolg, iligg, dble(zc(jvale-1+kterm)), &
                    dimag(zc(jvale-1+kterm))
            end if
        end if
!
    end do
!
!
!
!
!     -- IMPRESSION DES CARACTERISTIQUES DES EQUATIONS :
!     --------------------------------------------------
!
    if ((typimp .eq. ' ') .or. (typimp .eq. 'ASTER')) then
        write (ific, *) ' '
        write (ific, *) 'DESCRIPTION DES EQUATIONS :'
        write (ific, *) ' '
        write (ific, *) '   NUM_EQUA NOEUD    CMP'
        call jeveuo(nonu//'.NUME.DEEQ', 'L', vi=deeq)
        call jeveuo(nonu//'.NUME.REFN', 'L', vk24=refn)
        call jelira(nonu//'.NUME.DEEQ', 'LONMAX', n1)
        nomgd = refn(2) (1:8)
        call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jcmp)
        ASSERT(n1 .eq. 2*n)
        do k = 1, n
            nuno = deeq(2*(k-1)+1)
            nucmp = deeq(2*(k-1)+2)
            if (l_parallel_matrix) then
                if (zi(jprddl-1+k) .eq. rang) then
                    localOrGhost = 'DDL_LOCAL'
                else
                    localOrGhost = 'DDL_GHOST'
                end if
            end if
            if (nuno .gt. 0 .and. nucmp .gt. 0) then
                nono = int_to_char8(nuno)
                nocmp = zk8(jcmp-1+nucmp)
                write (ific, 1004) k, nono, nocmp, localOrGhost
            else if (nucmp .lt. 0) then
                ASSERT(nuno .gt. 0)
                nono = int_to_char8(nuno)
                nocmp = zk8(jcmp-1-nucmp)
                if (zi(jdelg-1+k) .eq. -1) then
                    write (ific, 1005) k, nono, nocmp, localOrGhost, ' LAGR1 BLOCAGE'
                else
                    ASSERT(zi(jdelg-1+k) .eq. -2)
                    write (ific, 1005) k, nono, nocmp, localOrGhost, ' LAGR2 BLOCAGE'
                end if
            else
                ASSERT(nuno .eq. 0 .and. nucmp .eq. 0)
                nono = ' '
                nocmp = ' '
                if (zi(jdelg-1+k) .eq. -1) then
                    write (ific, 1005) k, nono, nocmp, localOrGhost, ' LAGR1 RELATION LINEAIRE '
                else
                    ASSERT(zi(jdelg-1+k) .eq. -2)
                    write (ific, 1005) k, nono, nocmp, localOrGhost, ' LAGR2 RELATION LINEAIRE '
                end if
            end if
        end do
    end if
!
!     --- FIN IMPRESSION
    write (ific, *) '% --------------------------------------------'//&
     &              '----------------------------------------------'
!
!
!
1001 format(2i12, 1(1x, 1pe23.15))
1002 format(2i12, 2(1x, 1pe23.15, 1pe23.15))
1003 format(3a12, 1x, 1a14)
1004 format(i12, 2(1x, a8), a10)
1005 format(i12, 2(1x, a8), a10, a)
!
    call jedema()
end subroutine
