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
subroutine avdowh(nbvec, nbordr, nommat, nomcri, ncycl, &
                  jgdeq, grdvie, forvie, post, jdomel, &
                  jnrupt)
! person_in_charge: van-xuan.tran at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/limend.h"
#include "asterfort/rccome.h"
#include "asterfort/rcpare.h"
#include "asterfort/rcvale.h"
#include "asterfort/renrfa.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbvec, nbordr, ncycl(nbvec)
!    real(kind=8) :: gdeq(nbvec*nbordr)
!    real(kind=8) :: nrupt(nbvec*nbordr), domel(nbvec*nbordr)
    integer(kind=8) :: jgdeq, jnrupt, jdomel
    aster_logical :: post
    character(len=8) :: nommat
    character(len=16) :: nomcri, forvie, grdvie
! ----------------------------------------------------------------------
! BUT: CALCULER LE DOMMAGE ELEMENTAIRE A PARTIR D'UNE COURBE
!      GRANDEUR EQ - VIE POUR TOUS LES CYCLES
!      ELEMETAIRES DE CHAQUE VECTEUR NORMAL.
! ----------------------------------------------------------------------
! ARGUMENTS :
!  NBVEC    IN   I  : NOMBRE DE VECTEURS NORMAUX.
!  NBORDR   IN   I  : NOMBRE DE NUMEROS D'ORDRE.
!  NOMMAT   IN   K  : NOM DU MATERIAU.
!  NOMCRI   IN   K  : NOM DU CRITERE.
!  NCYCL    IN   I  : NOMBRE DE CYCLES ELEMENTAIRES POUR TOUS LES
!                     VECTEURS NORMAUX.
!  JGDEQ     IN   I  : ADDRESSE VECTEUR CONTENANT LES VALEURS DE LA GRANDEUR
!                     EQUIVALENTE (SIGEQ OU EPSEQ), POUR TOUS LES SOUS
!                     CYCLES DE CHAQUE VECTEUR NORMAL.
!  JDOMEL    OUT  I  : ADDRESSE VECTEUR CONTENANT LES VALEURS DES DOMMAGES
!                     ELEMENTAIRES, POUR TOUS LES SOUS CYCLES
!                     DE CHAQUE VECTEUR NORMAL.
!  JNRUPT    OUT  I  : ADDRESSE VECTEUR CONTENANT LES NOMBRES DE CYCLES
!                     ELEMENTAIRES, POUR TOUS LES SOUS CYCLES
!                     DE CHAQUE VECTEUR NORMAL.
! ----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: ivect, icycl, adrs, i
    integer(kind=8) :: icodre(1)
    character(len=16) :: kbid
    character(len=8) :: nomgrd
    aster_logical :: limit
!     ------------------------------------------------------------------
!
!234567                                                              012
!
    call jemarq()
!
! INITITIALISATION
    do i = 1, nbvec*nbordr
        zr(jdomel+i) = 0
    end do
!
    if (.not. post) then
        call rccome(nommat, 'FATIGUE', icodre(1))
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_24')
        end if
    end if
!
    if (nomcri(1:16) .eq. 'FATESOCI_MODI_AV') then
        call rcpare(nommat, 'FATIGUE', 'MANSON_COFFIN', icodre(1))
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_89', sk=nomcri(1:16))
        end if
!
        do ivect = 1, nbvec
            do icycl = 1, ncycl(ivect)
                adrs = (ivect-1)*nbordr+icycl
!
                call rcvale(nommat, 'FATIGUE', 1, 'EPSI    ', zr(jgdeq+adrs), &
                            1, 'MANSON_COFFIN', zr(jnrupt+adrs), icodre(1), 1)
!
                call limend(nommat, zr(jgdeq+adrs), 'MANSON_COFFIN', kbid, limit)
                if (limit) then
                    zr(jnrupt+adrs) = r8maem()
                else
                    call rcvale(nommat, 'FATIGUE', 1, 'EPSI    ', zr(jgdeq+adrs), &
                                1, 'MANSON_COFFIN', zr(jnrupt+adrs), icodre(1), 1)
                end if
!
                zr(jdomel+adrs) = 1.0d0/zr(jnrupt+adrs)
                zr(jnrupt+adrs) = nint(zr(jnrupt+adrs))
!
            end do
        end do
!
    elseif ((nomcri(1:14) .eq. 'MATAKE_MODI_AV') .or. (nomcri(1:16) &
                                                       .eq. 'DANG_VAN_MODI_AV')) then
        call rcpare(nommat, 'FATIGUE', 'WOHLER', icodre(1))
        if (icodre(1) .eq. 1) then
            call utmess('F', 'FATIGUE1_90', sk=nomcri(1:16))
        end if
!
        do ivect = 1, nbvec
            do icycl = 1, ncycl(ivect)
                adrs = (ivect-1)*nbordr+icycl
!
                call limend(nommat, zr(jgdeq+adrs), 'WOHLER', kbid, limit)
                if (limit) then
                    zr(jnrupt+adrs) = r8maem()
                else
                    call rcvale(nommat, 'FATIGUE', 1, 'SIGM    ', zr(jgdeq+adrs), &
                                1, 'WOHLER  ', zr(jnrupt+adrs), icodre(1), 1)
                end if
!
                zr(jdomel+adrs) = 1.0d0/zr(jnrupt+adrs)
                zr(jnrupt+adrs) = nint(zr(jnrupt+adrs))
!
            end do
        end do
!
    else if (nomcri(1:7) .eq. 'FORMULE') then
!
        do ivect = 1, nbvec
            do icycl = 1, ncycl(ivect)
                adrs = (ivect-1)*nbordr+icycl
!
                call limend(nommat, zr(jgdeq+adrs), grdvie, forvie, limit)
!
                if (limit) then
                    zr(jnrupt+adrs) = r8maem()
                else
!
                    if (grdvie .eq. 'WOHLER') then
                        nomgrd = 'SIGM    '
!
                        call rcvale(nommat, 'FATIGUE', 1, nomgrd, zr(jgdeq+adrs), &
                                    1, grdvie, zr(jnrupt+adrs), icodre(1), 1)
                    end if
!
                    if (grdvie .eq. 'MANSON_COFFIN') then
                        nomgrd = 'EPSI    '
                        call rcvale(nommat, 'FATIGUE', 1, nomgrd, zr(jgdeq+adrs), &
                                    1, grdvie, zr(jnrupt+adrs), icodre(1), 1)
!
                    end if
!
                    if (grdvie .eq. 'FORM_VIE') then
                        call renrfa(forvie, zr(jgdeq+adrs), zr(jnrupt+adrs), icodre(1))
                    end if
!
                    zr(jdomel+adrs) = 1.0d0/zr(jnrupt+adrs)
                    zr(jnrupt+adrs) = nint(zr(jnrupt+adrs))
!
                end if
!
!
            end do
        end do
!
!
!
    end if
!
    call jedema()
!
end subroutine
