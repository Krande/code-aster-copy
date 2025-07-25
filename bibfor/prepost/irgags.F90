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
subroutine irgags(ncmpmx, nomcmp, nomsym, nbchs, nomchs, &
                  nbcmps, nomgds, ipcmps)
    implicit none
!
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
    integer(kind=8) :: ncmpmx, nbchs, nbcmps(*), ipcmps(*)
    character(len=*) :: nomcmp(*), nomsym, nomchs(*), nomgds(*)
!
!       RECHERCHE DE LA GRANDEUR SUPERTAB REPRESENTANT UNE GRANDEUR ASTER
!      ENTREE:
!         NCMPMX: NOMBRE MAXI DE CMP DE LA GRANDEUR NOMGD
!         NOMCMP: NOMS DES CMP
!         NOMSYM: NOM SYMBOLIQUE
!      SORTIE:
!         NBCHS : NOMBRE DE GRANDEURS SUPERTAB (=1)
!         NOMCHS: NOM DE L'INFORMATION SUPERTAB
!         NBCMPS: NOMBRE DE COMPOSANTES DE CHAQUE GRANDEUR SUPERTAB
!         NOMGDS: NOM DES GRANDEURS SUPERTAB
!         IPCMPS: POSITION DES COMPOSANTES SUPERTAB DANS LA COMPOSANTE
!                 ASTER
!-----------------------------------------------------------------------
    integer(kind=8) :: nbdepl, nbtemp, nbvari, nbsigm, nbepsm, nbflu, nbpres
    integer(kind=8) :: i, icmas, ient, ires, iva, ivar, j
    integer(kind=8), pointer :: videpl(:) => null()
    integer(kind=8), pointer :: viepsm(:) => null()
    integer(kind=8), pointer :: vitemp(:) => null()
    integer(kind=8), pointer :: viflu(:) => null()
    integer(kind=8), pointer :: visigm(:) => null()
    integer(kind=8), pointer :: vipres(:) => null()
    integer(kind=8), pointer :: vivari(:) => null()
    AS_ALLOCATE(vi=videpl, size=6)
    AS_ALLOCATE(vi=vitemp, size=1)
    AS_ALLOCATE(vi=visigm, size=6)
    AS_ALLOCATE(vi=viepsm, size=6)
    AS_ALLOCATE(vi=viflu, size=3)
    AS_ALLOCATE(vi=vipres, size=1)
    AS_ALLOCATE(vi=vivari, size=ncmpmx)
!
!  --- INITIALISATIONS ----
!
    nbdepl = 0
    nbflu = 0
    nbtemp = 0
    nbepsm = 0
    nbsigm = 0
    nbpres = 0
    nbvari = 0
!
    iva = 0
    nbchs = 0
!
!  --- RECHERCHE DES GRANDEURS SUPERTAB ASSOCIEES A LA GRANDEUR ASTER---
!
    do icmas = 1, ncmpmx
        if (nomcmp(icmas) .eq. 'DX') then
            nbdepl = nbdepl+1
            videpl(1) = icmas
        else if (nomcmp(icmas) .eq. 'DY') then
            nbdepl = nbdepl+1
            videpl(2) = icmas
        else if (nomcmp(icmas) .eq. 'DZ') then
            nbdepl = nbdepl+1
            videpl(3) = icmas
        else if (nomcmp(icmas) .eq. 'DRX') then
            nbdepl = nbdepl+1
            videpl(4) = icmas
        else if (nomcmp(icmas) .eq. 'DRY') then
            nbdepl = nbdepl+1
            videpl(5) = icmas
        else if (nomcmp(icmas) .eq. 'DRZ') then
            nbdepl = nbdepl+1
            videpl(6) = icmas
!
        else if (nomcmp(icmas) .eq. 'FLUX') then
            nbflu = nbflu+1
            viflu(1) = icmas
        else if (nomcmp(icmas) .eq. 'FLUY') then
            nbflu = nbflu+1
            viflu(2) = icmas
        else if (nomcmp(icmas) .eq. 'FLUZ') then
            nbflu = nbflu+1
            viflu(3) = icmas
!
        else if (nomcmp(icmas) .eq. 'TEMP') then
            nbtemp = nbtemp+1
            vitemp(1) = icmas
!
        else if (nomcmp(icmas) .eq. 'PRES') then
            nbpres = nbpres+1
            vipres(1) = icmas
!
        else if (nomcmp(icmas) .eq. 'SIXX') then
            nbsigm = nbsigm+1
            visigm(1) = icmas
        else if (nomcmp(icmas) .eq. 'SIXY') then
            nbsigm = nbsigm+1
            visigm(2) = icmas
        else if (nomcmp(icmas) .eq. 'SIYY') then
            nbsigm = nbsigm+1
            visigm(3) = icmas
        else if (nomcmp(icmas) .eq. 'SIXZ') then
            nbsigm = nbsigm+1
            visigm(4) = icmas
        else if (nomcmp(icmas) .eq. 'SIYZ') then
            nbsigm = nbsigm+1
            visigm(5) = icmas
        else if (nomcmp(icmas) .eq. 'SIZZ') then
            nbsigm = nbsigm+1
            visigm(6) = icmas
!
!
        else if (nomcmp(icmas) .eq. 'EPXX') then
            nbepsm = nbepsm+1
            viepsm(1) = icmas
        else if (nomcmp(icmas) .eq. 'EPXY') then
            nbepsm = nbepsm+1
            viepsm(2) = icmas
        else if (nomcmp(icmas) .eq. 'EPYY') then
            nbepsm = nbepsm+1
            viepsm(3) = icmas
        else if (nomcmp(icmas) .eq. 'EPXZ') then
            nbepsm = nbepsm+1
            viepsm(4) = icmas
        else if (nomcmp(icmas) .eq. 'EPYZ') then
            nbepsm = nbepsm+1
            viepsm(5) = icmas
        else if (nomcmp(icmas) .eq. 'EPZZ') then
            nbepsm = nbepsm+1
            viepsm(6) = icmas
!
        else
            nbvari = nbvari+1
            vivari(nbvari) = icmas
        end if
    end do
!
!  --- RECHERCHE DES GRANDEURS SUPERTAB ASSOCIEES A LA GRANDEUR ASTER---
!
    if (nbdepl .ne. 0) then
        if (nomsym .eq. 'DEPL' .or. nomsym .eq. 'VITE' .or. nomsym .eq. 'ACCE') then
            nbchs = nbchs+1
            nomchs(nbchs) = 'DEPL'
            nbcmps(nbchs) = nbdepl
            nomgds(nbchs) = 'DEPL'
            if (nomsym .eq. 'DEPL') nomgds(nbchs) = 'DEPL'
            if (nomsym .eq. 'VITE') nomgds(nbchs) = 'VITE'
            if (nomsym .eq. 'ACCE') nomgds(nbchs) = 'ACCE'
        end if
    end if
    if (nbflu .ne. 0) then
        nbchs = nbchs+1
        nomchs(nbchs) = 'FLUX'
        nbcmps(nbchs) = nbflu
        nomgds(nbchs) = 'FLUX'
    end if
    if (nbtemp .ne. 0) then
        nbchs = nbchs+1
        nomchs(nbchs) = 'TEMP'
        nbcmps(nbchs) = nbtemp
        nomgds(nbchs) = 'TEMP'
    end if
    if (nbsigm .ne. 0) then
        nbchs = nbchs+1
        nomchs(nbchs) = 'SIGM'
        nbcmps(nbchs) = nbsigm
        nomgds(nbchs) = 'SIGM'
    end if
    if (nbepsm .ne. 0) then
        nbchs = nbchs+1
        nomchs(nbchs) = 'EPSI'
        nbcmps(nbchs) = nbepsm
        nomgds(nbchs) = 'EPSI'
    end if
    if (nbpres .ne. 0) then
        nbchs = nbchs+1
        nomchs(nbchs) = 'PRES'
        nbcmps(nbchs) = nbpres
        nomgds(nbchs) = 'PRES'
    end if
    if (nbvari .ne. 0) then
        ient = nbvari/6
        ires = nbvari-(ient*6)
        if (ient .ne. 0) then
            do ivar = 1, ient
                nomchs(nbchs+ivar) = 'VARI'
                nbcmps(nbchs+ivar) = 6
                nomgds(nbchs+ivar) = 'VARI'
            end do
        end if
        if (ires .ne. 0) then
            nomchs(nbchs+ient+1) = 'VARI'
            nbcmps(nbchs+ient+1) = ires
            nomgds(nbchs+ient+1) = 'VARI'
            nbchs = nbchs+ient+1
        else
            nbchs = nbchs+ient
        end if
    end if
!
!  ---- POSITIONS DES COMPOSANTES SUPERTAB DANS LA GRANDEUR ASTER ----
!
!   -- Je ne pense pas que ce soit une bonne idee d'eclater
!      une grandeur ASTER en plusieurs grandeurs IDEAS.
!      DEPL_R et SIEF_R contiennent des composantes de meme nom et cela pose
!      probleme. Voir issue24844.
!      => on suppose que les cas courants, la premiere grandeur est la bonne
!      => nbchs=1
!      Puis, on verifie que ipcmps(k) > 0 (c'est un indice dans des tableaux)
    nbchs = 1
    do i = 1, nbchs
        if (nomchs(i) .eq. 'DEPL') then
            do j = 1, nbcmps(i)
                ipcmps((i-1)*ncmpmx+j) = videpl(j)
            end do
        else if (nomchs(i) .eq. 'PRES') then
            do j = 1, nbcmps(i)
                ipcmps((i-1)*ncmpmx+j) = vipres(j)
            end do
        else if (nomchs(i) .eq. 'FLUX') then
            do j = 1, nbcmps(i)
                ipcmps((i-1)*ncmpmx+j) = viflu(j)
            end do
        else if (nomchs(i) .eq. 'TEMP') then
            do j = 1, nbcmps(i)
                ipcmps((i-1)*ncmpmx+j) = vitemp(j)
            end do
        else if (nomchs(i) .eq. 'SIGM') then
            do j = 1, nbcmps(i)
                ipcmps((i-1)*ncmpmx+j) = visigm(j)
            end do
        else if (nomchs(i) .eq. 'EPSI') then
            do j = 1, nbcmps(i)
                ipcmps((i-1)*ncmpmx+j) = viepsm(j)
            end do
        else if (nomchs(i) .eq. 'VARI') then
            iva = iva+1
            do j = 1, nbcmps(i)
                ipcmps((i-1)*ncmpmx+j) = vivari((iva-1)*6+j)
            end do
        end if
!
        do j = 1, nbcmps(i)
            ASSERT(ipcmps((i-1)*ncmpmx+j) .gt. 0)
        end do
!
    end do
!
    AS_DEALLOCATE(vi=videpl)
    AS_DEALLOCATE(vi=visigm)
    AS_DEALLOCATE(vi=viepsm)
    AS_DEALLOCATE(vi=vitemp)
    AS_DEALLOCATE(vi=vipres)
    AS_DEALLOCATE(vi=viflu)
    AS_DEALLOCATE(vi=vivari)
end subroutine
