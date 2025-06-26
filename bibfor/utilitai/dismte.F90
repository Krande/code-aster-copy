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
subroutine dismte(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(TYPE_ELEM)
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismtm.h"
#include "asterfort/indiis.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/teattr.h"
!
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi
    character(len=*) :: nomobz, repkz
    character(len=32) :: repk
    character(len=16) :: nomob
! ----------------------------------------------------------------------
!    IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE TYPE TYPE_ELEM (K16)
!    OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    integer(kind=8) :: ibid
    character(len=8) :: nomtm
    character(len=16) :: nophen, nomodl
    integer(kind=8) :: ite, nbphen, nbtm, ico, iphen, nbmodl, imodl, iamodl, ii
    integer(kind=8) :: nbopt, iopt, ioptte, nrig, irig
    parameter(nrig=5)
    character(len=16) :: optrig(nrig)
    integer(kind=8), pointer :: optte(:) => null()
    data optrig/'RIGI_ACOU', 'RIGI_THER', 'RIGI_MECA', 'RIGI_MECA_TANG',&
     &     'FULL_MECA'/
!
!
!
    call jemarq()
    nomob = nomobz
    repk = ' '
    ierd = 0
!
    call jenonu(jexnom('&CATA.TE.NOMTE', nomob), ite)
    call jeveuo('&CATA.TE.TYPEMA', 'L', ibid)
    nomtm = zk8(ibid-1+ite)
!
!
    if (questi .eq. 'PHENOMENE') then
!     --------------------------------------
        call jelira('&CATA.PHENOMENE', 'NOMUTI', nbphen)
        call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbtm)
!
        ico = 0
        do iphen = 1, nbphen
            call jenuno(jexnum('&CATA.PHENOMENE', iphen), nophen)
            call jelira('&CATA.'//nophen, 'NMAXOC', nbmodl)
            do imodl = 1, nbmodl
                call jeveuo(jexnum('&CATA.'//nophen, imodl), 'L', iamodl)
                ii = indiis(zi(iamodl), ite, 1, nbtm)
                if (ii .gt. 0) then
                    repk = nophen
                    goto 30
!
                end if
            end do
        end do
        ierd = 1
30      continue
!
    else if (questi .eq. 'FORMULATION') then
!
        call teattr('C', 'FORMULATION', repk, ierd, typel=nomob)
        if (ierd .eq. 1) then
            repk = ' '
            ierd = 0
        end if
!
!
!
    else if (questi .eq. 'MODELISATION') then
!     --------------------------------------
        call jelira('&CATA.PHENOMENE', 'NOMUTI', nbphen)
        call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbtm)
!
        ico = 0
        do iphen = 1, nbphen
            call jenuno(jexnum('&CATA.PHENOMENE', iphen), nophen)
            call jelira('&CATA.'//nophen, 'NMAXOC', nbmodl)
            do imodl = 1, nbmodl
                call jeveuo(jexnum('&CATA.'//nophen, imodl), 'L', iamodl)
                ii = indiis(zi(iamodl), ite, 1, nbtm)
                if (ii .gt. 0) then
                    call jenuno(jexnum('&CATA.'//nophen(1:13)//'.MODL', imodl), nomodl)
                    repk = nomodl
                    ico = ico+1
                end if
            end do
            if (ico .gt. 0) then
                if (ico .eq. 1) then
                    goto 60
!
                else if (ico .gt. 1) then
                    repk = '#PLUSIEURS'
                    goto 60
!
                end if
            end if
        end do
        ierd = 1
60      continue
!
!
    else if ((questi .eq. 'PHEN_MODE')) then
!     --------------------------------------
        call jelira('&CATA.PHENOMENE', 'NOMUTI', nbphen)
        call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbtm)
!
        ico = 0
        do iphen = 1, nbphen
            call jenuno(jexnum('&CATA.PHENOMENE', iphen), nophen)
            call jelira('&CATA.'//nophen, 'NMAXOC', nbmodl)
            do imodl = 1, nbmodl
                call jeveuo(jexnum('&CATA.'//nophen, imodl), 'L', iamodl)
                ii = indiis(zi(iamodl), ite, 1, nbtm)
                if (ii .gt. 0) then
                    call jenuno(jexnum('&CATA.'//nophen(1:13)//'.MODL', imodl), nomodl)
                    repk = nophen//nomodl
                    ico = ico+1
                end if
            end do
        end do
        if (ico .gt. 1) then
            repk = '#PLUSIEURS'
        else if (ico .eq. 0) then
            repk = '#AUCUN'
        end if
!
!
    else if (questi .eq. 'NOM_TYPMAIL') then
!     --------------------------------------
        call jeveuo('&CATA.TE.TYPEMA', 'L', ibid)
        repk = zk8(ibid-1+ite)
!
!
    else if (questi .eq. 'TYPE_TYPMAIL') then
!     --------------------------------------
        call dismtm(questi, nomtm, repi, repk, ierd)
!
!
    else if (questi .eq. 'NBNO_TYPMAIL') then
!     --------------------------------------
        call dismtm(questi, nomtm, repi, repk, ierd)
!
!
    else if (questi .eq. 'DIM_TOPO') then
!     --------------------------------------
        call dismtm(questi, nomtm, repi, repk, ierd)
!
!
    else if (questi .eq. 'DIM_GEOM') then
!     --------------------------------------
        call jeveuo('&CATA.TE.DIM_GEOM', 'L', ibid)
        repi = zi(ibid-1+ite)
!
!
    else if (questi .eq. 'CALC_RIGI') then
!     --------------------------------------
        repk = 'NON'
        call jeveuo('&CATA.TE.OPTTE', 'L', vi=optte)
        call jelira('&CATA.OP.NOMOPT', 'NOMMAX', nbopt)
        do irig = 1, nrig
            call jenonu(jexnom('&CATA.OP.NOMOPT', optrig(irig)), iopt)
            ASSERT(iopt .gt. 0)
            ioptte = optte((ite-1)*nbopt+iopt)
            if (ioptte .eq. 0) goto 90
            repk = 'OUI'
            goto 100
!
90          continue
        end do
100     continue
!
!
    else
        ierd = 1
    end if
!
    repkz = repk
    call jedema()
end subroutine
