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

subroutine refe81(nomres, basmod, raid, mass, amor, &
                  mailla)
    implicit none
!  P. RICHARD     DATE 13/07/90
!-----------------------------------------------------------------------
!  BUT : < CREATION DU REFE ET DU DESC POUR OP0081 >
!
!        - RECUPERER LES NOMS UTILISATEUR DES CONCEPTS ASSOCIES AUX
!          MATRICES ASSEMBLEES ET BASE MODALE CONSIDEREES
!        - EFFECTUER QUELQUES CONTROLES ET DETERMINER
!          OPTION DE CALCUL MATRICES PROJETEES
!-----------------------------------------------------------------------
!
! NOMRES /I/ : NOM UTILISATEUR DU RESULTAT
! BASMOD /O/ : NOM UT DE LA BASE MODALE DE PROJECTION
! RAID   /O/ : NOM UT DE LA MATRICE RAIDEUR A PROJETER
! MASS   /O/ : NOM UT DE LA MATRICE DE MASSE A PROJETER
! AMOR   /O/ : NOM UT DE LA MATRICE D'AMORTISSEMENT A PROJETER
! MAILLA /O/ : NOM UT DU MAILLAGE EN AMONT
!
!
!
!
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=24) :: valk(2), typbas
    character(len=8) :: nomres, mailla, basmod, maillb, bl8, lintf
    character(len=14) :: numddl, numbis, numter
    character(len=19) :: raid, mass, amor
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iadref, ioc, iret, lddesc
    integer(kind=8) ::  nbval
    integer(kind=8), pointer :: idc_desc(:) => null()
    character(len=24), pointer :: idc_refe(:) => null()
!-----------------------------------------------------------------------
    data bl8/'        '/
!-----------------------------------------------------------------------
    call jemarq()
!
! --- RECUPERATION BASE MODALE OBLIGATOIRE
!
    call getvid(' ', 'BASE_MODALE', scal=basmod, nbret=nbval)
!
! --- RECUPERATION MATRICE RAIDEUR EN ARGUMENT
!
    raid = bl8
    call getvid(' ', 'MATR_RIGI', nbval=0, nbret=ioc)
    ioc = -ioc
    if (ioc .eq. 0) then
        call dismoi('REF_RIGI_PREM', basmod, 'RESU_DYNA', repk=raid, arret='C', &
                    ier=iret)
    else if (ioc .eq. 1) then
        call getvid(' ', 'MATR_RIGI', scal=raid, nbret=iret)
    else
        call utmess('F', 'ALGORITH14_14')
    end if
!
! --- RECUPERATION MATRICE MASSE EN ARGUMENT
!
    mass = bl8
    call getvid(' ', 'MATR_MASS', nbval=0, nbret=ioc)
    ioc = -ioc
    if (ioc .eq. 0) then
        call dismoi('REF_MASS_PREM', basmod, 'RESU_DYNA', repk=mass, arret='C', &
                    ier=iret)
    else if (ioc .eq. 1) then
        call getvid(' ', 'MATR_MASS', scal=mass, nbret=iret)
    else
        call utmess('F', 'ALGORITH14_15')
    end if
!
! --- RECUPERATION MATRICE AMORTISSEMENT EN ARGUMENT
!
    amor = bl8
    call getvid(' ', 'MATR_AMOR', nbval=0, nbret=ioc)
    ioc = -ioc
    if (ioc .eq. 0) then
        call dismoi('REF_AMOR_PREM', basmod, 'RESU_DYNA', repk=amor, arret='C', &
                    ier=iret)
    else if (ioc .eq. 1) then
        call getvid(' ', 'MATR_AMOR', scal=amor, nbret=iret)
    else
        call utmess('F', 'ALGORITH14_16')
    end if
!
    call dismoi('NUME_DDL', basmod, 'RESU_DYNA', repk=numddl, arret='C', &
                ier=iret)
    call dismoi('REF_INTD_PREM', basmod, 'RESU_DYNA', repk=lintf, arret='C', &
                ier=iret)
    call dismoi('TYPE_BASE', basmod, 'RESU_DYNA', repk=typbas, arret='C', &
                ier=iret)
!
!
! --- RECUPERATION MAILLAGE
!
    if (lintf .ne. bl8) then
        call jeveuo(lintf//'.IDC_REFE', 'L', vk24=idc_refe)
        mailla = idc_refe(1)
    else
        call dismoi('NOM_MAILLA', numddl, 'NUME_DDL', repk=mailla)
    end if
!
! --- TRAITEMENT COHERENCE MATRICE ASSEMBLEES
!     SI MASSF ET RAIDF NON BL8
!
    if (raid .eq. bl8 .or. mass .eq. bl8) goto 10
!
    call dismoi('NOM_NUME_DDL', raid, 'MATR_ASSE', repk=numddl)
    call dismoi('NOM_MAILLA', numddl, 'NUME_DDL', repk=maillb)
    call dismoi('NOM_NUME_DDL', mass, 'MATR_ASSE', repk=numbis)
!
    if (amor .ne. bl8) then
        call dismoi('NOM_NUME_DDL', amor, 'MATR_ASSE', repk=numter)
    end if
!
! --- CONTROLE DE LA COHERENCE DES MATRICES ASSEMBLEES
!
    if (numddl .ne. numbis) then
        valk(1) = mass
        valk(2) = raid
        call utmess('F', 'ALGORITH14_21', nk=2, valk=valk)
    end if
!
    if (amor .ne. bl8) then
        if (numddl .ne. numter) then
            valk(1) = amor
            valk(2) = raid
            call utmess('F', 'ALGORITH14_22', nk=2, valk=valk)
        end if
    end if
!
    if (mailla .ne. maillb) then
        valk(1) = maillb
        valk(2) = mailla
        call utmess('F', 'ALGORITH14_23', nk=2, valk=valk)
    end if
!
10  continue
!
! --- REMPLISSAGE DU .REFE
!
    call wkvect(nomres//'.MAEL_REFE', 'G V K24', 2, iadref)
    zk24(iadref) = basmod
    zk24(iadref+1) = mailla
!
! --- REMPLISSAGE DU .DESC
!
    call wkvect(nomres//'.MAEL_DESC', 'G V I', 3, lddesc)
    if (lintf .ne. bl8) then
        call jeveuo(lintf//'.IDC_DESC', 'L', vi=idc_desc)
        zi(lddesc) = idc_desc(2)
        zi(lddesc+1) = idc_desc(3)
        zi(lddesc+2) = idc_desc(4)
    end if
!
    call jedema()
end subroutine
