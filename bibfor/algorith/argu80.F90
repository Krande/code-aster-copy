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

subroutine argu80(nomres)
    implicit none
!
!***********************************************************************
!    P. RICHARD     DATE 28/03/91
!-----------------------------------------------------------------------
!  BUT : RECUPERER LES ARGUMENTS D'APPEL (SAUF LES DIAMETRES ET LE
!        NOMBRE DE MODES A CALCULER) ET CREATION DES OBJETS
!        CORRESPONDANTS
!        VERIFICATION DES PROPRIETES DE REPETITIVITE SUR LE MAILLAGE
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM UTILISATEUR DU CONCEPT RESULTAT
!
!
!
!
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/verecy.h"
#include "asterfort/wkvect.h"
!
    character(len=24) :: valk
    character(len=8) :: droite, gauche, axe, typd, typg, typa
    character(len=8) :: nomres, intf
    character(len=8) :: blanc
    character(len=8) :: kar8
    integer(kind=8) :: vali
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ibaxe, ibid, lddnbs, lddnin, lddtyp
    integer(kind=8) :: nbsec, ndist, numa, numd, numg, nveri
    real(kind=8) :: dist, prec
    character(len=24), pointer :: cycl_refe(:) => null()
    character(len=8), pointer :: idc_type(:) => null()
!-----------------------------------------------------------------------
    data blanc/'  '/
!-----------------------------------------------------------------------
!
!
!-------------CREATION DES OBJETS DE LA SDD RESULTAT--------------------
!
    call jemarq()
    call wkvect(nomres//'.CYCL_NUIN', 'G V I', 3, lddnin)
    call wkvect(nomres//'.CYCL_TYPE', 'G V K8', 1, lddtyp)
    call wkvect(nomres//'.CYCL_NBSC', 'G V I', 1, lddnbs)
!
!--------------------RECUPERATION DES CONCEPTS AMONTS-------------------
!
    call jeveuo(nomres//'.CYCL_REFE', 'L', vk24=cycl_refe)
    intf = cycl_refe(2) (1:8)
!
!----------RECUPERATION NOM DES INTERFACES DE LIAISON-------------------
!
    call getvtx('LIAISON', 'DROITE', iocc=1, scal=kar8, nbret=ibid)
    droite = kar8
    call getvtx('LIAISON', 'GAUCHE', iocc=1, scal=kar8, nbret=ibid)
    gauche = kar8
    call getvtx('LIAISON', 'AXE', iocc=1, nbval=0, nbret=ibaxe)
    if (ibaxe .eq. -1) then
        call getvtx('LIAISON', 'AXE', iocc=1, scal=kar8, nbret=ibid)
        axe = kar8
    else
        axe = ' '
    end if
!
!   RECUPERATION DES NUMEROS D'INTERFACE
!
!   INTERFACE DE DROITE OBLIGATOIRE
!
    call jenonu(jexnom(intf//'.IDC_NOMS', droite), numd)
    if (numd .eq. 0) then
        valk = droite
        call utmess('F', 'ALGORITH15_85', sk=valk)
    end if
!
!   INTERFACE DE GAUCHE OBLIGATOIRE
!
    call jenonu(jexnom(intf//'.IDC_NOMS', gauche), numg)
    if (numg .eq. 0) then
        valk = gauche
        call utmess('F', 'ALGORITH15_86', sk=valk)
    end if
!
!   INTERFACE AXE FACULTATIVE
!
    if (axe .ne. '        ') then
        call jenonu(jexnom(intf//'.IDC_NOMS', axe), numa)
        if (numa .eq. 0) then
            valk = axe
            call utmess('F', 'ALGORITH15_87', sk=valk)
        end if
    else
        numa = 0
    end if
!
    zi(lddnin) = numd
    zi(lddnin+1) = numg
    zi(lddnin+2) = numa
!
!   RECUPERATION DES TYPES DES INTERFACES
!
    call jeveuo(intf//'.IDC_TYPE', 'L', vk8=idc_type)
    typd = idc_type(numd)
    typg = idc_type(numg)
    if (numa .gt. 0) then
        typa = idc_type(numa)
    else
        typa = typd
    end if
!
!  VERIFICATIONS SUR LES TYPES INTERFACES
!
    if (typg .ne. typd .or. typa .ne. typd) then
        call utmess('F', 'ALGORITH15_88')
    end if
!
    if (typd .ne. 'MNEAL   ' .and. typd .ne. 'CRAIGB  ') then
        if (typd .ne. 'AUCUN   ' .and. typd .ne. 'CB_HARMO') then
            valk = typd
            call utmess('F', 'ALGORITH15_89', sk=valk)
        end if
    end if
!
! STOCKAGE TYPE INTERFACE
!
    zk8(lddtyp) = typd
!
!  RECUPERATION DU NOMBRE DE SECTEURS
!
    call getvis(blanc, 'NB_SECTEUR', iocc=1, scal=nbsec, nbret=ibid)
    if (nbsec .lt. 2) then
        vali = nbsec
        call utmess('F', 'ALGORITH15_59', si=vali)
    end if
!
    zi(lddnbs) = nbsec
!
!---------------VERIFICATION DE LA REPETITIVITE SUR MAILLAGE------------
!
    call getfac('VERI_CYCL', nveri)
    call getvr8('VERI_CYCL', 'PRECISION', iocc=1, scal=prec, nbret=ibid)
    call getvr8('VERI_CYCL', 'DIST_REFE', iocc=1, scal=dist, nbret=ndist)
    if (nveri .eq. 0) prec = 1.d-3
    if (ndist .eq. 0) then
!     --- AU CAS OU LA DISTANCE DE REFERENCE N'EST PAS DONNEE,ON DEVRAIT
!         LA LIRE DANS LA SD MAILLAGE (VOIR COMMANDE LIRE_MAILLAGE).
!         CE TRAVAIL N'ETANT PAS ACCOMPLI, ON MET DIST < 0 AFIN DE
!         SIGNIFIER A VERECY DE TRAVAILLER COMME AVANT
        dist = -1.d0
    end if
    call verecy(intf, numd, numg, nbsec, prec, &
                dist)
!
    call jedema()
end subroutine
