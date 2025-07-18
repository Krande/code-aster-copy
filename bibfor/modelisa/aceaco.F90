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
subroutine aceaco(nomu, noma, lmax, locagb, locamb, &
                  nbocc)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8rddg.h"
#include "asterfort/alcart.h"
#include "asterfort/angvx.h"
#include "asterfort/assert.h"
#include "asterfort/exisd.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: lmax, nbocc
    aster_logical :: locagb, locamb
    character(len=8) :: nomu, noma
!                          AFFE_CARA_ELEM
!
!     AFFECTATION DES CARACTERISTIQUES POUR L'ELEMENT COQUE
!
! ----------------------------------------------------------------------
!  IN
!     NOMU   : NOM UTILISATEUR DE LA COMMANDE
!     NOMA   : NOM DU MAILLAGE
!     LMAX   : NOMBRE MAXIMUM DE MAILLE AFFECTEES PAR AFFE_CARA_ELEM
!     LOCAGB : EXISTANCE DE GRILLE
!     LOCAMB : EXISTANCE DE MEMBRANE
!     NBOCC  : NOMBRE D'OCCURENCES DU MOT CLE COQUE
! ----------------------------------------------------------------------
    integer(kind=8) :: nvec, i, ioc, jdcc, jdls, jdvc, jdccf, jdvcf
    integer(kind=8) :: na, nco, ncr, nex, ng, nin, nk, nv, nvf, nexf
    integer(kind=8) :: iret
    aster_logical :: lcartf
    real(kind=8) :: ang(2), epa, kappa, correc, rigi, excent
    real(kind=8) :: vect(3), xiner
    character(len=8) :: inert, korrec, epaf, excf
    character(len=19) :: cartco, cartcf
    character(len=24) :: tmpnco, tmpvco, tmpncf, tmpvcf
!-----------------------------------------------------------------------
!
! --- CONSTRUCTION DES CARTES ET ALLOCATION
    call jemarq()
!
!     CARTE POUR LES VALEURS REELLES
    cartco = nomu//'.CARCOQUE'
    call exisd('CARTE', cartco, iret)
    if (iret .eq. 0) then
        call alcart('G', cartco, noma, 'CACOQU_R')
    end if
    tmpnco = cartco//'.NCMP'
    tmpvco = cartco//'.VALV'
    call jeveuo(tmpnco, 'E', jdcc)
    call jeveuo(tmpvco, 'E', jdvc)
!     LES NOMS DES GRANDEURS REELLES
! PLAQUE : 'EP','ALPHA','BETA','KAPPA','CTOR', 'EXCENT','INERTIE'
! COQUE  : 'EP','ALPHA','BETA',       ,'CTOR', 'EXCENT','INERTIE'
    zk8(jdcc) = 'EP'
    zk8(jdcc+1) = 'ALPHA'
    zk8(jdcc+2) = 'BETA'
    zk8(jdcc+3) = 'KAPPA'
    zk8(jdcc+4) = 'C_METR'
    zk8(jdcc+5) = 'CTOR'
    zk8(jdcc+6) = 'EXCENT'
    zk8(jdcc+7) = 'INERTIE'
!
!     CARTE POUR LES FONCTIONS
    cartcf = nomu//'.CARCOQUF'
    call exisd('CARTE', cartcf, iret)
    lcartf = .false.
    if (iret .eq. 0) then
! ------ DOIT-ON CREER LA CARTE DE FONCTION
        do ioc = 1, nbocc
            call getvid('COQUE', 'EPAIS_FO', iocc=ioc, scal=epaf, nbret=nvf)
            call getvid('COQUE', 'EXCENTREMENT_FO', iocc=ioc, scal=excf, nbret=nexf)
            if (nvf+nexf .ne. 0) then
                lcartf = .true.
                goto 110
            end if
        end do
110     continue
!
!        CARTE POUR LES NOMS DES FONCTIONS
        if (lcartf) then
            call alcart('V', cartcf, noma, 'CACOQU_F')
        end if
    else
        lcartf = .true.
    end if
!     SI LA CARTE EXISTE
    if (lcartf) then
        tmpncf = cartcf//'.NCMP'
        tmpvcf = cartcf//'.VALV'
        call jeveuo(tmpncf, 'E', jdccf)
        call jeveuo(tmpvcf, 'E', jdvcf)
!        LES NOMS DES FONCTIONS
        zk8(jdccf) = 'EP'
        zk8(jdccf+1) = 'EXCENT'
    end if
!
    call wkvect('&&TMPCOQUE', 'V V K24', lmax, jdls)
!
! --- LECTURE DES VALEURS ET AFFECTATION DANS : CARTCO OU CARTCF
    do ioc = 1, nbocc
        ang(1) = 0.d0
        ang(2) = 0.d0
        correc = 0.d0
        kappa = 0.d0
        excent = 0.d0
        xiner = 0.d0
        inert = 'NON'
        call getvem(noma, 'GROUP_MA', 'COQUE', 'GROUP_MA', ioc, &
                    lmax, zk24(jdls), ng)
        call getvr8('COQUE', 'EPAIS', iocc=ioc, scal=epa, nbret=nv)
        call getvid('COQUE', 'EPAIS_FO', iocc=ioc, scal=epaf, nbret=nvf)
        call getvr8('COQUE', 'ANGL_REP', iocc=ioc, nbval=2, vect=ang, &
                    nbret=na)
        call getvr8('COQUE', 'VECTEUR', iocc=ioc, nbval=3, vect=vect, &
                    nbret=nvec)
        call getvr8('COQUE', 'A_CIS', iocc=ioc, scal=kappa, nbret=nk)
        call getvtx('COQUE', 'MODI_METRIQUE', iocc=ioc, scal=korrec, nbret=nco)
        call getvr8('COQUE', 'COEF_RIGI_DRZ', iocc=ioc, scal=rigi, nbret=ncr)
        call getvr8('COQUE', 'EXCENTREMENT', iocc=ioc, scal=excent, nbret=nex)
        call getvid('COQUE', 'EXCENTREMENT_FO', iocc=ioc, scal=excf, nbret=nexf)
        call getvtx('COQUE', 'INER_ROTA', iocc=ioc, scal=inert, nbret=nin)
!        EPAIS EST OBLIGATOIRE : ASSERT SI PAS LA
        if (nv .ne. 0) then
            zr(jdvc) = epa
            if (lcartf) zk8(jdvcf) = '&&ACEACO'
        else if (nvf .ne. 0) then
            zr(jdvc) = 0.0d0
            if (lcartf) zk8(jdvcf) = epaf
        else
            ASSERT(.false.)
        end if
        zr(jdvc+1) = ang(1)
        zr(jdvc+2) = ang(2)
        if (nvec .ne. 0) then
            call angvx(vect, ang(1), ang(2))
            zr(jdvc+1) = ang(1)*r8rddg()
            zr(jdvc+2) = ang(2)*r8rddg()
        end if
        zr(jdvc+3) = kappa
        if (korrec .eq. 'OUI') correc = 1.d0
        zr(jdvc+4) = correc
        zr(jdvc+5) = rigi
!
        zr(jdvc+6) = excent
        if (lcartf) zk8(jdvcf+1) = '&&ACEACO'
        if (nexf .ne. 0) then
            zr(jdvc+6) = 0.0d0
            if (lcartf) zk8(jdvcf+1) = excf
        end if
!
        if ((nex+nexf) .ne. 0 .and. nin .eq. 0) inert = 'OUI'
        if (inert .eq. 'OUI') xiner = 1.d0
        zr(jdvc+7) = xiner
!
! ---    "GROUP_MA" = TOUTES LES MAILLES DE LA LISTE DE GROUPES MAILLES
        if (ng .gt. 0) then
            do i = 1, ng
                call nocart(cartco, 2, 8, groupma=zk24(jdls+i-1))
            end do
            if (lcartf) then
                do i = 1, ng
                    call nocart(cartcf, 2, 2, groupma=zk24(jdls+i-1))
                end do
            end if
        end if
!
    end do
!
    call jedetr('&&TMPCOQUE')
!     SI NI GRILLE NI MEMBRANE
    if ((.not. locagb) .and. (.not. locamb)) then
        call jedetr(tmpnco)
        call jedetr(tmpvco)
        if (lcartf) then
            call jedetr(tmpncf)
            call jedetr(tmpvcf)
        end if
    end if
!
    call jedema()
end subroutine
