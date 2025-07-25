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
subroutine arg126(nomres)
    implicit none
!  P. RICHARD     DATE 13/10/93
!-----------------------------------------------------------------------
!  BUT:      < RECUPERATION DES ARGUMENTS POUR OP0126 >
!
!  RECUPERER LES ARGUMENTS UTILISATEUR POUR LA DEFINITION DU MODELE
!  GENERALISE. DEFINITION DES SOUS-STRUCTURES ET DES LIAISONS ENTRE
!  LES SOUS-STRUCTURES.
!
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM UTILISATEUR DU RESULTAT
!
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8dgrd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mgutdm.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: vali(3)
!
!
!
    character(len=8) :: nomres, lintf, nomsst, mclcou, nomcou
    character(len=16) :: clesst, clenom, clerot, clemcl, cletra, clelia, clel(4)
    character(len=24) :: repsst, nommcl, rotsst, famli, trasst
    character(len=24) :: valk(3)
    integer(kind=8) :: nbsst, i, j, ioc, ibid, ldnmcl, ldrot, nblia, ldlid, iret, ldtra
    real(kind=8) :: pi
!
!-----------------------------------------------------------------------
    data clesst, clenom/'SOUS_STRUC', 'NOM'/
    data clerot, clemcl/'ANGL_NAUT', 'MACR_ELEM_DYNA'/
    data clelia, cletra/'LIAISON', 'TRANS'/
    data clel/'SOUS_STRUC_1', 'SOUS_STRUC_2', 'INTERFACE_1',&
     &                    'INTERFACE_2'/
!-----------------------------------------------------------------------
!
    call jemarq()
    pi = 4.d+00*atan2(1.d+00, 1.d+00)
!
!-----TRAITEMENT DEFINITION SOUS-STRUCTURES-----------------------------
!
    call getfac(clesst, nbsst)
!
    if (nbsst .lt. 2) then
        vali(1) = 2
        vali(2) = nbsst
        call utmess('F', 'ALGORITH11_92', ni=2, vali=vali)
    end if
!
    repsst = nomres//'      .MODG.SSNO'
    nommcl = nomres//'      .MODG.SSME'
    rotsst = nomres//'      .MODG.SSOR'
    trasst = nomres//'      .MODG.SSTR'
!
    call jecreo(repsst, 'G N K8')
    call jeecra(repsst, 'NOMMAX', nbsst)
    call jecrec(nommcl, 'G V K8', 'NU', 'CONTIG', 'CONSTANT', &
                nbsst)
    call jecrec(rotsst, 'G V R', 'NU', 'CONTIG', 'CONSTANT', &
                nbsst)
    call jecrec(trasst, 'G V R', 'NU', 'CONTIG', 'CONSTANT', &
                nbsst)
    do i = 1, nbsst
        call jecroc(jexnum(nommcl, i))
        call jecroc(jexnum(rotsst, i))
        call jecroc(jexnum(trasst, i))
    end do
    call jeecra(nommcl, 'LONT', nbsst)
    call jeecra(rotsst, 'LONT', 3*nbsst)
    call jeecra(trasst, 'LONT', 3*nbsst)
!
!
!-----BOUCLE SUR LES SOUS-STRUCTURES DEFINIES-------------------------
!
    do i = 1, nbsst
        call getvtx(clesst, clenom, iocc=i, nbval=0, nbret=ioc)
        ioc = -ioc
        if (ioc .ne. 1) then
            vali(1) = 1
            vali(2) = ioc
            call utmess('F', 'ALGORITH11_93', ni=2, vali=vali)
        else
            call getvtx(clesst, clenom, iocc=i, scal=nomsst, nbret=ibid)
        end if
        call jecroc(jexnom(repsst, nomsst))
!
        call getvid(clesst, clemcl, iocc=i, nbval=0, nbret=ioc)
        ioc = -ioc
        if (ioc .ne. 1) then
            valk(1) = nomsst
            vali(1) = ioc
            vali(2) = 1
            call utmess('F', 'ALGORITH11_94', sk=valk(1), ni=2, vali=vali)
        else
            call getvid(clesst, clemcl, iocc=i, scal=mclcou, nbret=ibid)
        end if
        call jenonu(jexnom(repsst, nomsst), ibid)
        call jeveuo(jexnum(nommcl, ibid), 'E', ldnmcl)
        zk8(ldnmcl) = mclcou
!
!  TRAITEMENT DES ROTATIONS
!
        call jenonu(jexnom(repsst, nomsst), ibid)
        call jeveuo(jexnum(rotsst, ibid), 'E', ldrot)
        call getvr8(clesst, clerot, iocc=i, nbval=0, nbret=ioc)
        ioc = -ioc
        if (ioc .eq. 0) then
            do j = 1, 3
                zr(ldrot+j-1) = 0.d+00
            end do
        else if (ioc .eq. 3) then
            call getvr8(clesst, clerot, iocc=i, nbval=3, vect=zr(ldrot), &
                        nbret=ibid)
            do j = 1, 3
                zr(ldrot+j-1) = zr(ldrot+j-1)*r8dgrd()
            end do
        else
            valk(1) = nomsst
            vali(1) = ioc
            vali(2) = 3
            call utmess('F', 'ALGORITH11_95', sk=valk(1), ni=2, vali=vali)
        end if
!
!  TRAITEMENT DES TRANSLATIONS SI INTRODUIT PAR L'UTILISATEUR
!
!
        call jenonu(jexnom(repsst, nomsst), ibid)
        call jeveuo(jexnum(trasst, ibid), 'E', ldtra)
        call getvr8(clesst, cletra, iocc=i, nbval=0, nbret=ioc)
        ioc = -ioc
        if (ioc .eq. 0) then
            do j = 1, 3
                zr(ldtra+j-1) = 0.d+00
            end do
        else if (ioc .eq. 3) then
            call getvr8(clesst, cletra, iocc=i, nbval=3, vect=zr(ldtra), &
                        nbret=ibid)
        else
            valk(1) = nomsst
            vali(1) = ioc
            vali(2) = 3
            call utmess('F', 'ALGORITH11_96', sk=valk(1), ni=2, vali=vali)
        end if
!
!
    end do
!
!-----RECUPERATION DU NOMBRE DE LIAISONS DEFINIES-----------------------
!
    call getfac(clelia, nblia)
    if (nblia .eq. 0) then
        vali(1) = nblia
        vali(2) = 1
        call utmess('F', 'ALGORITH11_97', ni=2, vali=vali)
    end if
!
    famli = nomres//'      .MODG.LIDF'
    call jecrec(famli, 'G V K8', 'NU', 'DISPERSE', 'CONSTANT', &
                nblia)
    call jeecra(famli, 'LONMAX', 5)
!
!-----BOUCLE SUR LES LIAISONS------------------------------------------
!
    do i = 1, nblia
        call jecroc(jexnum(famli, i))
        call jeveuo(jexnum(famli, i), 'E', ldlid)
!
!  BOUCLE SUR LES SOUS-STRUCTURES DE LA LIAISON
!
        do j = 1, 2
            call getvtx(clelia, clel(j), iocc=i, nbval=0, nbret=ioc)
            ioc = -ioc
            if (ioc .ne. 1) then
                vali(1) = i
                vali(2) = ioc
                vali(3) = 1
                valk(1) = clel(j)
                call utmess('F', 'ALGORITH11_98', sk=valk(1), ni=3, vali=vali)
            else
                call getvtx(clelia, clel(j), iocc=i, scal=nomcou, nbret=ibid)
!
!  VERIFICATION EXISTANCE DE LA SOUS-STRUCTURE
!
                call jenonu(jexnom(repsst, nomcou), iret)
                if (iret .eq. 0) then
                    vali(1) = i
                    valk(1) = nomcou
                    call utmess('F', 'ALGORITH11_99', sk=valk(1), si=vali(1))
                end if
                zk8(ldlid+(j-1)*2) = nomcou
            end if
        end do
!
!  BOUCLE SUR LES INTERFACES
!
        do j = 3, 4
            call getvtx(clelia, clel(j), iocc=i, nbval=0, nbret=ioc)
            ioc = -ioc
            if (ioc .ne. 1) then
                vali(1) = i
                vali(2) = ioc
                vali(3) = 1
                valk(1) = clel(j)
                call utmess('F', 'ALGORITH11_98', sk=valk(1), ni=3, vali=vali)
            else
                call getvtx(clelia, clel(j), iocc=i, scal=nomcou, nbret=ibid)
            end if
!
!  VERIFICATION DE L'EXISTANCE DE L'INTERFACE
!
            nomsst = zk8(ldlid+(j-3)*2)
            call mgutdm(nomres, nomsst, ibid, 'NOM_LIST_INTERF', ibid, &
                        lintf)
            if (lintf(1:2) .eq. ' ') then
                call utmess('F', 'ALGORITH12_3', sk=nomsst)
            end if
            call jenonu(jexnom(lintf//'.IDC_NOMS', nomcou), iret)
            if (iret .eq. 0) then
                vali(1) = i
                valk(1) = nomsst
                valk(2) = '   '
                valk(3) = nomcou
                call utmess('F', 'ALGORITH12_2', nk=3, valk=valk, si=vali(1))
            end if
            zk8(ldlid+(j-3)*2+1) = nomcou
        end do
!  ON INITIALISE L'ORDONANCEMENT A NON
        zk8(ldlid+4) = 'NON'
    end do
!
!
!
    call jedema()
end subroutine
