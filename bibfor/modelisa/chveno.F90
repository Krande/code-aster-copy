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
subroutine chveno(valeType, meshZ, modelZ)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8prem.h"
#include "asterfort/chbord.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/orilma.h"
#include "asterfort/ornorm.h"
#include "asterfort/utmamo.h"
#include "asterfort/utmess.h"
#include "asterfort/utmotp.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    character(len=4), intent(in) :: valeType
    character(len=*), intent(in) :: meshZ, modelZ
!
!      OPERATEURS :     AFFE_CHAR_MECA ET AFFE_CHAR_MECA_C
!                                      ET AFFE_CHAR_MECA_F
!                       DEFI_CONTACT
!
!     VERIFICATION DES NORMALES AUX MAILLES SURFACIQUES EN 3D
!     ET LINEIQUES EN 2D
!     V1 : ON VERIFIE QUE LES NORMALES SONT HOMOGENES
!     V2 : ON VERIFIE QUE LES NORMALES SONT SORTANTES
!
!-----------------------------------------------------------------------
    integer(kind=8) :: nbt
    parameter(nbt=5)
    integer(kind=8) :: ier, iret, zero
    integer(kind=8) :: imfac, nbmfac, n, geomDime, ndim1, vali
    integer(kind=8) :: iocc, nocc, ic, nbmc, iobj, nbobj, iCell, impb, nbCell
    integer(kind=8) :: cellNume, idtyma, cellTypeNume, nbmapr, nbmabo, ntrait
    integer(kind=8) :: jcoor
    integer(kind=8) :: if1, if2, if3, imf1, imf2, ipres, idnor, idtan
    integer(kind=8) :: norien, norie1, norie2, jlima, nbmamo, nconex
    real(kind=8) :: dnor
    aster_logical :: reorie, mcfl(nbt)
    character(len=8) :: mot, mesh, model, cellTypeName, algo
    character(len=16) :: mcft(nbt), motfac, valmc(4), typmc(4)
    character(len=19) :: limamo
    character(len=24) :: grmama, nogr, cellName
    character(len=24) :: valk(2)
    character(len=24), pointer :: objet(:) => null()
    integer(kind=8), pointer :: listCellNume(:) => null(), typmail(:) => null()
    integer(kind=8), pointer :: listCellPrin(:) => null(), listCellBord(:) => null()
    integer(kind=8), pointer :: listCellBord2(:) => null()
!
    data mcft/'FACE_IMPO', 'PRES_REP', 'FORCE_COQUE',&
     &            'EFFE_FOND', 'ZONE'/
!
!     LA NORMALE DOIT ETRE SORTANTE:
    data mcfl/.true., .true., .false.,&
     &            .true., .true./
!     ------------------------------------------------------------------
!
!     INITIALISATIONS
    ier = 0
    zero = 0
    reorie = .false.
    algo = ''
!
!     NOMBRE DE MOTS-CLES FACTEUR A VERIFIER
    nbmfac = nbt
!
    mesh = meshZ
    model = modelZ
    grmama = mesh//'.GROUPEMA'
    limamo = '&&CHVENO.MAIL_MODEL'
!
    call getvtx(' ', 'VERI_NORM', scal=mot, nbret=n)
    if (mot .eq. 'NON') nbmfac = 0
!
    geomDime = 0
    call dismoi('DIM_GEOM', modelZ, 'MODELE', repi=geomDime)
!
    call jeveuo(mesh//'.COORDO    .VALE', 'L', jcoor)
!
    do imfac = 1, nbmfac
        motfac = mcft(imfac)
        call getfac(motfac, nocc)
        do iocc = 1, nocc
!         POUR CERTAINS MOTS-CLES, IL NE FAUT TESTER QUE
!         POUR CERTAINS CHARGEMENTS
            if (motfac .eq. 'FACE_IMPO') then
                ipres = utmotp(valeType, motfac, iocc, 'PRES')
                idnor = utmotp(valeType, motfac, iocc, 'DNOR')
                idtan = utmotp(valeType, motfac, iocc, 'DTAN')
                if (ipres .eq. 0 .and. idnor .eq. 0 .and. idtan .eq. 0) goto 200
                if (idnor .ne. 0) then
                    if (valeType .eq. 'REEL') then
                        call getvr8(motfac, 'DNOR', iocc=iocc, scal=dnor, nbret=n)
                        if (abs(dnor) .le. r8prem()) goto 200
                    end if
                end if
            else if (motfac .eq. 'FORCE_COQUE') then
                ipres = utmotp(valeType, motfac, iocc, 'PRES')
                if1 = utmotp(valeType, motfac, iocc, 'F1  ')
                if2 = utmotp(valeType, motfac, iocc, 'F2  ')
                if3 = utmotp(valeType, motfac, iocc, 'F3  ')
                imf1 = utmotp(valeType, motfac, iocc, 'MF1 ')
                imf2 = utmotp(valeType, motfac, iocc, 'MF2 ')
                if (ipres .eq. 0 .and. if1 .eq. 0 .and. if2 .eq. 0 .and. if3 .eq. 0 .and. imf1 &
                    .eq. 0 .and. imf2 .eq. 0) goto 200
            end if
!
            if (motfac .eq. 'ZONE') then
                nbmc = 4
                valmc(1) = 'GROUP_MA_ESCL'
                valmc(2) = 'GROUP_MA_MAIT'
                valmc(3) = 'MAILLE_ESCL'
                valmc(4) = 'MAILLE_MAIT'
                typmc(1) = 'GROUP_MA'
                typmc(2) = 'GROUP_MA'
                typmc(3) = 'MAILLE'
                typmc(4) = 'MAILLE'
            else
                nbmc = 2
                valmc(1) = 'GROUP_MA'
                valmc(2) = 'MAILLE'
                typmc(1) = 'GROUP_MA'
                typmc(2) = 'MAILLE'
            end if
!
! ---     RECUPERATION DE LA DIMENSION DU PROBLEME
!
            do ic = 1, nbmc
                call getvtx(motfac, valmc(ic), iocc=iocc, nbval=0, nbret=nbobj)
                if (nbobj .eq. 0) goto 210
!
                nbobj = -nbobj
                AS_ALLOCATE(vk24=objet, size=nbobj)
                call getvem(meshZ, typmc(ic), motfac, valmc(ic), iocc, &
                            nbobj, objet, nbobj)
                if (typmc(ic) .eq. 'GROUP_MA') then
                    do iobj = 1, nbobj
                        nogr = objet(iobj)
                        if (motfac .eq. 'ZONE') then

                            call getvtx(motfac, 'ALGO_CONT', iocc=1, scal=algo, nbret=iret)
!
! ---             RECUPERATION DU NOMBRE DE MAILLES DU GROUP_MA :
!                 ---------------------------------------------
                            call jelira(jexnom(grmama, nogr), 'LONUTI', nbCell)
                            call jeveuo(jexnom(grmama, nogr), 'L', vi=listCellNume)
!
                            do iCell = 1, nbCell
                                cellNume = listCellNume(iCell)
                                cellName = int_to_char8(cellNume)
!
! ---               NUMERO DE LA MAILLE
!                   ------------------
                                call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
                                cellTypeNume = typmail(cellNume)
!
! ---               TYPE DE LA MAILLE :
!                   -----------------
                                call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
!
! ---               CAS D'UNE MAILLE POINT
!                   ----------------------
                                if (cellTypeName(1:3) .eq. 'POI') then
!                     ON SAUTE
                                    goto 211
!
! ---               CAS D'UNE MAILLE SEG
!                   --------------------
                                else if (cellTypeName(1:3) .eq. 'SEG') then
                                    ndim1 = 2
                                    if (geomDime .ne. ndim1) then
!                       ON SAUTE
                                        goto 211
                                    end if
!
                                end if
                            end do
!
! ---           FIN DE BOUCLE SUR LES MAILLES DU GROUP_MA
!
                        end if
                        norie1 = 0
                        norie2 = 0
                        call jelira(jexnom(grmama, nogr), 'LONUTI', nbCell)
                        call jeveuo(jexnom(grmama, nogr), 'L', vi=listCellNume)
!
                        if (mcfl(ic) .and. (nbCell .gt. 0)) then
!
                            call wkvect('&&CHVENO.MAILLE_BORD', 'V V I', nbCell, vi=listCellBord2)
                            call chbord(modelZ, nbCell, listCellNume, listCellBord2, nbmapr, &
                                        nbmabo)
                            if (nbmapr .eq. nbCell .and. nbmabo .eq. 0) then
                                call ornorm(mesh, listCellNume, nbCell, reorie, norie1, nconex)
                                if (nconex .gt. 1) then
                                    call utmess('F', 'MESH3_99')
                                end if
                            elseif ((nbmapr .eq. 0 .and. &
                                     nbmabo .eq. nbCell) .or. (motfac .eq. &
                                                               'ZONE')) then
                                if (motfac .eq. 'ZONE') then
                                    nbmamo = 0
                                    jlima = 1
                                else
                                    call utmamo(model, nbmamo, limamo)
                                    call jeveuo(limamo, 'L', jlima)
                                end if
                                call orilma(mesh, geomDime, listCellNume, nbCell, norie1, &
                                            ntrait, reorie, nbmamo, zi(jlima))
                                if ((algo .eq. 'LAC') .and. (ntrait .ne. 0)) then
                                    call utmess('A', 'CONTACT2_20')
                                end if
                                call jedetr(limamo)
                            elseif (nbmapr .eq. 0 .and. nbmabo .eq. 0) &
                                then
                                call ornorm(mesh, listCellNume, nbCell, reorie, norie1, nconex)
                                if (nconex .gt. 1) then
                                    call utmess('F', 'MESH3_99')
                                end if
                            else
                                call wkvect('&&CHVENO.PRIN', 'V V I', nbmapr, vi=listCellPrin)
                                call wkvect('&&CHVENO.BORD', 'V V I', nbmabo, vi=listCellBord)
                                nbmapr = 0
                                nbmabo = 0
                                do impb = 1, nbCell
                                    if (listCellBord2(impb) .eq. 0) then
                                        nbmapr = nbmapr+1
                                        listCellPrin(nbmapr) = listCellNume(impb)
                                    else
                                        nbmabo = nbmabo+1
                                        listCellBord(nbmabo) = listCellNume(impb)
                                    end if
                                end do
                                call ornorm(mesh, listCellPrin, nbmapr, reorie, norie1, nconex)
                                if (nconex .gt. 1) then
                                    call utmess('F', 'MESH3_99')
                                end if
                                call orilma(mesh, geomDime, listCellBord, nbmabo, norie1, &
                                            ntrait, reorie, 0, [0])
                                call jedetr('&&CHVENO.PRIN')
                                call jedetr('&&CHVENO.BORD')
                            end if
                            call jedetr('&&CHVENO.MAILLE_BORD')
                        else
                            call ornorm(mesh, listCellNume, nbCell, reorie, norie2, nconex)
                            if (nconex .gt. 1) then
                                call utmess('F', 'MESH3_99')
                            end if
                        end if
                        norien = norie1+norie2
                        if (norien .ne. 0) then
                            ier = ier+1
                            valk(1) = nogr
                            vali = norien
                            call utmess('E', 'MODELISA8_56', sk=valk(1), si=vali)
                        end if
                    end do
!
! ----------CAS DES MAILLES :
!           ---------------
                else
                    AS_ALLOCATE(vi=listCellNume, size=nbobj)
                    do iobj = 1, nbobj
                        cellName = objet(iobj)
                        cellNume = char8_to_int(cellName)
                        listCellNume(iobj) = cellNume
                        if (motfac .eq. 'ZONE') then
                            call jeveuo(mesh//'.TYPMAIL', 'L', idtyma)
                            cellTypeNume = zi(idtyma+cellNume-1)
                            call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
!
! ---             CAS D'UNE MAILLE POINT
!                 -----------------------
                            if (cellTypeName(1:3) .eq. 'POI') then
!                   ON SAUTE
                                goto 211
!
! ---             CAS D'UNE MAILLE SEG
!                 --------------------
                            else if (cellTypeName(1:3) .eq. 'SEG') then
                                ndim1 = 2
                                if (geomDime .ne. ndim1) then
!                     ON SAUTE
                                    goto 211
                                end if
                            end if
!
                        end if
                    end do
                    norie1 = 0
                    norie2 = 0
                    if (mcfl(ic)) then
                        call wkvect('&&CHVENO.MAILLE_BORD', 'V V I', nbobj, vi=listCellBord2)
                        call chbord(modelZ, nbobj, listCellNume, listCellBord2, nbmapr, &
                                    nbmabo)
                        if (nbmapr .eq. nbobj .and. nbmabo .eq. 0) then
                            call ornorm(mesh, listCellNume, nbobj, reorie, norie1, nconex)
                            if (nconex .gt. 1) then
                                call utmess('F', 'MESH3_99')
                            end if
                        elseif ((nbmapr .eq. 0 .and. nbmabo .eq. nbobj) &
                                .or. (motfac .eq. 'ZONE')) then
                            if (motfac .eq. 'ZONE') then
                                nbmamo = 0
                                jlima = 1
                            else
                                call utmamo(model, nbmamo, limamo)
                                call jeveuo(limamo, 'L', jlima)
                            end if
                            call orilma(mesh, geomDime, listCellNume, nbobj, norie1, &
                                        ntrait, reorie, nbmamo, zi(jlima))
                            call jedetr(limamo)
                        else if (nbmapr .eq. 0 .and. nbmabo .eq. 0) then
                            call ornorm(mesh, listCellNume, nbobj, reorie, norie1, nconex)
                            if (nconex .gt. 1) then
                                call utmess('F', 'MESH3_99')
                            end if
                        else
                            call wkvect('&&CHVENO.PRIN', 'V V I', nbmapr, vi=listCellPrin)
                            call wkvect('&&CHVENO.BORD', 'V V I', nbmabo, vi=listCellBord)
                            nbmapr = 0
                            nbmabo = 0
                            do impb = 1, nbobj
                                if (listCellBord2(impb) .eq. 0) then
                                    nbmapr = nbmapr+1
                                    listCellPrin(nbmapr) = listCellNume(impb)
                                else
                                    nbmabo = nbmabo+1
                                    listCellBord(nbmabo) = listCellNume(impb)
                                end if
                            end do
                            call ornorm(mesh, listCellPrin, nbmapr, reorie, norie1, nconex)
                            if (nconex .gt. 1) then
                                call utmess('F', 'MESH3_99')
                            end if
                            call orilma(mesh, geomDime, listCellBord, nbmabo, norie1, &
                                        ntrait, reorie, 0, [0])
                            call jedetr('&&CHVENO.PRIN')
                            call jedetr('&&CHVENO.BORD')
                        end if
                        call jedetr('&&CHVENO.MAILLE_BORD')
                    else
                        call ornorm(mesh, listCellNume, nbobj, reorie, norie2, nconex)
                        if (nconex .gt. 1) then
                            call utmess('F', 'MESH3_99')
                        end if
                    end if
                    norien = norie1+norie2
                    if (norien .ne. 0) then
                        ier = ier+1
                        valk(1) = cellName
                        call utmess('E', 'MODELISA8_57', sk=valk(1))
                    end if
                    AS_DEALLOCATE(vi=listCellNume)
                end if
211             continue

                AS_DEALLOCATE(vk24=objet)
210             continue
            end do
200         continue
        end do
    end do
!
    if (ier .ne. 0) then
        call utmess('F', 'MODELISA4_24')
    end if
!
end subroutine
