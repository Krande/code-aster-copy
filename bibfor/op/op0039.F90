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
subroutine op0039()
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/irmfac.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/ulaffe.h"
#include "asterfort/ulexis.h"
#include "asterfort/ultype.h"
#include "asterfort/ulopen.h"
#include "asterfort/utmess.h"
#include "asterfort/asmpi_info.h"
!
! --------------------------------------------------------------------------------------------------
!
! IMPR_RESU
!
! --------------------------------------------------------------------------------------------------
!
    mpi_int :: mrank, msize
    character(len=16), parameter :: keywf = 'RESU'
    integer(kind=8) :: keywfNb, keywfIocc
    integer(kind=8) :: fileUnit, fileVersion
    character(len=1) :: fileType
    character(len=3) :: fichierUnique
    character(len=8) :: fileFormat
    character(len=16) :: fileName, resuType
    integer(kind=8) :: nbMesh, nbResu, nbField, nbRet, nbNodeCmp
    integer(kind=8) :: nn, rank, nbprocs
    real(kind=8) :: fileVersionR
    real(kind=8), parameter :: eps = 1.0d-6
    character(len=8) :: model, mesh, resultMesh, proc0, ispar, nameType
    character(len=19) :: result
    aster_logical :: lResu, lMesh, lfichUniq, lNomCas
!
! ------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
!
! - Get rank of processor
!
    call asmpi_info(rank=mrank, size=msize)
    rank = to_aster_int(mrank)
    nbprocs = to_aster_int(msize)
!
    call getvtx(' ', 'PROC0', scal=proc0, nbret=nbRet)
    if (nbRet .ne. 1) then
        proc0 = 'OUI'
    end if
    lfichUniq = .false._1
    lNomCas = ASTER_FALSE
    result = ' '
!
!XX if (nbIndexes .eq. 0) then
    call getfac(keywf, keywfNb)
! ----- Check for parallel mesh vs restriction on proc0
    call getvid(keywf, 'MAILLAGE', iocc=1, scal=mesh, nbret=nbMesh)
    if (nbMesh .eq. 0) then
        call getvid(keywf, 'RESULTAT', iocc=1, scal=result, nbret=nbResu)
        if (nbResu .ne. 0) then
            call dismoi('NOM_MAILLA', result, 'RESULTAT', repk=mesh)
        else
            call getvid(keywf, 'CHAM_GD', iocc=1, scal=result, nbret=nbField)
            ASSERT(nbField .eq. 1)
            call dismoi('NOM_MAILLA', result, 'CHAMP', repk=mesh)
        end if
    end if
    call dismoi('PARALLEL_MESH', mesh, 'MAILLAGE', repk=ispar)
!       if parallel mesh, print on all processors
    if (ispar .eq. 'NON' .and. proc0 .eq. 'OUI' .and. rank .ne. 0) then
!           if only on 'proc0' enabled and rank != 0
        call utmess('I', 'APPELMPI_3', si=rank)
        goto 999
    end if
! ----- Check keyword NOEUD_CMP
    do keywfIocc = 1, keywfNb
        call getvtx(keywf, 'NOEUD_CMP', iocc=keywfIocc, nbval=0, nbret=nbNodeCmp)
        if (nbNodeCmp .ne. 0) then
            nn = nbNodeCmp/2
            if (2*nn .ne. nbNodeCmp) then
                call utmess('F', 'RESULT3_65')
            end if
        end if
    end do
! ----- Get model
    model = ' '
    call getvid(' ', 'MODELE', scal=model, nbret=nbRet)
! ----- Get format of file
    call getvtx(' ', 'FORMAT', scal=fileFormat, nbret=nbRet)
! ----- Check consistency of mesh for IDEAS
    if (fileFormat .eq. 'IDEAS') then
        call getvid(keywf, 'RESULTAT', iocc=1, scal=result, nbret=nbResu)
        call getvid(keywf, 'MAILLAGE', iocc=1, scal=mesh, nbret=nbMesh)
        if (nbResu*nbMesh .gt. 0) then
            call dismoi('NOM_MAILLA', result, 'RESULTAT', repk=resultMesh)
            if (resultMesh .ne. mesh) then
                call utmess('A', 'RESULT3_74')
            end if
        end if
    end if
! ----- Get version of file
    fileVersion = 0
    if (fileFormat .eq. 'IDEAS') then
        fileVersion = 5
        call getvis(' ', 'VERSION', scal=fileVersion, nbret=nbRet)
    else if (fileFormat .eq. 'GMSH') then
        fileVersion = 1
        fileVersionR = 1.d0
        call getvr8(' ', 'VERSION', scal=fileVersionR, nbret=nbRet)
        if (fileVersionR .gt. 1.0d0-eps .and. fileVersionR .lt. 1.0d0+eps) then
            fileVersion = 1
        else if (fileVersionR .gt. 1.2d0-eps .and. fileVersionR .lt. 1.2d0+eps) then
            fileVersion = 2
        end if
    end if
! ----- Get logical unit of file
    fileUnit = 0
    call getvis(' ', 'UNITE', scal=fileUnit, nbret=nbRet)
! ----- Open file
    fileName = 'F_'//fileFormat
    if (.not. ulexis(fileUnit)) then
        if (fileFormat .eq. 'MED') then
            call ulaffe(fileUnit, ' ', fileName, 'NEW', 'O')
        else
            call ulopen(fileUnit, ' ', fileName, 'NEW', 'O')
        end if
    elseif (fileFormat .eq. 'MED') then
        call ultype(fileUnit, fileType)
        if (fileType .ne. 'B' .and. fileType .ne. 'L') then
            call utmess('A', 'RESULT3_12')
        end if
    end if
    if (fileFormat .eq. 'MED') then
        if (nbprocs .gt. 1) then
            call getvtx(' ', 'FICHIER_UNIQUE', scal=fichierUnique, nbret=nbRet)
            if (fichierUnique .eq. 'OUI' .and. ispar == "OUI") then
                lfichUniq = .true._1
                if (proc0 .eq. 'OUI') call utmess('F', 'RESULT3_36')
            end if
        end if
        call ultype(fileUnit, fileType)
        if (fileType .ne. 'B' .and. fileType .ne. 'L') then
            call utmess('A', 'RESULT3_12')
        end if
    end if
! ----- Check consistency for GMSH
    if (fileFormat .eq. 'GMSH') then
        lMesh = ASTER_FALSE
        lResu = ASTER_FALSE
        do keywfIocc = 1, keywfNb
            call getvid(keywf, 'MAILLAGE', iocc=keywfIocc, scal=mesh, nbret=nbMesh)
            call getvid(keywf, 'RESULTAT', iocc=keywfIocc, scal=result, nbret=nbResu)
            call getvid(keywf, 'CHAM_GD', iocc=keywfIocc, scal=result, nbret=nbField)
            if (nbResu .ne. 0 .or. nbField .ne. 0) then
                lResu = ASTER_TRUE
                goto 220
            end if
            if (nbMesh .ne. 0) then
                lMesh = ASTER_TRUE
            end if
220         continue
        end do
        if (lMesh .and. lResu) then
            call utmess('F', 'RESULT3_68')
        end if
    end if
! ----- Check consistency for ASTER
    if (fileFormat .eq. 'ASTER') then
        lResu = ASTER_FALSE
        do keywfIocc = 1, keywfNb
            call getvid(keywf, 'RESULTAT', iocc=keywfIocc, scal=result, nbret=nbResu)
            call getvid(keywf, 'CHAM_GD', iocc=keywfIocc, scal=result, nbret=nbField)
            if (nbResu .ne. 0 .or. nbField .ne. 0) then
                lResu = ASTER_TRUE
                exit
            end if
        end do
        if (lResu) then
            call utmess('A', 'RESULT3_67')
        end if
    end if
! ----- Loop on factor keywords
    do keywfIocc = 1, keywfNb
        if (fileFormat .eq. 'MED') then
            call gettco(result, resuType)
            if (resuType .eq. "MULT_ELAS") then
                call getvtx(keywf, 'TYPE_NOM', iocc=keywfIocc, nbval=0, nbret=nbRet)
                if (nbRet .ne. 0) then
                    call getvtx(keywf, 'TYPE_NOM', iocc=keywfIocc, scal=nameType)
                    if (nameType .eq. "NOM_CAS") then
                        lNomCas = ASTER_TRUE
                    end if
                end if
            end if
        end if
        call irmfac(keywfIocc, fileFormat, fileUnit, fileVersion, model, lfichUniq, lNomCas)
    end do
    if (fileFormat .ne. 'MED') then
        flush (fileUnit)
        if (fileUnit .ne. 6) then
            call ulopen(-fileUnit, ' ', ' ', ' ', ' ')
        end if
    end if
!XX endif
999 continue
    call jedema()
end subroutine
