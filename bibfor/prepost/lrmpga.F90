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

subroutine lrmpga(fileUnit, ligrel, MEDFieldName, nbCell, pgmail, &
                  pgmmil, spmmil, ntypel, npgmax, indpg, &
                  numpt, numord, option, param, nomaas)
!
! person_in_charge: nicolas.sellenet at edf.fr
!     LECTURE FICHIER MED - LOCALISATION POINTS DE GAUSS
!     -    -            -                  -         --
!-----------------------------------------------------------------------
!     IN :
!       NROFIC : UNITE LOGIQUE DU FICHIER MED
!       LIGREL : NOM DU LIGREL
!       MEDFieldName : NOM DU CHAMP MED
!       NBMA   : NOMBRE DE MAILLES DU MAILLAGE
!       NTYPEL : NOMBRE TOTAL DE TYPES DE MAILLE (=27)
!       NPGMAX : NOMBRE DE PG MAX (=27)
!       NUMPT  : NUMERO DE PAS DE TEMPS EVENTUEL
!       NUMORD : NUMERO D'ORDRE EVENTUEL DU CHAMP
!
!     OUT:
!       PGMAIL : NOMBRE DE POINTS DE GAUSS PAR MAILLE (ASTER)
!       PGMMIL : NOMBRE DE POINTS DE GAUSS PAR MAILLE (MED)
!     IN/OUT:
!       INDPG  : TABLEAU D'INDICES DETERMINANT L'ORDRE DES POINTS
!                DE GAUSS DANS UN ELEMENT DE REFERENCE :
!                INDPG(K,I_MED)=I_ASTER :
!                   - K = NUM DU TYPE DE MAILLE
!                   - I_ASTER = NUMERO LOCAL DU PG DANS L'ELEMENT
!                               DE REFERENCE ASTER
!                   - I_MED  = NUMERO LOCAL DU PG DANS L'ELEMENT
!                               DE REFERENCE MED
!
!-----------------------------------------------------------------------
!
    use as_med_module, only: as_med_open
    implicit none
!
#include "jeveux.h"
#include "asterfort/as_mfdfin.h"
#include "asterfort/as_mfdncn.h"
#include "asterfort/as_mfdonp.h"
#include "asterfort/as_mfdonv.h"
#include "asterfort/as_mlcnlc.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lrcmpr.h"
#include "asterfort/lrvcpg.h"
#include "asterfort/modat2.h"
#include "asterfort/typele.h"
#include "asterfort/ulisog.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: fileUnit, nbCell, ntypel, npgmax, numpt, numord
    integer(kind=8) :: pgmail(nbCell), pgmmil(nbCell), spmmil(nbCell), indpg(ntypel, npgmax)
    character(len=8) :: param
    character(len=19) :: ligrel
    character(len=24) :: option
    character(len=*) :: MEDFieldName
    character(len=8) :: nomaas
!
    character(len=6) :: nompro
    parameter(nompro='LRMPGA')
!
    integer(kind=8), parameter :: nbCellType = 19
    integer(kind=8), parameter :: MED_ACC_RDONLY = 0
    integer(kind=8), parameter :: MED_CELL = 0
    integer(kind=8), parameter :: MED_COMPACT_STMODE = 2
!
    integer(kind=8) :: ifm, nivinf, nbCmp, nbProfile, jnopro
    med_idt :: MEDFileIden
    integer(kind=8) :: codret, nbLocalizations, iret, igrel
    integer(kind=8) :: iTypeCellInField, nbgrel, jnonpg, iProfile, lgproa, codre2
    integer(kind=8) :: typeElemNume, i, iCellType, nbStep
    integer(kind=8) :: AsterNbPg, MEDNbPg, nbElem, lonmax
    integer(kind=8) :: nufgpg, jnoloc
    integer(kind=8) :: modelDime, jvGrel
    integer(kind=8) :: nbTypeCellInField, iElem
    integer(kind=8) :: jngaok
    integer(kind=8) :: ipg, ipgm, jperm
    integer(kind=8) :: npr, nbValues, l_fapg
    integer(kind=8) :: iopt, imod, jvModeLoc, igrd, jvDescrigd, nec, nbsp
    integer(kind=8) :: typeCellNume, jadproa, jtypma, jmedtoaster, ima, imed, itypma
    character(len=8) :: typma
    integer(kind=8), parameter :: MEDIterMesh = 1
!
    character(len=1) :: fileState
    character(len=8) :: cellTypeIden, fapg, elref, typeCellName
    character(len=16) :: nofgpg
    character(len=24) :: liel, nolipr, nlnbpg
    character(len=64) :: profileName, localizationName, MEDMeshName
    character(len=64) :: k24Dummy1, k24Dummy2
    character(len=200) :: MEDFileName
    character(len=255) :: fileName
    character(len=16), pointer :: cname(:) => null()
    character(len=16), pointer :: cunit(:) => null()
    character(len=8), pointer :: typema(:) => null()
    character(len=8), pointer :: asterCellType(:) => null()
    integer(kind=8), pointer :: MEDCellType(:) => null()
    !character(len=8), pointer :: asterElemType(:) => null()
    integer(kind=8), pointer :: tmfpg(:) => null()
    aster_logical :: l_parallel_mesh
    integer(kind=8), pointer :: repe(:) => null()
!
    integer(kind=8), parameter :: MED_GEOMETRY_TYPE(nbCellType) = &
                                  (/1, 102, 103, 104, &
                                    203, 204, 206, 207, &
                                    208, 209, 304, 305, &
                                    306, 308, 310, 313, &
                                    315, 320, 327/)
    character(len=8), parameter :: AsterAllCellType(nbCellType) = &
                                   (/'POI1    ', 'SEG2    ', 'SEG3    ', 'SEG4    ', &
                                     'TRIA3   ', 'QUAD4   ', 'TRIA6   ', 'TRIA7   ', &
                                     'QUAD8   ', 'QUAD9   ', 'TETRA4  ', 'PYRAM5  ', &
                                     'PENTA6  ', 'HEXA8   ', 'TETRA10 ', 'PYRAM13 ', &
                                     'PENTA15 ', 'HEXA20  ', 'HEXA27  '/)
!
!-----------------------------------------------------------------------
!
    call jemarq()
!
    call infniv(ifm, nivinf)
!
    if (nivinf .gt. 1) then
        write (ifm, 101) 'DEBUT DE '//nompro
    end if
    l_parallel_mesh = isParallelMesh(nomaas)

    call dismoi('DIM_GEOM', ligrel(1:8), 'MODELE', repi=modelDime)
    if (.not. (modelDime .eq. 2 .or. modelDime .eq. 3)) then
        call utmess('F', 'MODELISA2_6')
    end if
!
!  ======================================
!  == 1 : EXPLOITATION DU FICHIER MED  ==
!  ======================================
!
!  == 1.1. INITIALISATIONS
    do i = 1, nbCell
        pgmail(i) = 0
    end do

! - Get name of file from logical unit
    call ulisog(fileUnit, fileName, fileState)

! - Get name of MED file
    if (fileName(1:1) .eq. ' ') then
        call codent(fileUnit, 'G', cellTypeIden)
        MEDFileName = 'fort.'//cellTypeIden
    else
        MEDFileName = fileName(1:200)
    end if

! - Open MED file
    call as_med_open(MEDFileIden, MEDFileName, MED_ACC_RDONLY, codret)

! - Nombre de localisations dans le fichier MED
    nbLocalizations = 0
    call as_mlcnlc(MEDFileIden, nbLocalizations, iret)
!
!  == 1.3. A PARTIR DU NOM DU CHAMP MED ET DE L'INDICE DU PAS DE TEMPS,
!      ON RECUPERE POUR CHAQUE TYPE DE MAILLE PRESENT:
!      - LE NOM DU TYPE GEOMETRIQUE   : ZK8(JTYMED)
!      - UN MARQUEUR DE CORRESPONDANCE : ZI(JNGAOK)
!      REMARQUE: L'INDICE DU PAS DE TEMPS EST OBTENU EN PARCOURANT
!      LA LISTE DES PAS DE TEMPS (BOUCLE 39)
    AS_ALLOCATE(vk8=asterCellType, size=nbCellType)
    AS_ALLOCATE(vi=MEDCellType, size=nbCellType)
    call wkvect('&&LRMPGA_TYPGEO_OKPG_MED', 'V V I', nbCellType, jngaok)
    call wkvect('&&LRMPGA_TYPGEO_NOMLOC', 'V V K80', 2*nbCellType, jnoloc)
    nbTypeCellInField = 0
    if (nivinf .gt. 1) then
        write (ifm, 201) MEDFieldName
    end if

! - Get number of components in field
    call as_mfdncn(MEDFileIden, MEDFieldName, nbCmp, iret)

! - Get parameters about field
    AS_ALLOCATE(vk16=cname, size=nbCmp)
    AS_ALLOCATE(vk16=cunit, size=nbCmp)
    call as_mfdfin(MEDFileIden, MEDFieldName, MEDMeshName, nbStep, cname(1), &
                   cunit(1), iret)
    AS_DEALLOCATE(vk16=cname)
    AS_DEALLOCATE(vk16=cunit)

! - Get profiles
    if (nbStep .gt. 0) then
        do iCellType = 1, nbCellType

! --------- Current cell
            call codent(iCellType, 'G', cellTypeIden)

! --------- Get number of profiles
            call as_mfdonp(MEDFileIden, &
                           MEDFieldName, numpt, numord, &
                           MED_CELL, MED_GEOMETRY_TYPE(iCellType), &
                           MEDIterMesh, MEDMeshName, &
                           k24Dummy1, k24Dummy2, &
                           nbProfile, codret)

! --------- Get profiles
            if (nbProfile .ne. 0) then
                nolipr = '&&LRMPGA'//cellTypeIden//'PR'
                nlnbpg = '&&LRMPGA'//cellTypeIden//'PG'
                call wkvect(nolipr, 'V V K80', 2*nbProfile, jnopro)
                call wkvect(nlnbpg, 'V V I', nbProfile, jnonpg)

                do iProfile = 1, nbProfile
                    call as_mfdonv(MEDFileIden, MEDFieldName, &
                                   MED_CELL, MED_GEOMETRY_TYPE(iCellType), MEDMeshName, &
                                   numpt, numord, iProfile, profileName, MED_COMPACT_STMODE, &
                                   npr, localizationName, MEDNbPg, nbValues, iret)
                    ASSERT(iret .eq. 0)
                    zk80(jnopro+2*iProfile-2) = profileName
                    zk80(jnopro+2*iProfile-1) = localizationName
                    zi(jnonpg+iProfile-1) = MEDNbPg
                end do
            else
                call as_mfdonv(MEDFileIden, MEDFieldName, &
                               MED_CELL, MED_GEOMETRY_TYPE(iCellType), MEDMeshName, &
                               numpt, numord, 1, profileName, MED_COMPACT_STMODE, &
                               npr, localizationName, MEDNbPg, nbValues, iret)
                ASSERT(iret .eq. 0)
                nolipr = '&&LRMPGA'//cellTypeIden//'PR'
                nlnbpg = '&&LRMPGA'//cellTypeIden//'PG'
                call wkvect(nolipr, 'V V K80', 2, jnopro)
                call wkvect(nlnbpg, 'V V I', 1, jnonpg)
                zk80(jnopro) = ' '
                zk80(jnopro+1) = localizationName
                zi(jnonpg) = -1
            end if

            if (nbValues .gt. 0) then
                nbTypeCellInField = nbTypeCellInField+1
                asterCellType(nbTypeCellInField) = AsterAllCellType(iCellType)
                !asterElemType(nbTypeCellInField) = AsterAllElemType(iCellType)
                zk80(jnoloc+2*nbTypeCellInField-2) = nolipr
                zk80(jnoloc+2*nbTypeCellInField-1) = nlnbpg
                zi(jngaok+nbTypeCellInField-1) = 0
                MEDCellType(nbTypeCellInField) = MED_GEOMETRY_TYPE(iCellType)
            else
                call jedetr(nolipr)
                call jedetr(nlnbpg)
                MEDNbPg = 0
            end if
!
        end do
    end if
    if (nivinf .gt. 1) then
        write (ifm, *) ' '
    end if
!
    if (nbTypeCellInField .eq. 0) then
        call utmess('F', 'MED_77', sk=MEDFieldName, si=fileUnit)
    end if
!

    call jeveuo(ligrel//'.REPE', 'L', vi=repe)
    liel = ligrel//'.LIEL'
    call jelira(liel, 'NMAXOC', nbgrel)
!
!  =========================================
!  == 2 : EXPLOITATION DES DONNEES ASTER  ==
!  =========================================
!
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), iopt)
    call jeveuo('&CATA.TE.TYPEMA', 'L', vk8=typema)

!   ON PARCOURT LES GROUPES D'ELEMENTS PRESENTS DANS LE MODELE
    do igrel = 1, nbgrel

! ----- Type of finite element on current GREL
        typeElemNume = typele(ligrel, igrel)

! ----- Geometric support for this GREL
        typeCellName = typema(typeElemNume)

! ----- Access to GREL
        call jeveuo(jexnum(liel, igrel), 'L', jvGrel)
        call jelira(jexnum(liel, igrel), 'LONMAX', lonmax)
        nbElem = lonmax-1
!
!       MODE LOCAL ASSOCIE AU PARAMETRE PARAM DE L'OPTION OPTION
        imod = modat2(iopt, typeElemNume, param)

        if (imod .eq. 0) then
            AsterNbPg = 0
        else
            call jeveuo(jexnum('&CATA.TE.MODELOC', imod), 'L', jvModeLoc)
!           CHAMP ELGA
            ASSERT(zi(jvModeLoc-1+1) .eq. 3)
!
            igrd = zi(jvModeLoc-1+2)
            call jeveuo(jexnum('&CATA.GD.DESCRIGD', igrd), 'L', jvDescrigd)
            nec = zi(jvDescrigd-1+3)

!           NUMERO ET NOM DE LA FAMILLE GLOBALE DE PTS GAUSS
            nufgpg = zi(jvModeLoc-1+4+nec+1)
            call jenuno(jexnum('&CATA.TM.NOFPG', nufgpg), nofgpg)
            elref = nofgpg(1:8)
            fapg = nofgpg(9:16)

!           NOMBRE DE PG : NBPG
            call jeveuo('&CATA.TM.TMFPG', 'L', vi=tmfpg)
            AsterNbPg = tmfpg(nufgpg)

            l_fapg = len(trim(adjustl(fapg)))

            call jeveuo(nomaas(1:8)//'.TYPMAIL', 'L', jtypma)
            itypma = zi(jtypma+zi(jvGrel)-1)

            call wkvect('&&LRMPGA.MED_TO_ASTER', 'V V I', nbcell, jmedtoaster)
            imed = 0
            do ima = 1, nbcell
                if (zi(jtypma+ima-1) .eq. itypma) then
                    imed = imed+1
                    if (repe(2*(ima-1)+1) .eq. igrel) then
                        zi(jmedtoaster+imed-1) = ima
                    end if
                end if
            end do

! --------- Looking for MED cell
            do iTypeCellInField = 1, nbTypeCellInField
                if (asterCellType(iTypeCellInField) .eq. typeCellName) then
!
!                   VERIFICATION DU NOMBRE DE PG ASTER/MED
!                   COMPARAISON DES COORDONNEES DES PG ASTER/MED
                    call wkvect('&&LRMPGA_PERMUT', 'V V I', AsterNbPg, jperm)

! ----------------- Get number of profiles on this MED cell
                    nolipr = zk80(jnoloc+2*iTypeCellInField-2) (1:24)
                    call jeveuo(nolipr, 'L', jnopro)
                    call jelira(nolipr, 'LONMAX', lonmax)
                    nbProfile = lonmax/2

! ----------------- Acces to number of integration points for these profiles
                    nlnbpg = zk80(jnoloc+2*iTypeCellInField-1) (1:24)
                    call jeveuo(nlnbpg, 'L', jnonpg)

                    do iProfile = 1, nbProfile
! --------------------- Acces to number of integration points for these profiles
                        profileName = zk80(jnopro+2*iProfile-2) (1:64)
                        localizationName = zk80(jnopro+2*iProfile-1) (1:64)
                        MEDNbPg = zi(jnonpg+iProfile-1)

                        ! on ne compare que le meme profile du meme type element
                        if (.not. l_parallel_mesh) then
                            if ((fapg(1:l_fapg) .ne. localizationName(9:8+l_fapg)) .and. &
                                (localizationName .ne. ' ')) then
                                if (localizationName(1:17) .ne. "NOM_LOC_GAUSS_001") goto 999
                            end if
                        end if

! --------------------- Check consistency of integration points between Aster and MED
                        call lrvcpg(MEDFileIden, MEDNbPg, AsterNbPg, &
                                    typeCellName, MEDCellType(iTypeCellInField), &
                                    elref, fapg, nbLocalizations, localizationName, zi(jperm), &
                                    typeCellNume, nbsp, codret)

                        if (codret .ne. 4) then
                            if (profileName .ne. ' ') then
                                call lrcmpr(MEDFileIden, profileName, '&&LRMPGA.TMP', lgproa, &
                                            codre2)
                                call jeveuo('&&LRMPGA.TMP', 'L', jadproa)
                                do i = 1, lgproa
                                    imed = zi(jadproa+i-1)
                                    ima = zi(jmedtoaster+imed-1)
                                    if (ima .gt. 0) then
                                        pgmail(ima) = AsterNbPg
                                        pgmmil(ima) = MEDNbPg
                                        spmmil(ima) = nbsp
                                    end if
                                end do
                                call jedetr('&&LRMPGA.TMP')
                            else
                                do iElem = 1, nbElem
                                    pgmail(zi(jvGrel+iElem-1)) = AsterNbPg
                                    pgmmil(zi(jvGrel+iElem-1)) = MEDNbPg
                                    spmmil(zi(jvGrel+iElem-1)) = nbsp
                                end do
                            end if
                        else
                            codret = 4
                        end if
!
!                       SI LE NBRE PT GAUSS INCORRECT ET PAS DE <F>,
!                       NBPG=0 : RIEN A ECRIRE DANS LRCMVE
                        if (codret .eq. 4) then
                            AsterNbPg = 0
!                           SI PERMUTATIONS AU NIVEAU DES PG ASTER/MED :
                        else if (codret .eq. 1) then
!  ===>                     REMPLISSAGE DU TABLEAU INDPG: CAS OU L'ON A
!                           UNE PERMUTATION DANS LES PG MED/ASTER
                            if (zi(jngaok+iTypeCellInField-1) .eq. 0) then
                                do ipgm = 1, AsterNbPg
                                    indpg(typeCellNume, ipgm) = zi(jperm+ipgm-1)
                                end do
                            else
                                do ipgm = 1, AsterNbPg
                                    ASSERT(indpg(typeCellNume, ipgm) .eq. zi(jperm+ipgm-1))
                                end do
                            end if
                            zi(jngaok+iTypeCellInField-1) = 1
                        else
!  ===>                     SINON REMPLISSAGE DU TABLEAU INDPG: CAS OU L'ON A :
!                            - ABSENCE DE LOCALISATION
!                            - L UN DES PG MED N A PAS ETE IDENTIFIE A UN PG ASTER
!                            - LES PG ASTER/MED CORRESPONDENT
                            if (zi(jngaok+iTypeCellInField-1) .eq. 0) then
                                do ipg = 1, AsterNbPg
                                    indpg(typeCellNume, ipg) = ipg
                                end do
                            else
                                do ipg = 1, AsterNbPg
                                    if (indpg(typeCellNume, ipg) .ne. 0) then
                                        ASSERT(indpg(typeCellNume, ipg) .eq. ipg)
                                    else
                                        indpg(typeCellNume, ipg) = ipg
                                    end if
                                end do
                            end if
                            zi(jngaok+iTypeCellInField-1) = 1
                        end if
999                     continue
                    end do
                    call jedetr('&&LRMPGA_PERMUT')
                end if
            end do
            call jedetr('&&LRMPGA.MED_TO_ASTER')
        end if
    end do
!
!   DESTRUCTION DES TABLEAUX TEMPORAIRES
    do iCellType = 1, 2*nbTypeCellInField
        call jedetr(zk80(jnoloc+iCellType-1))
    end do
    AS_DEALLOCATE(vk8=asterCellType)
    AS_DEALLOCATE(vi=MEDCellType)
    call jedetr('&&LRMPGA_TYPGEO_NBPG_MED')
    call jedetr('&&LRMPGA_TYPGEO_OKPG_MED')
    call jedetr('&&LRMPGA_TYPGEO_NOMLOC')
!
    call jedema()
!
    if (nivinf .gt. 1) then
        write (ifm, 101) 'FIN DE '//nompro
    end if
!
101 format(/, 10('='), a, 10('='),/)
201 format('POUR LE CHAMP MED ', a,&
    &     /, 'MAILLE ! NBRE DE PTS DE GAUSS',&
    &       ' ! ELREFE ASTER ASSOCIE ! NOM LOCALISATION',&
    &     /, 72('-'))
!
end subroutine
