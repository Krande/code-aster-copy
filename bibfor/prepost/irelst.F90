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
subroutine irelst(nofimd, chanom, nochmd, typech, nomaas, &
                  nomamd, nbimpr, caimpi, caimpk, sdcarm, &
                  carael)
!
    use as_med_module, only: as_med_open
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/as_mmhcow.h"
#include "asterfort/as_mmhcyw.h"
#include "asterfort/as_msecre.h"
#include "asterfort/as_msense.h"
#include "asterfort/as_msesei.h"
#include "asterfort/as_msevac.h"
#include "asterfort/as_msmcre.h"
#include "asterfort/as_msmsmi.h"
#include "asterfort/assert.h"
#include "asterfort/elref2.h"
#include "asterfort/irmaes.h"
#include "asterfort/jedetr.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/lrmtyp.h"
#include "asterfort/uteref.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8) :: nomaas, typech, sdcarm, carael
    character(len=*) :: nofimd
    character(len=19) :: chanom
    character(len=64) :: nomamd, nochmd
    integer(kind=8) :: nbimpr, caimpi(10, nbimpr)
    character(len=80) :: caimpk(3, nbimpr)
!
! --------------------------------------------------------------------------------------------------
!
!  IMPR_RESU - IMPRESSION DES ELEMENTS DE STRUCTURE AU FORMAT MED
!
! --------------------------------------------------------------------------------------------------
!
! IN  :
!   NOFIMD  K*   NOM DU FICHIER MED
!   CHANOM  K19  NOM DU CHAMP A IMPRIMER
!   TYPECH  K8   TYPE DU CHAMP
!   NOMAAS  K8   NOM DU MAILLAGE ASTER A COMPLETER DANS LE FICHIER MED
!   NOMAMD  K*   NOM DU MAILLAGE MED
!   NBIMPR  I    NOMBRE D'IMPRESSIONS
!   CAIMPI  I*   ENTIERS POUR CHAQUE IMPRESSION
!   CAIMPK  K80* CARACTERES POUR CHAQUE IMPRESSION
!   SDCARM  K*   SD_CARA_ELEM EN CHAM_ELEM_S
!
! --------------------------------------------------------------------------------------------------
!
    med_idt :: idfimd
    integer(kind=8), parameter :: lgmax = 1000, edleaj = 1
    integer(kind=8), parameter :: edmail = 0, edcart = 0, edfuin = 0, ednoda = 0, edtyre = 6
    !
    integer(kind=8) :: inimpr, nbcouc, nbsect, nummai, ntypef, codret
    integer(kind=8) :: nbnoso, nbnoto, nbrepg, ndim, nbelr, nbgamm, nbpg, nbtyp, nbsp
    integer(kind=8) :: medcel, nbmssu, nbattc, prespr, ino, inimp2
    integer(kind=8) :: tymaas, tymamd, connex(9)
   integer(kind=8) :: imasup, jmasup, nbmasu, nbmsmx, nvtymd, edcar2, nbattv, dimest, nbnosu, tygems
    !
    integer(kind=8)           :: nnotyp(MT_NTYMAX), typgeo(MT_NTYMAX), renumd(MT_NTYMAX), modnum(MT_NTYMAX)
    integer(kind=8)           :: nuanom(MT_NTYMAX, MT_NNOMAX), numnoa(MT_NTYMAX, MT_NNOMAX)
    character(len=8)  :: lielrf(MT_NBFAMX)
    character(len=8)  :: nomtyp(MT_NTYMAX)
    !
    real(kind=8)      :: refcoo(3*lgmax), gscoo(3*lgmax), wg(lgmax)
    !
    character(len=8)  :: saux08
    character(len=16) :: nomtef, nomfpg, nocoo2(3), uncoo2(3), valk(3)
    character(len=64) :: nomasu, nomaes, nomas2
    !
    character(len=64), parameter :: atepai = 'SCALE', atangv = 'ANGLE'
    character(len=64), parameter :: nocoqu = 'SP_SHELL', nompmf = 'SP_BEAM'
    character(len=64), parameter :: notuya = 'SP_PIPE'
    !
    character(len=200) :: desmed
    !
    aster_logical    :: newest, okgr, okcq, oktu, okpf
    integer(kind=8), pointer :: nv_type_med(:) => null()
    !
    character(len=16), parameter :: nocoor(3) = (/'X               ', &
                                                  'Y               ', &
                                                  'Z               '/)
    character(len=16), parameter :: uncoor(3) = (/'INCONNU         ', &
                                                  'INCONNU         ', &
                                                  'INCONNU         '/)

!
! --------------------------------------------------------------------------------------------------
!
    call as_med_open(idfimd, nofimd, edleaj, codret)
    if (codret .ne. 0) then
        saux08 = 'mfiope'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
    !
    !   RELECTURE DES ELEMENTS DE STRUCTURES DEJA PRESENTS
    nbmasu = 0
    call as_msense(idfimd, nbmasu, codret)
    if (codret .ne. 0) then
        saux08 = 'msmnsm'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
    nbmsmx = nbmasu+10
    call wkvect('&&IRELST.MAIL_SUPP', 'V V K80', nbmsmx, jmasup)
    AS_ALLOCATE(vi=nv_type_med, size=nbmsmx)
    if (nbmasu .ne. 0) then
        do imasup = 1, nbmasu
            call as_msmsmi(idfimd, imasup, nomasu, ndim, desmed, &
                           edcar2, nocoo2, uncoo2, codret)
            if (codret .ne. 0) then
                saux08 = 'msmsmi'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
            zk80(jmasup+imasup-1) = nomasu
            !
            call as_msesei(idfimd, imasup, nomaes, nvtymd, dimest, &
                           nomasu, medcel, nbnosu, nbmssu, tygems, &
                           nbattc, prespr, nbattv, codret)
            if (codret .ne. 0) then
                saux08 = 'msesei'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
            nv_type_med(imasup) = nvtymd
        end do
    end if
    !
    desmed = ' '
    !
    call lrmtyp(nbtyp, nomtyp, nnotyp, typgeo, renumd, &
                modnum, nuanom, numnoa)
    !
    ! CREATION DES ELEMENTS DE STRUCTURES DANS LE FICHIER MED
    ! UN ELEMENT DE STRUCTURE EST DEFINIT PAR UNE PAIRE :
    !   TYPE ELEMENT (COQUE, TUYAU, ...) + TYPE MAILLE
    newest = ASTER_FALSE
    cinimpr: do inimpr = 1, nbimpr
        ntypef = caimpi(1, inimpr)
        nbpg = caimpi(2, inimpr)
        nbsp = caimpi(3, inimpr)
        nbcouc = caimpi(4, inimpr)
        nbsect = caimpi(5, inimpr)
        nummai = caimpi(6, inimpr)
        tymaas = caimpi(8, inimpr)
        tymamd = caimpi(9, inimpr)
        !
        call jenuno(jexnum('&CATA.TE.NOMTE', ntypef), nomtef)
        !
        call elref2(nomtef, MT_NBFAMX, lielrf, nbelr)
        ASSERT(nbelr .gt. 0)
        !
        call uteref(chanom, typech, ntypef, nomtef, .false._1, &
                    nomfpg, nbnoso, nbnoto, nbrepg, ndim, refcoo, &
                    gscoo, wg, nochmd, codret)
        !
        nomasu = ' '
        nomas2 = ' '
        ! CAS : GRILLE, COQUE, TUYAU, PMF
        okgr = (nummai .eq. 0) .and. (nbcouc .eq. 1) .and. (nbsect .eq. 0) .and. (nbsp .eq. 1)
        okcq = (nummai .eq. 0) .and. (nbcouc .ge. 1) .and. (nbsect .eq. 0) .and. &
               (nbsp .eq. 3*nbcouc)
        oktu = (nummai .eq. 0) .and. (nbcouc .ge. 1) .and. (nbsect .ge. 1)
        okpf = (nummai .ne. 0) .and. (nbcouc .eq. 0) .and. (nbsect .eq. 0)
        if (okgr) then
            ! CAS D'UNE GRILLE            : nocoqu=SP_SHELL   (8)
            !   Il faut distinguer les éléments supports TRIA3 et QUA4
            nomasu(1:8) = nocoqu
            nomas2(1:8) = nocoqu
            if (nomfpg(1:3) .eq. 'TR3') then
                nomas2(9:10) = '_3'
            else if (nomfpg(1:3) .eq. 'QU4') then
                nomas2(9:10) = '_4'
            else
                valk(1) = 'GRILLE'
                valk(2) = nomfpg(1:8)
                valk(3) = nocoqu(1:16)
                call utmess('F', 'MED2_20', nk=3, valk=valk)
            end if
        else if (okcq) then
            ! CAS D'UNE COQUE               : nocoqu=SP_SHELL   (8)
            !   Il faut distinguer les éléments supports TRIA3 et QUA4
            nomasu(1:8) = nocoqu
            nomas2(1:8) = nocoqu
            if (nomfpg(1:3) .eq. 'TR3') then
                nomas2(9:10) = '_3'
            else if (nomfpg(1:3) .eq. 'QU4') then
                nomas2(9:10) = '_4'
            else
                valk(1) = 'COQUE'
                valk(2) = nomfpg(1:8)
                valk(3) = nocoqu(1:16)
                call utmess('F', 'MED2_20', nk=3, valk=valk)
            end if
        else if (oktu) then
            ! CAS D'UN TUYAU                : notuya=SP_PIPE    (7)
            !   Il faut distinguer les éléments supports SE3 et SE4
            nomasu(1:7) = notuya
            nomas2(1:7) = notuya
            nbgamm = 0
            if (nomfpg(1:3) .eq. 'SE3') then
                nbgamm = 3
                nomas2(8:9) = '_3'
            else if (nomfpg(1:3) .eq. 'SE4') then
                nbgamm = 4
                nomas2(8:9) = '_4'
            else
                valk(1) = 'TUYAU'
                valk(2) = nomfpg(1:8)
                valk(3) = notuya(1:16)
                call utmess('F', 'MED2_20', nk=3, valk=valk)
            end if
        else if (okpf) then
            ! CAS D'UNE PMF                 : nompmf= SP_BEAM   (7)
            nomasu(1:7) = nompmf
            nomas2(1:7) = nompmf
        else
            cycle cinimpr
        end if
        nomasu(9:11) = nomfpg(1:3)
        do inimp2 = 1, nbimpr
            if (caimpk(3, inimp2) .eq. nomasu) then
                caimpk(3, inimpr) = nomasu
                caimpi(9, inimpr) = caimpi(9, inimp2)
                cycle cinimpr
            end if
        end do
        do imasup = 1, nbmasu
            if (zk80(jmasup+imasup-1) .eq. nomasu) then
                caimpk(3, inimpr) = zk80(jmasup+imasup-1)
                caimpi(9, inimpr) = nv_type_med(imasup)
                cycle cinimpr
            end if
        end do
        ! DEFINITION DU MAILLAGE SUPPORT MED
        call as_msmcre(idfimd, nomas2, ndim, desmed, edcart, nocoor, uncoor, codret)
        if (codret .ne. 0) then
            saux08 = 'msmcre'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
        ! DEFINITION DES NOEUDS DU MAILLAGE SUPPORT MED
        call as_mmhcow(idfimd, nomas2, refcoo, edfuin, nbnoto, codret)
        if (codret .ne. 0) then
            saux08 = 'mmhcow'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
        ! CREATION DE LA CONNECTIVITE
        ASSERT(nbnoto .le. 9)
        if (modnum(tymaas) .eq. 0) then
            do ino = 1, nbnoto
                connex(ino) = ino
            end do
        else
            do ino = 1, nbnoto
                connex(ino) = nuanom(tymaas, ino)
            end do
        end if
        ! DEFINITION DE LA MAILLE DU MAILLAGE SUPPORT
        call as_mmhcyw(idfimd, nomas2, connex, nbnoto, edfuin, 1, edmail, tymamd, ednoda, codret)
        if (codret .ne. 0) then
            saux08 = 'mmhcyw'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
        ! SAUVEGARDE DE L'ELEMENT DE STRUCTURE
        nbmasu = nbmasu+1
        if (nbmasu .gt. nbmsmx) then
            nbmsmx = nbmsmx+10
            call juveca('&&IRELST.MAIL_SUPP', nbmsmx)
            call jeveuo('&&IRELST.MAIL_SUPP', 'E', jmasup)
        end if
        zk80(jmasup+nbmasu-1) = nomasu
        !
        nvtymd = -9999
        call as_msecre(idfimd, nomasu, ndim, nomas2, edmail, tymamd, nvtymd, codret)
        ASSERT(nvtymd .ne. -9999)
        if (codret .ne. 0) then
            saux08 = 'msecre'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
        !
        if (nomasu(1:8) .eq. nocoqu) then
            ! 1 attribut variable dimension 2 : epaisseur excentrement
            call as_msevac(idfimd, nomasu, atepai, edtyre, 2, codret)
            if (codret .ne. 0) then
                saux08 = 'msevac'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
        else if (nomasu(1:7) .eq. notuya) then
            ! 1 attribut variable dimension 2 : Rmin, Rmax
            call as_msevac(idfimd, nomasu, atepai, edtyre, 2, codret)
            if (codret .ne. 0) then
                saux08 = 'msevac'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
            ! 1 attribut variable dimension 3 ou 4 : angle
            call as_msevac(idfimd, nomasu, atangv, edtyre, nbgamm, codret)
            if (codret .ne. 0) then
                saux08 = 'msevac'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
        else if (nomasu(1:7) .eq. nompmf) then
            ! 1 attribut variable dimension 1 : angle
            call as_msevac(idfimd, nomasu, atangv, edtyre, 1, codret)
            if (codret .ne. 0) then
                saux08 = 'msevac'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
        !
        ! MODIFICATION DU TYPE MED A IMPRIMER
        caimpi(9, inimpr) = nvtymd
        caimpk(3, inimpr) = nomasu
        newest = ASTER_TRUE
        !
    end do cinimpr
    ! AJOUT DES MAILLES "STRUCTURES" AU MAILLAGE
    if (newest) then
        call irmaes(idfimd, nomaas, nomamd, nbimpr, caimpi, &
                    modnum, nuanom, nomtyp, nnotyp, sdcarm, &
                    carael)
    end if
    !
    call as_mficlo(idfimd, codret)
    if (codret .ne. 0) then
        saux08 = 'mficlo'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
    !
    call jedetr('&&IRELST.MAIL_SUPP')
    AS_DEALLOCATE(vi=nv_type_med)
!
end subroutine
