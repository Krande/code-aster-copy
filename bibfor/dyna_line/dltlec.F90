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
! aslint: disable=W1504
!
subroutine dltlec(result, model, numeDOF, materField, mate, &
                  caraElem, jvMatr, masse, rigid, &
                  amort, lamort, nbLoad, nbVectAsse, listLoad, &
                  loadNameJv, loadInfoJv, loadFuncJv, jvVectAsse, jvVectFunc, &
                  nbWave, jvLoadWave, solveu, iinteg, t0, &
                  nume, numrep, ds_inout)
!
    use listLoad_type
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/chpver.h"
#include "asterfort/codent.h"
#include "asterfort/cresol.h"
#include "asterfort/dismoi.h"
#include "asterfort/dltp0.h"
#include "asterfort/focste.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtdscr.h"
#include "asterfort/nmarnr.h"
#include "asterfort/nmdoch.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/nonlinDSInOutCreate.h"
!
    character(len=19), intent(out) :: listLoad, solveu
    character(len=24), intent(out) :: loadNameJv, loadInfoJv, loadFuncJv
    integer(kind=8), intent(out) :: nbVectAsse, nbLoad, nbWave
    integer(kind=8), intent(out) :: jvLoadWave, jvVectAsse, jvVectFunc
    character(len=8), intent(out) :: result
    character(len=8), intent(out) :: masse, rigid, amort
    aster_logical, intent(out) :: lamort
    integer(kind=8), intent(out) :: jvMatr(3)
    character(len=24), intent(out) :: model, numeDOF, mate, caraElem, materField
    type(NL_DS_InOut), intent(out) :: ds_inout
    integer(kind=8), intent(out) :: nume, numrep, iinteg
    real(kind=8), intent(out) :: t0
!
! --------------------------------------------------------------------------------------------------
!
!       DYNAMIQUE LINEAIRE TRANSITOIRE - LECTURE DES DONNEES
!
! --------------------------------------------------------------------------------------------------
!
!      OUT RESULT : NOM UTILISATEUR DU RESULTAT DE STAT_NON_LINE
!      OUT MODELE : NOM DU MODELE
!      OUT NUMEDD : NUME_DDL DE LA MATR_ASSE RIGID
!      OUT MATERI : NOM DU CHAMP DE MATERIAU
!      OUT MATE   : NOM DU CHAMP DE MATERIAU CODE
!      OUT CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
!      OUT MASSE  : MATRICE DE MASSE
!      OUT RIGID  : MATRICE DE RIGIDITE
!      OUT AMORT  : MATRICE D'AMORTISSEMENT
!      OUT LAMORT : LOGIQUE INDIQUANT SI IL Y A AMORTISSEMENT
!      OUT IMAT   : TABLEAU D'ADRESSES POUR LES MATRICES
!      OUT NCHAR  : NOMBRE D'OCCURENCES DU MOT CLE CHARGE
!      OUT NVECA  : NOMBRE D'OCCURENCES DU MOT CLE VECT_ASSE
!      OUT LISCHA : INFO SUR LES CHARGES
!      OUT CHARGE : LISTE DES CHARGES
!      OUT INFOCH : INFO SUR LES CHARGES
!      OUT FOMULT : LISTE DES FONC_MULT ASSOCIES A DES CHARGES
!      OUT IAADVE : ADRESSE
!      OUT IAADVE : ADRESSE
!      OUT NONDP  : NOMBRE D'ONDES PLANES
!      OUT IONDP  : ADRESSE
!      OUT SOLVEU : NOM DU SOLVEUR
!      OUT IINTEG : TYPE D'INTEGRATION
!                   1 : NEWMARK
!                   2 : WILSON
!                   3 : DIFF_CENTRE
!                   4 : ADAPT
!      OUT T0     : INSTANT INITIAL
!      OUT NUME   : NUMERO D'ORDRE DE REPRISE
!      OUT NUMREP : NUMERO DE REUSE POUR LA TABLE PARA_CALC
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jvEner
    character(len=1), parameter :: jvBase = "V"
    character(len=24) :: loadType
    integer(kind=8) :: nbLoadKeyword, iLoadKeyword, n1, nbCine
    integer(kind=8) :: iaux, ibid, nbRet
    integer(kind=8) :: indic, iLoad
    real(kind=8) :: rval
    character(len=16) :: method, k16bid
    character(len=19) :: channo

    type(ListLoad_Prep) :: listLoadPrep
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    character(len=24), pointer :: listLoadName(:) => null(), listLoadFunc(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    solveu = '&&COMDLT.SOLVEUR'
    listLoad = '&&COMDLT.LISTLOAD'
    model = ' '
    materField = ' '
    caraElem = ' '
    jvLoadWave = 0
    jvVectAsse = 0
    jvVectFunc = 0

!
!====
! 2. LES DONNEES DU CALCUL
!====
!
! 2.1. ==> LE CONCEPT RESULTAT CREE PAR LA COMMANDE
!
    call getres(result, k16bid, k16bid)
!
! 2.3. ==> CALCUL DES ENERGIES
!
    call wkvect('&&COMDLT.ENER      .VALE', 'V V R', 6, jvEner)
!
! - Get matrices
    rigid = ' '
    masse = ' '
    amort = ' '
    call getvid(' ', 'MATR_RIGI', scal=rigid, nbret=nbRet)
    call getvid(' ', 'MATR_MASS', scal=masse, nbret=nbRet)
    call getvid(' ', 'MATR_AMOR', scal=amort, nbret=nbRet)
    lamort = nbRet .gt. 0

    jvMatr = 0
    call mtdscr(rigid)
    call jeveuo(rigid//'           .&INT', 'E', jvMatr(1))
    call mtdscr(masse)
    call jeveuo(masse//'           .&INT', 'E', jvMatr(2))
    if (lamort) then
        call mtdscr(amort)
        call jeveuo(amort//'           .&INT', 'E', jvMatr(3))
    end if

! - Get loads
    nbVectAsse = 0
    nbLoad = 0
    nbCine = 0
    call getfac('EXCIT', nbLoadKeyword)
    if (nbLoadKeyword .gt. 0) then
        do iLoadKeyword = 1, nbLoadKeyword
            call getvid('EXCIT', 'VECT_ASSE', iocc=iLoadKeyword, scal=channo, nbret=iaux)
            if (iaux .eq. 1) then
                nbVectAsse = nbVectAsse+1
            end if
            call getvid('EXCIT', 'CHARGE', iocc=iLoadKeyword, scal=channo, nbret=iaux)
            if (iaux .eq. 1) then
                nbLoad = nbLoad+1
            end if
            call gettco(channo, loadType)
            if (loadType .eq. 'CHAR_CINE_MECA') then
                nbCine = nbCine+1
            end if
        end do
!
! 3.1.2. ==> LISTE DE VECT_ASSE DECRIVANT LE CHARGEMENT
!
        if (nbVectAsse .ne. 0) then
            call wkvect('&&COMDLT.LIFONCT', 'V V K24', nbVectAsse, jvVectFunc)
            call wkvect('&&COMDLT.ADVECASS', 'V V I  ', nbVectAsse, jvVectAsse)
            indic = 0
            do iLoadKeyword = 1, nbVectAsse
                indic = indic+1
10              continue
                call getvid('EXCIT', 'VECT_ASSE', iocc=indic, scal=channo, nbret=iaux)
                if (iaux .eq. 0) then
                    indic = indic+1
                    goto 10
                end if
                call chpver('F', channo, 'NOEU', 'DEPL_R', ibid)
                call jeveuo(channo//'.VALE', 'L', zi(jvVectAsse+iLoadKeyword-1))
                call getvid('EXCIT', 'FONC_MULT', iocc=indic, &
                            scal=zk24(jvVectFunc+iLoadKeyword-1), nbret=iaux)
                if (iaux .eq. 0) then
                    call getvid('EXCIT', 'ACCE', iocc=indic, &
                                scal=zk24(jvVectFunc+iLoadKeyword-1), nbret=iaux)
                    if (iaux .eq. 0) then
                        rval = 1.d0
                        call getvr8('EXCIT', 'COEF_MULT', iocc=indic, scal=rval, nbret=iaux)
                        zk24(jvVectFunc+iLoadKeyword-1) = '&&COMDLT.F_'
                        call codent(iLoadKeyword, 'G', zk24(jvVectFunc+iLoadKeyword-1) (12:19))
                        call focste(zk24(jvVectFunc+iLoadKeyword-1), 'INST', rval, 'V')
                    end if
                end if
            end do
        end if

! ----- Get loads/BC and create list of loads datastructure
        if (nbLoad .ne. 0) then
            call getvid(' ', 'MODELE', scal=model, nbret=iaux)
            if (iaux .eq. 0) then
                call utmess('F', 'DYNALINE1_24')
            end if
            listLoadPrep%model = model(1:8)
            call nmdoch(listLoadPrep, listLoad, jvBase)
        end if
    else
        nbVectAsse = 0
        nbLoad = 0
    end if

! - Access to load object
    if (nbLoad .eq. 0) then
        loadNameJv = listLoad(1:19)//'.LCHA'
        loadInfoJv = listLoad(1:19)//'.INFC'
        loadFuncJv = listLoad(1:19)//'.FCHA'
    else
        loadNameJv = listLoad(1:19)//'.LCHA'
        loadInfoJv = listLoad(1:19)//'.INFC'
        loadFuncJv = listLoad(1:19)//'.FCHA'
        call jeveuo(loadNameJv, 'L', vk24=listLoadName)
        call jeveuo(loadInfoJv, 'L', vi=listLoadInfo)
        call jeveuo(loadFuncJv, 'L', vk24=listLoadFunc)
    end if
!
! 3.2. ==> TEST DE LA PRESENCE DE CHARGES DE TYPE 'ONDE_PLANE'
!
    nbWave = 0
    if (nbLoad .ne. 0) then
        do iLoad = 1, nbLoad
            if (listLoadInfo(1+nbLoad+iLoad) .eq. 6) then
                nbWave = nbWave+1
            end if
        end do
    end if
!
! 3.3. ==> RECUPERATION DES DONNEES DE CHARGEMENT PAR ONDE PLANE
!
    if (nbWave .eq. 0) then
        call wkvect('&&COMDLT.ONDP', 'V V K8', 1, jvLoadWave)
    else
        call wkvect('&&COMDLT.ONDP', 'V V K8', nbWave, jvLoadWave)
        indic = 0
        do iLoad = 1, nbLoad
            if (listLoadInfo(1+nbLoad+iLoad) .eq. 6) then
                indic = indic+1
                zk8(jvLoadWave+indic-1) = listLoadName(iLoad) (1:8)
            end if
        end do
    end if

! - Main parameters from rigid matrix
    call dismoi('NOM_NUME_DDL', rigid, 'MATR_ASSE', repk=numeDOF)
    call dismoi('NOM_MODELE', rigid, 'MATR_ASSE', repk=model)

! - Main parameters from command
    if (materField .eq. ' ') call getvid(' ', 'CHAM_MATER', scal=materField, nbret=n1)
    if (caraElem .eq. ' ') call getvid(' ', 'CARA_ELEM', scal=caraElem, nbret=n1)
    if (materField .ne. ' ') then
        call rcmfmc(materField, mate)
    end if
!
! 4.2. ==> LECTURE DES PARAMETRES DU MOT CLE FACTEUR SOLVEUR ---
!
    call cresol(solveu)
!
! 4.3. ==> TYPE D'INTEGRATION
!
    call getvtx('SCHEMA_TEMPS', 'SCHEMA', iocc=1, scal=method, nbret=n1)
!
    if (method .eq. 'NEWMARK') then
        iinteg = 1
    else
        if (method .eq. 'WILSON') then
            iinteg = 2
        else
            if (method .eq. 'DIFF_CENTRE') then
                iinteg = 3
            else
                if (method .eq. 'ADAPT_ORDRE2') then
                    iinteg = 4
                end if
            end if
        end if
    end if
!
! 4.4. ==> L'INSTANT INITIAL ET SON NUMERO D'ORDRE SI REPRISE
!
    call dltp0(t0, nume)

!
! --- RECUPERATION NUMERO REUSE - TABLE PARA_CALC
!
    call nmarnr(result, 'PARA_CALC', numrep)
!
! 4.5. ==> pour OBSERVATION
!
! - Create input/output management datastructure
!
    call nonlinDSInOutCreate('VIBR', ds_inout)

end subroutine
