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
subroutine ccpoux(postCompPoux, &
                  listLoadZ, modelZ, &
                  resultIn, resultType, numeStore, &
                  nbParaIn, lpain, lchin, &
                  iret)
!
    use postComp_type
    use loadMecaCompute_type
    use loadMecaCompute_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/focste.h"
#include "asterfort/fointe.h"
#include "asterfort/fozero.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "LoadTypes_type.h"
!
    type(POST_COMP_POUX), intent(in) :: postCompPoux
    character(len=*), intent(in) :: modelZ, listLoadZ
    character(len=8), intent(in) :: resultIn
    character(len=16), intent(in) :: resultType
    integer(kind=8), intent(in)  :: numeStore
    integer(kind=8), intent(in) :: nbParaIn
    character(len=8), intent(in) :: lpain(100)
    character(len=24), intent(inout) :: lchin(100)
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP
!
! POUX beams
!
! --------------------------------------------------------------------------------------------------
!
! IN  :
!   RESUIN  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT IN
!   TYPESD  K16  TYPE DE LA STRUCTURE DE DONNEES RESULTAT
!   NORDRE  I    NUMERO D'ORDRE COURANT
!   NBCHRE  I    NOMBRE DE CHARGES REPARTIES (POUTRES)
!   IOCCUR  I    NUMERO D'OCCURENCE OU SE TROUVE LE CHARGE REPARTIE
!   KCHARG  K19  NOM DE L'OBJET JEVEUX CONTENANT LES CHARGES
!   MODELE  K8   NOM DU MODELE
!   NBPAIN  I    NOMBRE DE PARAMETRES IN
!   LIPAIN  K8*  LISTE DES PARAMETRES IN
!   LICHIN  K8*  LISTE DES CHAMPS IN
!   SUROPT  K24
!
! OUT :
!   IRET    I    CODE RETOUR (0 SI OK, 1 SINON)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbCmpForc = 11, nbCmpGrav = 4
    real(kind=8) :: cmpZeroReal(nbCmpForc)
    complex(kind=8) :: cmpZeroCplx(nbCmpForc)
    character(len=8) :: cmpZeroFunc(nbCmpForc)
    character(len=8), parameter :: cmpName(nbCmpForc) = (/'FX   ', 'FY   ', 'FZ   ', &
                                                          'MX   ', 'MY   ', 'MZ   ', &
                                                          'BX   ', 'REP  ', 'ALPHA', &
                                                          'BETA ', 'GAMMA'/)
    character(len=8), parameter :: ncmppe(nbCmpGrav) = (/'G ', 'AG', 'BG', 'CG'/)
    character(len=16), parameter :: loadKeyword = "EXCIT"
    aster_logical :: exif1d, exigrav, hasLoad
    integer(kind=8) :: jvPara, jvMGamma, jvAcce, iEqua, nbEqua
    integer(kind=8) :: funcFromUser, coefFromUser, loadFromUser, iParaIn, ier
    real(kind=8) :: coefReal, funcVale, omega2
    real(kind=8) :: freq, inst
    complex(kind=8) :: coefCplx
    character(len=1) :: typeScal
    character(len=5) :: ch5
    character(len=6) :: scalType
    character(len=8) :: paraCurr, loadName, loadCommand
    character(len=8) :: loadFunc, physQuanName
    character(len=13) :: loadPreObject
    character(len=16) :: modeType
    character(len=19) :: chacce
    character(len=24) :: fieldInName, chdepl, modelLigrel
    character(len=24) :: loadFieldForce, loadFieldGrav
    character(len=8), pointer :: listLoadFunc(:) => null()
    real(kind=8), pointer :: dispValeR(:) => null()
    complex(kind=8), pointer :: dispValeC(:) => null()
    character(len=8), pointer :: listLoadName(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    iret = 0

! - Initializations
    hasLoad = postCompPoux%loadIndx .ne. 0
    cmpZeroReal(1:nbCmpForc) = 0.d0
    cmpZeroFunc(1:nbCmpForc) = '&FOZERO'
    cmpZeroCplx(1:nbCmpForc) = (0.d0, 0.d0)
    call dismoi("NOM_LIGREL", modelZ, "MODELE", repk=modelLigrel)

! - Type of base
    modeType = ' '
    if (resultType .eq. 'MODE_MECA') then
        call rsadpa(resultIn, 'L', 1, 'TYPE_MODE', 1, 0, sjv=jvPara)
        modeType = zk16(jvPara)
    end if

! - Compute M.Gamma
    coefReal = 0.d0
    coefCplx = (0.d0, 0.d0)
    if ((resultType .eq. 'MODE_MECA' .and. modeType(1:8) .eq. 'MODE_DYN') .or. &
        (resultType .eq. 'MODE_ACOU')) then

! ----- Get displacements
        call rsexch('F', resultIn, 'DEPL', numeStore, chdepl, ier)
        call jelira(chdepl(1:19)//'.VALE', 'LONMAX', nbEqua)

! ----- Get pulsation
        call rsadpa(resultIn, 'L', 1, 'OMEGA2', numeStore, 0, sjv=jvPara)
        omega2 = zr(jvPara)

! ----- Compute M.Gamma
        call dismoi('NOM_GD', postCompPoux%chdynr, 'CHAMP', repk=physQuanName)
        call dismoi('TYPE_SCA', physQuanName, 'GRANDEUR', repk=scalType)
        call jeveuo(postCompPoux%chdynr(1:19)//'.VALE', 'E', jvMGamma)

        if (scalType .eq. 'R') then
            call jeveuo(chdepl(1:19)//'.VALE', 'L', vr=dispValeR)
            do iEqua = 1, nbEqua
                zr(jvMGamma-1+iEqua) = -omega2*dispValeR(iEqua)
            end do
        elseif (scalType .eq. 'C') then
            call jeveuo(chdepl(1:19)//'.VALE', 'L', vc=dispValeC)
            do iEqua = 1, nbEqua
                zc(jvMGamma-1+iEqua) = -omega2*dispValeC(iEqua)
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
        call jelibe(chdepl(1:19)//'.VALE')

    else if (resultType .eq. 'DYNA_TRANS') then
! ----- Get displacements
        call rsexch('F', resultIn, 'DEPL', numeStore, chdepl, ier)
        call jelira(chdepl(1:19)//'.VALE', 'LONMAX', nbEqua)

! ----- Get acceleration
        call rsexch(' ', resultIn, 'ACCE', numeStore, chacce, ier)

! ----- Compute M.Gamma
        call jeveuo(postCompPoux%chdynr(1:19)//'.VALE', 'E', jvMGamma)
        if (ier .eq. 0) then
            call jeveuo(chacce//'.VALE', 'L', jvAcce)
            do iEqua = 1, nbEqua
                zr(jvMGamma-1+iEqua) = zr(jvAcce-1+iEqua)
            end do
            call jelibe(chacce//'.VALE')
        else
            call utmess('A', 'CALCULEL3_1')
            do iEqua = 1, nbEqua
                zr(jvMGamma-1+iEqua) = 0.d0
            end do
        end if

    else if (resultType .eq. 'DYNA_HARMO') then
! ----- Get displacements
        call rsexch('F', resultIn, 'DEPL', numeStore, chdepl, ier)
        call jelira(chdepl(1:19)//'.VALE', 'LONMAX', nbEqua)

! ----- Access to M.Gamma vector
        call jeveuo(postCompPoux%chdynr(1:19)//'.VALE', 'E', jvMGamma)

! ----- Get acceleration
        call rsexch(' ', resultIn, 'ACCE', numeStore, chacce, ier)

! ----- Compute M.Gamma
        if (ier .eq. 0) then
            call jeveuo(chacce//'.VALE', 'L', jvAcce)
            do iEqua = 1, nbEqua
                zc(jvMGamma+iEqua) = zc(jvAcce+iEqua)
            end do
            call jelibe(chacce//'.VALE')
        else
            call utmess('A', 'CALCULEL3_1')
            do iEqua = 1, nbEqua
                zc(jvMGamma+iEqua) = (0.d0, 0.d0)
            end do
        end if
    end if

! - Get coefficient for load
    coefReal = 0.0d0
    coefCplx = (0.d0, 0.d0)
    loadName = " "
    loadFunc = " "
    loadFieldForce = " "
    exif1d = ASTER_FALSE
    loadFieldGrav = " "
    exigrav = ASTER_FALSE
    typeScal = ' '

    if (hasLoad) then
! ----- Get load from user
        loadName = ' '
        loadFunc = ' '
        call getvid(loadKeyword, 'CHARGE', iocc=postCompPoux%loadIndx, &
                    scal=loadName, nbret=loadFromUser)
        call getvid(loadKeyword, 'FONC_MULT', iocc=postCompPoux%loadIndx, &
                    scal=loadFunc, nbret=funcFromUser)
        call getvr8(loadKeyword, 'COEF_MULT', iocc=postCompPoux%loadIndx, &
                    scal=coefReal, nbret=coefFromUser)

        if (loadFromUser .eq. 0) then
! --------- Get load from result
            call jeveuo(listLoadZ(1:19)//'.LCHA', 'L', vk8=listLoadName)
            call jeveuo(listLoadZ(1:19)//'.FCHA', 'L', vk8=listLoadFunc)
            loadName = listLoadName(postCompPoux%loadIndx)
            loadFunc = listLoadFunc(postCompPoux%loadIndx)
            if (loadFunc(1:2) .eq. '&&') then
                ASSERT(loadFunc .eq. '&&NMDOME')
                coefReal = 1.d0
                call focste(loadFunc, 'TOUTRESU', coefReal, 'V')
            end if
            funcFromUser = 1
            coefFromUser = 0
        else
            if (funcFromUser+coefFromUser .ne. 0) then
                if ((resultType .ne. 'DYNA_HARMO') .and. (resultType .ne. 'DYNA_TRANS') .and. &
                    (resultType .ne. 'EVOL_ELAS')) then
                    call utmess('A', 'CALCULEL3_4')
                    iret = 1
                    goto 999
                end if
            end if
        end if

! ----- Get coefficient from user
        if (funcFromUser .ne. 0 .or. coefFromUser .ne. 0) then
            if (resultType .eq. 'DYNA_HARMO') then
                typeScal = 'C'
                call rsadpa(resultIn, 'L', 1, 'FREQ', numeStore, 0, sjv=jvPara)
                freq = zr(jvPara)
                if (funcFromUser .ne. 0) then
                    call fointe('F ', loadFunc, 1, ['FREQ'], [freq], funcVale, ier)
                    coefCplx = dcmplx(funcVale, 0.d0)
                else if (coefFromUser .ne. 0) then
                    coefCplx = dcmplx(coefReal, 1.d0)
                end if
            else if (resultType .eq. 'DYNA_TRANS') then
                typeScal = 'R'
                call rsadpa(resultIn, 'L', 1, 'INST', numeStore, 0, sjv=jvPara)
                inst = zr(jvPara)
                if (funcFromUser .ne. 0) then
                    call fointe('F ', loadFunc, 1, ['INST'], [inst], funcVale, ier)
                    coefReal = funcVale
                else if (coefFromUser .ne. 0) then

                else
                    call utmess('A', 'CALCULEL3_2')
                    iret = 1
                    goto 999
                end if
            else if (resultType .eq. 'EVOL_ELAS') then
                typeScal = 'R'
                call rsadpa(resultIn, 'L', 1, 'INST', numeStore, 0, sjv=jvPara)
                inst = zr(jvPara)
                if (funcFromUser .ne. 0) then
                    call fointe('F ', loadFunc, 1, ['INST'], [inst], funcVale, ier)
                    coefReal = funcVale
                else
                    call utmess('A', 'CALCULEL3_3')
                    iret = 1
                    goto 999
                end if
            end if
        end if

! ----- Detect FORCE_POUTRE/PESANTEUR
        call dismoi('TYPE_CHARGE', loadName, 'CHARGE', repk=loadCommand)
        loadPreObject = loadName(1:8)//'.CHME'
        call getMecaNeumField(LOAD_NEUM_FORC_BEAM, loadPreObject, loadFieldForce)
        call exisd('CHAMP_GD', loadFieldForce, ier)
        exif1d = ier .ne. 0
        call getMecaNeumField(LOAD_NEUM_GRAVITY, loadPreObject, loadFieldGrav)
        call exisd('CHAMP_GD', loadFieldGrav, ier)
        exigrav = ier .ne. 0
    end if

    ch5 = '.    '
    do iParaIn = 1, nbParaIn
        paraCurr = lpain(iParaIn)
        ch5 = '.    '

! ----- Create impedance input field (real)
        if ((paraCurr .eq. 'PCOEFFR') .and. (typeScal .eq. 'R')) then
            fieldInName = '&&MECHPO.    .COEFF'
            call mecact('V', fieldInName, 'MODELE', modelLigrel, 'IMPE_R', &
                        ncmp=1, nomcmp='IMPE', sr=coefReal)
            lchin(iParaIn) = fieldInName
        end if

! ----- Create impedance input field (complexe)
        if ((paraCurr .eq. 'PCOEFFC') .and. (typeScal .eq. 'C')) then
            fieldInName = '&&MECHPO.    .COEFF'
            call mecact('V', fieldInName, 'MODELE', modelLigrel, 'IMPE_C', &
                        ncmp=1, nomcmp='IMPE', sc=coefCplx)
            lchin(iParaIn) = fieldInName
        end if

! ----- Create load input field (gravity) - Set zero if no load
        if (paraCurr .eq. 'PPESANR') then
            if (exigrav) then
                lchin(iParaIn) = loadFieldGrav
            else
                call codent(iParaIn, 'D0', ch5(2:5))
                fieldInName = '&&MECHPO'//ch5//'.PESAN.DESC'
                lchin(iParaIn) = fieldInName
                call mecact('V', fieldInName, 'MODELE', modelLigrel, 'PESA_R  ', &
                            ncmp=nbCmpGrav, lnomcmp=ncmppe, vr=cmpZeroReal)
            end if
        end if

! ----- Create load input field (beam force, function) - Set zero if no load
        if (paraCurr .eq. 'PFF1D1D') then
            if (exif1d) then
                if (loadCommand(5:7) .eq. '_FO') then
                    lchin(iParaIn) = loadFieldForce
                else
                    call codent(iParaIn, 'D0', ch5(2:5))
                    fieldInName = '&&MECHPO'//ch5//'.P1D1D'
                    lchin(iParaIn) = fieldInName
                    call fozero(cmpZeroFunc(1))
                    call mecact('V', fieldInName, 'MODELE', modelLigrel, 'FORC_F  ', &
                                ncmp=nbCmpForc, lnomcmp=cmpName, vk=cmpZeroFunc)
                end if
            else
                call codent(iParaIn, 'D0', ch5(2:5))
                fieldInName = '&&MECHPO'//ch5//'.P1D1D'
                lchin(iParaIn) = fieldInName
                call fozero(cmpZeroFunc(1))
                call mecact('V', fieldInName, 'MODELE', modelLigrel, 'FORC_F  ', &
                            ncmp=nbCmpForc, lnomcmp=cmpName, vk=cmpZeroFunc)
            end if
        end if

! ----- Create load input field (beam force, real) - Set zero if no load
        if (paraCurr .eq. 'PFR1D1D') then
            if (exif1d) then
                if ((loadCommand(5:7) .eq. '_FO') .or. (loadCommand(5:7) .eq. '_RI')) then
                    call codent(iParaIn, 'D0', ch5(2:5))
                    fieldInName = '&&MECHPO'//ch5//'.P1D1D'
                    lchin(iParaIn) = fieldInName
                    call mecact('V', fieldInName, 'MODELE', modelLigrel, 'FORC_R  ', &
                                ncmp=nbCmpForc, lnomcmp=cmpName, vr=cmpZeroReal)
                else
                    lchin(iParaIn) = loadFieldForce
                end if
            else
                call codent(iParaIn, 'D0', ch5(2:5))
                fieldInName = '&&MECHPO'//ch5//'.P1D1D'
                lchin(iParaIn) = fieldInName
                call mecact('V', fieldInName, 'MODELE', modelLigrel, 'FORC_R  ', &
                            ncmp=nbCmpForc, lnomcmp=cmpName, vr=cmpZeroReal)
            end if
        end if

! ----- Create load input field (beam force, complex)  - Set zero if no load
        if (paraCurr .eq. 'PFC1D1D') then
            if (exif1d) then
                if (loadCommand(5:7) .eq. '_RI') then
                    lchin(iParaIn) = loadFieldForce
                else
                    call codent(iParaIn, 'D0', ch5(2:5))
                    fieldInName = '&&MECHPO'//ch5//'.P1D1D'
                    lchin(iParaIn) = fieldInName
                    call mecact('V', fieldInName, 'MODELE', modelLigrel, 'FORC_C  ', &
                                ncmp=nbCmpForc, lnomcmp=cmpName, vc=cmpZeroCplx)
                end if
            else
                call codent(iParaIn, 'D0', ch5(2:5))
                fieldInName = '&&MECHPO'//ch5//'.P1D1D'
                lchin(iParaIn) = fieldInName
                call mecact('V', fieldInName, 'MODELE', modelLigrel, 'FORC_C  ', &
                            ncmp=nbCmpForc, lnomcmp=cmpName, vc=cmpZeroCplx)

            end if
        end if

! ----- Create load input field (M.Gamma force)
        if (paraCurr .eq. 'PCHDYNR') then
            fieldInName = postCompPoux%chdynr(1:19)//'.VALE'
            call jeexin(fieldInName, ier)
            if (ier .eq. 0) then
                call codent(iParaIn, 'D0', ch5(2:5))
                fieldInName = '&&MECHPO'//ch5//'.PCHDY'
                call rsexch('F', resultIn, 'DEPL', numeStore, chdepl, ier)
                call copisd('CHAMP_GD', 'V', chdepl, fieldInName)
            end if
            lchin(iParaIn) = fieldInName
        end if

! ----- Create option input field
        if (paraCurr .eq. 'PSUROPT') then
            call codent(iParaIn, 'D0', ch5(2:5))
            fieldInName = '&&MECHPO'//ch5//'.SUR_OPTION'
            lchin(iParaIn) = fieldInName
            call mecact('V', fieldInName, 'MODELE', modelLigrel, 'NEUT_K24', &
                        ncmp=1, nomcmp='Z1', sk=postCompPoux%optionMass)
        end if
    end do
!
999 continue
!
    call jedema()
end subroutine
