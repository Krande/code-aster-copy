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
! person_in_charge: daniele.colombo at ifpen.fr
!
subroutine te0588(option, nomte)
!
    use Behaviour_module, only: behaviourOption
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/iselli.h"
#include "asterc/ismaem.h"
#include "asterfort/jevech.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/xasshm.h"
#include "asterfort/xcaehm.h"
#include "asterfort/xfnohm.h"
#include "asterfort/xhmddl.h"
#include "asterfort/xhmini.h"
#include "asterfort/xpeshm.h"
#include "jeveux.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_HM_*, D_PLAN_HM_* for XFEM
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          FORC_NODA, CHAR_MECA_PESA_R
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nno, imatuu, ndim, imate, iinstm, jcret
    integer(kind=8) :: dimmat, npi, npg, li, ibid, yaenrm
    integer(kind=8) :: codret, icodre(1)
    integer(kind=8) :: ipoids, ivf, idfde, igeom
    integer(kind=8) :: iinstp, ideplm, ideplp, icarcr, ipesa
    integer(kind=8) :: icontm, ivarip, ivarim, ivectu, icontp
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5), dimuel
    integer(kind=8) :: dimdef, dimcon, nbvari, nddls, nddlm
    integer(kind=8) :: nmec, np1, np2, nnos
    integer(kind=8) :: nnom
    real(kind=8) :: defgep(13), defgem(13)
    real(kind=8) :: dfdi(20, 3), dfdi2(20, 3)
    real(kind=8) :: drds(25, 11+5), drdsr(25, 11+5), dsde(11+5, 25)
    real(kind=8) :: r(25), sigbar(25), c(25), ck(25), cs(25)
    real(kind=8) :: angnau(3)
    real(kind=8), dimension(:, :), pointer :: work1 => null(), work2 => null(), b => null()
    character(len=3) :: modint
    character(len=8) :: typmod(2)
    character(len=16) :: phenom, elref
    real(kind=8) :: rho(1), rbid(1)
    aster_logical :: axi, fnoevo
    type(THM_DS) :: ds_thm
! =====================================================================
!  CETTE ROUTINE FAIT UN CALCUL EN HM AVEC XFEM
!  25 = (9 DEF MECA) + (9 DEF HEAV MECA) + 4 POUR P1 + 3 pour P1 HEAV
!  16 = 12 MECA + 4 POUR P1
! =====================================================================
!  POUR LES TABLEAUX DEFGEP ET DEFGEM ON A DANS L'ORDRE :
!                                      (PARTIE CLASSIQUE)
!                                      DX DY DZ
!                                      EPXX EPYY EPZZ EPXY EPXZ EPYZ
!                                      PRE1 P1DX P1DY P1DZ
!                                      (PARTIE ENRICHIE)
!                                      H1X  H1Y H1Z  H2X  H2Y  H2Z
!                                      H3X  H3Y H3Z
!                                      H1PRE1  H2PRE1  H3PRE1
!            EPSXY = RAC2/2*(DU/DY+DV/DX)
! =====================================================================
!    POUR LES CHAMPS DE CONTRAINTE
!                                      SIXX SIYY SIZZ SIXY SIXZ SIYZ
!                                      SIPXX SIPYY SIPZZ SIPXY SIPXZ SIPYZ
!                                      M11 FH11X FH11Y FH11Z
!
!        SIXY EST LE VRAI DE LA MECANIQUE DES MILIEUX CONTINUS
!        DANS EQUTHM ON LE MULITPLIERA PAR RAC2
! =====================================================================
!   POUR L'OPTION FORCNODA
!  SI LES TEMPS PLUS ET MOINS SONT PRESENTS
!  C'EST QUE L'ON APPELLE DEPUIS STAT NON LINE  : FNOEVO = VRAI
!  ET ALORS LES TERMES DEPENDANT DE DT SONT EVALUES
!  SI LES TEMPS PLUS ET MOINS NE SONT PAS PRESENTS
!  C'EST QUE L'ON APPELLE DEPUIS CALCNO  : FNOEVO = FAUX
!  ET ALORS LES TERMES DEPENDANT DE DT NE SONT PAS EVALUES
! =====================================================================
! AXI       AXISYMETRIQUE?
! TYPMOD    MODELISATION (D_PLAN, AXI, 3D ?)
! MODINT    METHODE D'INTEGRATION (CLASSIQUE,LUMPEE(D),REDUITE(R) ?)
! NNO       NB DE NOEUDS DE L'ELEMENT
! NNOS      NB DE NOEUDS SOMMETS DE L'ELEMENT
! NNOM      NB DE NOEUDS MILIEUX DE L'ELEMENT
! NDDLS     NB DE DDL SUR LES SOMMETS
! NDDLM     NB DE DDL SUR LES MILIEUX
! NPI       NB DE POINTS D'INTEGRATION DE L'ELEMENT
! NPG       NB DE POINTS DE GAUSS     POUR CLASSIQUE(=NPI)
!                 SOMMETS             POUR LUMPEE   (=NPI=NNOS)
!                 POINTS DE GAUSS     POUR REDUITE  (<NPI)
! NDIM      DIMENSION DE L'ESPACE
! DIMUEL    NB DE DDL TOTAL DE L'ELEMENT
! DIMCON    DIMENSION DES CONTRAINTES GENERALISEES ELEMENTAIRES
! DIMDEF    DIMENSION DES DEFORMATIONS GENERALISEES ELEMENTAIRES
! IVF       FONCTIONS DE FORMES QUADRATIQUES
! =====================================================================
    real(kind=8) :: dt
! =====================================================================
! DECLARATION POUR XFEM
!
    integer(kind=8) :: nfh, nfiss, jfisno, ddlc, contac
    integer(kind=8) :: ddld, ddlm, ddlp, nnop, nnops, nnopm
    integer(kind=8) :: enrmec(3), nenr, dimenr, enrhyd(3)
    integer(kind=8) :: jpintt, jcnset, jheavt, jpmilt, jheavn
    integer(kind=8) :: jlonch, jlst, jstno
    character(len=8) :: enr
    aster_logical :: lVect, lMatr, lVari, lSigm
    character(len=16) :: compor_copy(COMPOR_SIZE)
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: iCompor
    integer(kind=8) :: itabin(2), iSigm, iret
!
! --------------------------------------------------------------------------------------------------
!
    angnau = 0.d0
    imatuu = ismaem()
    ivectu = ismaem()
    icontp = ismaem()
    ivarip = ismaem()
    allocate (work1(11+5, 52*20))
    allocate (work2(25, 52*20))
    allocate (b(25, 52*20))

!
! - Get model of finite element
!
    call thmGetElemModel(ds_thm)
! INITIALISATION POUR XFEM
!
    call xhmini(nomte, nfh, ddld, ddlm, ddlp, nfiss, ddlc, contac)
    call xcaehm(ds_thm, nomte, axi, typmod, modint, &
                mecani, press1, press2, tempe, dimdef, &
                dimcon, nmec, np1, np2, ndim, &
                nno, nnos, nnom, npi, npg, &
                nddls, nddlm, dimuel, ipoids, ivf, &
                idfde, ddld, ddlm, ddlp, enrmec, nenr, &
                dimenr, nnop, nnops, nnopm, enrhyd, ddlc, nfh)
! =====================================================================
! --- PARAMETRES PROPRES A XFEM ---------------------------------------
! =====================================================================
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PLST', 'L', jlst)
    call jevech('PSTANO', 'L', jstno)
    call elref1(elref)
    call teattr('S', 'XFEM', enr, ibid)
    ASSERT(enr(1:2) .eq. 'XH')
    call jevech('PHEA_NO', 'L', jheavn)
!
! PARAMÈTRES PROPRES AUX ÉLÉMENTS 1D ET 2D QUADRATIQUES
!
    if ((ibid .eq. 0) .and. (enr(1:2) .eq. 'XH') .and. .not. iselli(elref)) then
        call jevech('PPMILTO', 'L', jpmilt)
    end if
! PARAMETRE PROPRE AU MULTI-HEAVISIDE
    if (nfiss .gt. 1) then
        call jevech('PFISNO', 'L', jfisno)
    end if
! =====================================================================
! --- DEBUT DES DIFFERENTES OPTIONS -----------------------------------
! =====================================================================
! --- 2. OPTIONS : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA -------------
! =====================================================================
    if ((option(1:9) .eq. 'RIGI_MECA') .or. (option(1:9) .eq. 'RAPH_MECA') .or. &
        (option(1:9) .eq. 'FULL_MECA')) then
! =====================================================================
! --- PARAMETRES EN ENTREE --------------------------------------------
! =====================================================================
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PMATERC', 'L', imate)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PCARCRI', 'L', icarcr)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PCONTMR', 'L', icontm)

! ---- Make copy of COMPOR map
        do iCompor = 1, COMPOR_SIZE
            compor_copy(iCompor) = compor(iCompor)
        end do

! ----- Force DEFO_LDC="MECANIQUE" for THM
        if (option(1:9) .eq. 'RIGI_MECA') then
            compor_copy(DEFO_LDC) = "MECANIQUE"
        end if

        read (compor_copy(NVAR), '(I16)') nbvari
! =====================================================================
! ----RECUPERATION DES ANGLES NAUTIQUES/EULER DEFINIS PAR AFFE_CARA_ELEM
! --- ORIENTATION DU MASSIF
! --- COORDONNEES DU BARYCENTRE ( POUR LE REPRE CYLINDRIQUE )
! --- CONVERSION DES ANGLES NAUTIQUES EN ANGLES D'EULER
! =====================================================================
!
        call getElemOrientation(ndim, nno, igeom, angnau)
! ----- Select objects to construct from option name
        call behaviourOption(option, compor_copy, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)

! ----- Output fields
        if (lMatr) then
            call jevech('PMATUNS', 'E', imatuu)
        end if
        if (lVect) then
            call jevech('PVECTUR', 'E', ivectu)
        end if
        if (lVari) then
            call jevech('PVARIPR', 'E', ivarip)
        end if
        if (lSigm) then
            call jevech('PCONTPR', 'E', icontp)
            call tecach('OOO', 'PCONTMR', 'L', iret=iret, nval=2, itab=itabin)
            do iSigm = 1, itabin(2)
                zr(icontp-1+iSigm) = zr(icontm-1+iSigm)
            end do
            call jevech('PCODRET', 'E', jcret)
        end if
! ----- Compute
        codret = 0
        dimmat = nddls*nnop
        if (option(1:9) .eq. 'RIGI_MECA') then
            call xasshm(ds_thm, &
                        nno, npg, npi, ipoids, ivf, &
                        idfde, igeom, zr(igeom), zr(icarcr), zr(ideplm), &
                        zr(ideplm), zr(icontm), zr(icontp), zr(ivarim), zr(ivarim), &
                        defgem, defgep, drds, drdsr, dsde, &
                        b, dfdi, dfdi2, r, sigbar, &
                        c, ck, cs, zr(imatuu), zr(ivectu), &
                        zr(iinstm), zr(iinstp), option, zi(imate), mecani, &
                        press1, press2, tempe, dimdef, dimcon, &
                        dimuel, nbvari, nddls, nddlm, nmec, &
                        np1, ndim, compor_copy, axi, modint, &
                        codret, nnop, nnops, nnopm, enrmec, &
                        dimenr, zi(jheavt), zi(jlonch), zi(jcnset), jpintt, &
                        jpmilt, jheavn, angnau, dimmat, enrhyd, nfiss, nfh, jfisno, &
                        work1, work2, lVect, lMatr, lVari, lSigm)
        else
            do li = 1, dimuel
                zr(ideplp+li-1) = zr(ideplm+li-1)+zr(ideplp+li-1)
            end do
            call xasshm(ds_thm, &
                        nno, npg, npi, ipoids, ivf, &
                        idfde, igeom, zr(igeom), zr(icarcr), zr(ideplm), &
                        zr(ideplp), zr(icontm), zr(icontp), zr(ivarim), zr(ivarip), &
                        defgem, defgep, drds, drdsr, dsde, &
                        b, dfdi, dfdi2, r, sigbar, &
                        c, ck, cs, zr(imatuu), zr(ivectu), &
                        zr(iinstm), zr(iinstp), option, zi(imate), mecani, &
                        press1, press2, tempe, dimdef, dimcon, &
                        dimuel, nbvari, nddls, nddlm, nmec, &
                        np1, ndim, compor_copy, axi, modint, &
                        codret, nnop, nnops, nnopm, enrmec, &
                        dimenr, zi(jheavt), zi(jlonch), zi(jcnset), jpintt, &
                        jpmilt, jheavn, angnau, dimmat, enrhyd, nfiss, nfh, jfisno, &
                        work1, work2, lVect, lMatr, lVari, lSigm)
        end if
        if (lSigm) then
            zi(jcret) = codret
        end if
! =====================================================================
! --- SUPRESSION DES DDLS HEAVISIDE SUPERFLUS -------------------------
! =====================================================================
        call xhmddl(ndim, nfh, nddls, dimuel, nnop, nnops, &
                    zi(jstno), .false._1, option, nomte, zr(imatuu), &
                    zr(ivectu), nddlm, nfiss, jfisno, .false._1, contac)
    end if
! =====================================================================
! --- 3. OPTION : CHAR_MECA_PESA_R ------------------------------------
! =====================================================================
    if (option .eq. 'CHAR_MECA_PESA_R') then
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PMATERC', 'L', imate)
        call jevech('PPESANR', 'L', ipesa)
        call jevech('PVECTUR', 'E', ivectu)
        call rccoma(zi(imate), 'THM_DIFFU', 1, phenom, icodre(1))
        call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                    ' ', phenom, 0, ' ', [0.d0], &
                    1, 'RHO', rho(1), icodre, 1)
!
!        INDICATEUR POUR SAVOIR SI ON A DE L'ENRICHISSEMENT
        yaenrm = enrmec(1)
!
        call xpeshm(nno, nnop, nnops, ndim, nddls, &
                    nddlm, npg, igeom, jpintt, jpmilt, jheavn, &
                    ivf, ipoids, idfde, ivectu, ipesa, &
                    zi(jheavt), zi(jlonch), zi(jcnset), rho(1), axi, &
                    yaenrm, nfiss, nfh, jfisno)
!
! =====================================================================
! --- SUPRESSION DES DDLS HEAVISIDE SUPERFLUS -------------------------
! =====================================================================
        call xhmddl(ndim, nfh, nddls, dimuel, nnop, nnops, &
                    zi(jstno), .false._1, option, nomte, rbid, &
                    zr(ivectu), nddlm, nfiss, jfisno, .false._1, contac)
    end if
! ======================================================================
! --- 4. OPTION : FORC_NODA --------------------------------------------
! ======================================================================
    if (option .eq. 'FORC_NODA') then
! ======================================================================
! --- PARAMETRES EN ENTREE ---------------------------------------------
! ======================================================================
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PSIEFR', 'L', icontm)
        call jevech('PMATERC', 'L', imate)

! ----- Not a transient computation (CALC_CHAMP only !)
        dt = 0.d0
        fnoevo = .false.

! ======================================================================
! --- PARAMETRES EN SORTIE ---------------------------------------------
! ======================================================================
        call jevech('PVECTUR', 'E', ivectu)
!
        call xfnohm(ds_thm, &
                    fnoevo, dt, nno, npg, ipoids, &
                    ivf, idfde, zr(igeom), zr(icontm), b, &
                    dfdi, dfdi2, r, zr(ivectu), zi(imate), &
                    mecani, press1, dimcon, nddls, nddlm, &
                    dimuel, nmec, np1, ndim, axi, &
                    dimenr, nnop, nnops, nnopm, igeom, &
                    jpintt, jpmilt, jheavn, zi(jlonch), zi(jcnset), zi(jheavt), &
                    enrmec, enrhyd, nfiss, nfh, jfisno)
!
! =====================================================================
! --- SUPRESSION DES DDLS HEAVISIDE SUPERFLUS -------------------------
! =====================================================================
        call xhmddl(ndim, nfh, nddls, dimuel, nnop, nnops, &
                    zi(jstno), .false._1, option, nomte, rbid, &
                    zr(ivectu), nddlm, nfiss, jfisno, .false._1, contac)
    end if
!
    deallocate (work1)
    deallocate (work2)
    deallocate (b)
end subroutine
