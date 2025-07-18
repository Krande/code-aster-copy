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
! aslint: disable=W0413
!
subroutine te0239(option, nomte)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/defgen.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/effi.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/matdtd.h"
#include "asterfort/mattge.h"
#include "asterfort/moytpg.h"
#include "asterfort/nmcomp.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_AXIS
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ndimLdc = 2
    character(len=4), parameter :: fami = "RIGI"
    integer(kind=8) :: nbcou, npge, icontm, ideplm, ivectu, icou, inte, icontp
    integer(kind=8) :: kpki, k1, k2, kompt, ivarim, ivarip, iinstm, iinstp, lgpg, ideplp
    integer(kind=8) :: icarcr, nbvari, jcret, codret
    real(kind=8) :: cisail, zic, coef, rhos, rhot, epsx3, gsx3, sgmsx3
    real(kind=8) :: zmin, hic, depsx3
    integer(kind=8) :: itab(8), jnbspi
    character(len=8) :: nompar, elrefe
    real(kind=8) :: tempm
    real(kind=8) :: dfdx(3), zero, un, deux
    real(kind=8) :: test, test2, eps, nu, h, cosa, sina, cour, r
    real(kind=8) :: jacp, kappa, correc
    real(kind=8) :: eps2d(4), deps2d(4), sigtdi(5), sigmtd(5)
    real(kind=8) :: x3
    real(kind=8) :: dtild(5, 5), dtildi(5, 5), dsidep(6, 6)
    real(kind=8) :: rtangi(9, 9), rtange(9, 9), sigm2d(4), sigp2d(4)
    real(kind=8) :: angmas(3)
    integer(kind=8) :: nno, kp, npg, i, j, k, imatuu, icaco, ndimv
    integer(kind=8) :: ivarix
    integer(kind=8) :: ipoids, ivf, idfdk, igeom, imate
    integer(kind=8) :: nbpar, cod, iret, ksp
    aster_logical :: testl1, testl2
    type(Behaviour_Integ) :: BEHinteg
    integer(kind=8), parameter :: nbres = 2
    character(len=16), pointer :: compor(:) => null()
    character(len=16), parameter :: nomres(nbres) = (/'E ', 'NU'/)
    integer(kind=8) :: valret(nbres)
    real(kind=8) :: valres(nbres)
    parameter(npge=3)
    data zero, un, deux/0.d0, 1.d0, 2.d0/
    character(len=16) :: defo_comp, rela_comp, rela_cpla
    aster_logical :: lVect, lMatr, lVari, lSigm
    character(len=8), parameter :: typmod(2) = (/'C_PLAN  ', '        '/)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    ivarip = 1
!
    eps = 1.d-3
    codret = 0

!   Angle du mot clef MASSIF de AFFE_CARA_ELEM, initialisé à r8nnem (on ne s'en sert pas)
    angmas = r8nnem()

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

!
    call elref1(elrefe)
    call elrefe_info(fami='RIGI', nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdk)

! - Get input fields
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCACOQU', 'L', icaco)
    call jevech('PMATERC', 'L', imate)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', icarcr)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=itab)
!      LGPG = MAX(ITAB(6),1)*ITAB(7) resultats faux sur Bull avec ifort
    if (itab(6) .le. 1) then
        lgpg = itab(7)
    else
        lgpg = itab(6)*itab(7)
    end if
    call jevech('PNBSP_I', 'L', jnbspi)

! - Properties of shell
    h = zr(icaco)
    kappa = zr(icaco+1)
    correc = zr(icaco+2)
    zmin = -h/2.d0

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndimLdc, typmod, option, &
                              compor, zr(icarcr), &
                              zr(iinstm), zr(iinstp), &
                              fami, zi(imate), &
                              BEHinteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Properties of behaviour
    rela_comp = compor(RELA_NAME)
    defo_comp = compor(DEFO)
    rela_cpla = compor(PLANESTRESS)
    read (compor(NVAR), '(I16)') nbvari

! - Some checks
    if (rela_cpla .eq. 'COMP_ELAS') then
        if (rela_comp .ne. 'ELAS') then
            call utmess('F', 'PLATE1_8')
        end if
    end if
    if (defo_comp(6:10) .eq. '_REAC') then
        call utmess('A', 'PLATE1_9', sk=defo_comp)
    end if

!--- NBRE DE  COUCHES ET LONG. MAX
    nbcou = zi(jnbspi-1+1)
    if (nbcou .le. 0) then
        call utmess('F', 'PLATE1_10')
    end if
    if (nbcou .gt. 10) then
        call utmess('F', 'PLATE1_11')
    end if
!---- EPAISSEUR DE CHAQUE COUCHE
    hic = h/nbcou
    ndimv = npg*npge*nbcou*nbvari
!
! - Get output fields
!
    if (lMatr) then
        call jevech('PMATUUR', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PVARIMP', 'L', ivarix)
        b_n = to_blas_int(ndimv)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
    end if
!
    call r8inir(81, 0.d0, rtange, 1)
    kpki = 0
!-- DEBUT DE BOUCLE D'INTEGRATION SUR LA SURFACE NEUTRE
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm1d(nno, zr(ipoids+kp-1), zr(idfdk+k), zr(igeom), dfdx, &
                    cour, jacp, cosa, sina)
        r = zero
!
        call r8inir(5, 0.d0, sigmtd, 1)
        call r8inir(25, 0.d0, dtild, 1)
!
!-- BOUCLE SUR LES POINTS D'INTEGRATION SUR LA SURFACE
!
        do i = 1, nno
            r = r+zr(igeom+2*i-2)*zr(ivf+k+i-1)
        end do
!
!===============================================================
!     -- RECUPERATION DE LA TEMPERATURE POUR LE MATERIAU:
!     -- SI LA TEMPERATURE EST CONNUE AUX NOEUDS :
        call moytpg('RIGI', kp, 3, '-', tempm, &
                    iret)
        nbpar = 1
        nompar = 'TEMP'
        call rcvalb('RIGI', kp, 1, '-', zi(imate), &
                    ' ', 'ELAS', nbpar, nompar, [tempm], &
                    2, nomres, valres, valret, 1)
!
        nu = valres(2)
        cisail = valres(1)/(un+nu)
!
!       ON EST EN AXIS:
        jacp = jacp*r
!
        test = abs(h*cour/deux)
        if (test .ge. un) correc = zero
        test2 = abs(h*cosa/(deux*r))
        if (test2 .ge. un) correc = zero
!
        testl1 = (test .le. eps .or. correc .eq. zero)
        testl2 = ( &
                 test2 .le. eps .or. correc .eq. zero .or. abs(cosa) .le. eps .or. abs(cour*r) &
                 .le. eps .or. abs(cosa-cour*r) .le. eps &
                 )
!
!-- DEBUT DE BOUCLE D'INTEGRATION DANS L'EPAISSEUR
!
        do icou = 1, nbcou
            do inte = 1, npge
                if (inte .eq. 1) then
                    zic = zmin+(icou-1)*hic
                    coef = 1.d0/3.d0
                else if (inte .eq. 2) then
                    zic = zmin+hic/2.d0+(icou-1)*hic
                    coef = 4.d0/3.d0
                else
                    zic = zmin+hic+(icou-1)*hic
                    coef = 1.d0/3.d0
                end if
!
                x3 = zic
!
                if (testl1) then
                    rhos = 1.d0
                else
                    rhos = 1.d0+x3*cour
                end if
                if (testl2) then
                    rhot = 1.d0
                else
                    rhot = 1.d0+x3*cosa/r
                end if
!
!           CALCULS DES COMPOSANTES DE DEFORMATIONS TRIDIMENSIONNELLES :
!           EPSSS, EPSTT, EPSSX3 (EN FONCTION DES DEFORMATIONS
!           GENERALISEES :ESS,KSS,ETT,KTT,GS)
!           DE L'INSTANT PRECEDANT ET DES DEFORMATIONS INCREMENTALES
!           DE L'INSTANT PRESENT
!
                call defgen(testl1, testl2, nno, r, x3, &
                            sina, cosa, cour, zr(ivf+k), dfdx, &
                            zr(ideplm), eps2d, epsx3)
                call defgen(testl1, testl2, nno, r, x3, &
                            sina, cosa, cour, zr(ivf+k), dfdx, &
                            zr(ideplp), deps2d, depsx3)
!
!
!           CONSTRUCTION DE LA DEFORMATION GSX3
!           ET DE LA CONTRAINTE SGMSX3
                gsx3 = 2.d0*(epsx3+depsx3)
                sgmsx3 = cisail*kappa*gsx3/2.d0
!
!           CALCUL DU NUMERO DU POINT D'INTEGRATION COURANT
                kpki = kpki+1
                k1 = 4*(kpki-1)
                k2 = lgpg*(kp-1)+(npge*(icou-1)+inte-1)*nbvari
                ksp = (icou-1)*npge+inte
!
                do i = 1, 4
                    sigm2d(i) = zr(icontm+k1+i-1)
                end do

! ------------- Set main parameters for behaviour (on point)
                call behaviourSetParaPoin(kp, ksp, BEHinteg)

! ------------- Integrate
!     INTEGRATION DE LA LOI DE COMPORTEMENT POUR LES COQUE_1D :
!     COQUE_AXIS : COMPORTEMENT C_PLAN
!     COQUE_AXIS, DIRECTION Y : CPLAN (EVT DEBORST)
!     DIRECTION Z : EPSZZ CONNU
                dsidep = 0.d0
                sigp2d = 0.d0
                cod = 0
                call nmcomp(BEHinteg, &
                            fami, kp, ksp, ndimLdc, typmod, &
                            zi(imate), compor, zr(icarcr), zr(iinstm), zr(iinstp), &
                            4, eps2d, deps2d, 4, sigm2d, &
                            zr(ivarim+k2), option, angmas, &
                            sigp2d, zr(ivarip+k2), 36, dsidep, cod)

                if (lSigm) then
                    do i = 1, 4
                        zr(icontp+k1+i-1) = sigp2d(i)
                    end do
                end if
!
!           COD=1 : ECHEC INTEGRATION LOI DE COMPORTEMENT
!           COD=3 : C_PLAN DEBORST SIGZZ NON NUL
                if (cod .ne. 0) then
                    if (codret .ne. 1) then
                        codret = cod
                    end if
                end if
!
!
                if (lMatr) then
!-- CALCULS DE LA MATRICE TANGENTE : BOUCLE SUR L'EPAISSEUR
!-- CONSTRUCTION DE LA MATRICE DTD (DTILD)
                    call matdtd(nomte, testl1, testl2, dsidep, cisail, &
                                x3, cour, r, cosa, kappa, &
                                dtildi)
                    do i = 1, 5
                        do j = 1, 5
                            dtild(i, j) = dtild(i, j)+dtildi(i, j)*0.5d0*hic*coef
                        end do
                    end do
                end if
!
                if (lVect) then
!-- CALCULS DES EFFORTS INTERIEURS : BOUCLE SUR L'EPAISSEUR
!
                    sigtdi(1) = zr(icontp-1+k1+1)/rhos
                    sigtdi(2) = x3*zr(icontp-1+k1+1)/rhos
                    sigtdi(3) = zr(icontp-1+k1+2)/rhot
                    sigtdi(4) = x3*zr(icontp-1+k1+2)/rhot
                    sigtdi(5) = sgmsx3/rhos
!
                    do i = 1, 5
                        sigmtd(i) = sigmtd(i)+sigtdi(i)*0.5d0*hic*coef
                    end do
                end if
!-- FIN DE BOUCLE SUR LES POINTS D'INTEGRATION DANS L'EPAISSEUR
            end do
        end do
!
        if (lVect) then
!-- CALCUL DES EFFORTS INTERIEURS
            call effi(nomte, sigmtd, zr(ivf+k), dfdx, jacp, &
                      sina, cosa, r, zr(ivectu))
        end if
        if (lMatr) then
!-- CONSTRUCTION DE LA MATRICE TANGENTE
            call mattge(nomte, dtild, sina, cosa, r, &
                        jacp, zr(ivf+k), dfdx, rtangi)
            do i = 1, 9
                do j = 1, 9
                    rtange(i, j) = rtange(i, j)+rtangi(i, j)
                end do
            end do
        end if
    end do
    if (lMatr) then
!-- STOCKAGE DE LA MATRICE TANGENTE
        kompt = 0
        do j = 1, 9
            do i = 1, j
                kompt = kompt+1
                zr(imatuu-1+kompt) = rtange(i, j)
            end do
        end do
    end if
!
    if (lSigm) then
        zi(jcret) = codret
    end if
end subroutine
