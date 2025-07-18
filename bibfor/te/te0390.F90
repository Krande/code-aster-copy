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
subroutine te0390(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gddyng.h"
#include "asterfort/gdecva.h"
#include "asterfort/gdfine.h"
#include "asterfort/gdfint.h"
#include "asterfort/gdjrg0.h"
#include "asterfort/gdliva.h"
#include "asterfort/gdmine.h"
#include "asterfort/gdmrig.h"
#include "asterfort/gdpetk.h"
#include "asterfort/gdsig.h"
#include "asterfort/gdstag.h"
#include "asterfort/jevech.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: POU_D_T_GD
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
    real(kind=8) :: nu, instam, instap
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: elrefe
    real(kind=8) :: en(3, 2), enprim(3, 2), valres(4), granc(6), grani(4)
    real(kind=8) :: rigi(18, 18), fint(6, 3), y0(3), x00(3, 3), x0k(3, 3)
    real(kind=8) :: x0pg(3), tetak(3, 3), tetag(3), tetapg(3), qim(3, 3)
    real(kind=8) :: qikm1(3, 3), qik(3, 3), rot0(3, 3), rotm(3, 3), rotkm1(3, 3)
    real(kind=8) :: rotk(3, 3), petik(3), petikm(3), gn(3), gm(3), pn(3), pm(3)
    real(kind=8) :: x0sk(3, 3), rmkm1(3, 3), rmk(3, 3), omkm1(3, 3)
    real(kind=8) :: ompkm1(3, 3), omk(3, 3), ompk(3, 3), x0sec(3), rgmkm(3)
    real(kind=8) :: rgmk(3), omgkm(3), ompgkm(3), omgk(3), ompgk(3)
    character(len=16) :: rela_comp, defo_comp
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8), parameter :: nbres = 3
    integer(kind=8) :: icodre(nbres)
    character(len=16), parameter :: nomres(nbres) = (/'E  ', 'NU ', 'RHO'/)
    integer(kind=8) :: i, iacckm, iaccp, ico, iddepl, idepde
    integer(kind=8) :: idepkm, idepm, idfdk, ifint, igeom, imat, imate
    integer(kind=8) :: imatuu, instmr, instpr, ipoids, iret, iromk, iromkm
    integer(kind=8) :: istady, ivarim, ivarip, ivf, ivitkm, ivitp, j
    integer(kind=8) :: jcret, jefint, k0, k1, k2, k3
    integer(kind=8) :: k4, k5, k6, k7, kc, kp, ks
    integer(kind=8) :: iorien, lsig, lsigma, ne, nno
    integer(kind=8) :: nord, npg, codret
    real(kind=8) :: a, ajacob, alfnmk, ay, az, delnmk, demi
    real(kind=8) :: deux, e, g, pas, pjacob, r8bid, rho
    real(kind=8) :: stoudy, un, xiy, xiz, xjx
    integer(kind=8), parameter :: nb_cara = 6
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8), parameter :: noms_cara(nb_cara) = (/'A1   ', 'IY1  ', 'IZ1  ', &
                                                          'AY1  ', 'AZ1  ', 'JX1  '/)
!
! --------------------------------------------------------------------------------------------------
!
    call elref1(elrefe)
    r8bid = 0.d0
    demi = 5.d-1
    un = 1.d0
    deux = 2.d0
    rigi = 0.d0
    fint = 0.d0
!
!* STOUDY VAUT: 1., SI L'ON EST EN DYNAMIQUE
!*              0., SI L'ON EST EN STATIQUE
!
    call tecach('NNO', 'PSTADYN', 'L', iret, iad=istady)
    if (istady .eq. 0) then
        stoudy = 0.d0
    else
        stoudy = zr(istady)
    end if
!
    call elrefe_info(fami='RIGI', nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk)
    nord = 6*nno
!
    ico = 0
    do kp = 1, npg
        do ne = 1, nno
            ico = ico+1
            en(ne, kp) = zr(ivf-1+ico)
            enprim(ne, kp) = zr(idfdk-1+ico)
        end do
    end do
!
! - Get input fields
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCAORIE', 'L', iorien)
    call jevech('PVARIMP', 'L', ivarim)
    call jevech('PINSTMR', 'L', instmr)
    call jevech('PINSTPR', 'L', instpr)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLMR', 'L', idepm)

! ---- LA PRESENCE DU CHAMP DE DEPLACEMENT A L INSTANT T+
! ---- DEVRAIT ETRE CONDITIONNE  PAR L OPTION (AVEC RIGI_MECA_TANG
! ---- CA N A PAS DE SENS).
! ---- CEPENDANT CE CHAMP EST INITIALISE A 0 PAR LA ROUTINE NMMATR.
    call jevech('PDEPLPR', 'L', idepde)
    call jevech('PDDEPLA', 'L', iddepl)
!
!
! - Properties of behaviour
!
    rela_comp = compor(RELA_NAME)
    defo_comp = compor(DEFO)
!
! - Some checks
!
    if (rela_comp(1:4) .ne. 'ELAS') then
        call utmess('F', 'POUTRE0_17')
    end if
    if (defo_comp .ne. 'GROT_GDEP') then
        call utmess('F', 'POUTRE0_41')
    end if
!
! - Elastic properties
!
    call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                ' ', 'ELAS', 0, '  ', [r8bid], &
                2, nomres, valres, icodre, 1)
    call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                ' ', 'ELAS', 0, '  ', [r8bid], &
                1, nomres(3), valres(3), icodre(3), 0)
    if (icodre(3) .ne. 0) then
        if (stoudy .gt. demi) then
            call utmess('F', 'POUTRE0_42')
        else
            valres(3) = 0.d0
        end if
    end if
    e = valres(1)
    nu = valres(2)
    rho = valres(3)
    g = e/(deux*(un+nu))
!
!     --- RECUPERATION DES CARACTERISTIQUES GENERALES DES SECTIONS ---
!     --- LA SECTION EST SUPPOSEE CONSTANTE ---
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
!
    a = vale_cara(1)
    xiy = vale_cara(2)
    xiz = vale_cara(3)
    ay = vale_cara(4)
    az = vale_cara(5)
    xjx = vale_cara(6)
    granc(1) = e*a
    granc(2) = g*a/ay
    granc(3) = g*a/az
    granc(4) = g*xjx
    granc(5) = e*xiy
    granc(6) = e*xiz
!
!     --- RECUPERATION DES ORIENTATIONS INITIALES Y0(1), Y0(2), Y0(3)

    y0(1) = zr(iorien)
    y0(2) = zr(iorien+1)
    y0(3) = zr(iorien+2)
!
! - Select objects to construct from option name
!
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)
!
! - Get output fields
!
    if (lMatr) then
        call jevech('PMATUNS', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', jefint)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', lsigma)
        call jevech('PCODRET', 'E', jcret)
    end if
!
    k0 = igeom-1
    k1 = idepm-1
    k2 = idepde-1
    k3 = iddepl-1
!
    do ne = 1, nno
        do kc = 1, 3
            k0 = k0+1
            k1 = k1+1
            k2 = k2+1
            k3 = k3+1
            x00(kc, ne) = zr(k0)
            if (option .eq. 'RIGI_MECA_TANG') then
                x0k(kc, ne) = zr(k0)+zr(k1)
            else
                x0k(kc, ne) = zr(k0)+zr(k1)+zr(k2)
            end if
        end do
        do kc = 1, 3
            k1 = k1+1
            k2 = k2+1
            k3 = k3+1
            qim(kc, ne) = zr(k1)
            if (option .eq. 'RIGI_MECA_TANG') then
                qik(kc, ne) = 0.d0
                tetak(kc, ne) = 0.d0
            else
                qik(kc, ne) = zr(k2)
                tetak(kc, ne) = zr(k3)
            end if
        end do
    end do
!
    if (stoudy .gt. demi) then
!* ON TRAITE UN PROBLEME DYNAMIQUE
        instam = zr(instmr)
        instap = zr(instpr)
        pas = instap-instam
        grani(1) = rho*xjx
        grani(2) = rho*xiy
        grani(3) = rho*xiz
        grani(4) = rho*a
!* PARAMETRES ALPHA ET DELTA DE NEWMARK
        alfnmk = zr(istady+1)
        delnmk = zr(istady+2)
!
        call jevech('PDEPKM1', 'L', idepkm)
        call jevech('PVITKM1', 'L', ivitkm)
        call jevech('PACCKM1', 'L', iacckm)
        call jevech('PVITPLU', 'L', ivitp)
        call jevech('PACCPLU', 'L', iaccp)
        call jevech('PROMKM1', 'L', iromkm)
        call jevech('PROMK', 'L', iromk)
        k1 = idepkm-1
        k2 = ivitkm-1
        k3 = iacckm-1
        k4 = ivitp-1
        k5 = iaccp-1
        k6 = iromkm-1
        k7 = iromk-1
        do ne = 1, nno
            do kc = 1, 3
                k1 = k1+1
                k2 = k2+1
                k3 = k3+1
                k4 = k4+1
                k5 = k5+1
                k6 = k6+1
                k7 = k7+1
                x0sk(kc, ne) = zr(k5)
            end do
            do kc = 1, 3
                k1 = k1+1
                k2 = k2+1
                k3 = k3+1
                k4 = k4+1
                k5 = k5+1
                k6 = k6+1
                k7 = k7+1
                qikm1(kc, ne) = zr(k1)
                omkm1(kc, ne) = zr(k2)
                ompkm1(kc, ne) = zr(k3)
                omk(kc, ne) = zr(k4)
                ompk(kc, ne) = zr(k5)
                rmkm1(kc, ne) = zr(k6)
                rmk(kc, ne) = zr(k7)
            end do
        end do
    end if
!
!* BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, npg
        call gdjrg0(kp, nno, enprim, x00, y0, &
                    ajacob, rot0)
        pjacob = zr(ipoids-1+kp)*ajacob
!*** LECTURE, DANS 'PVARIMR', DU VECTEUR-COURBURE A L'ITER. PRECEDENTE
        call gdliva(kp, zr(ivarim), petikm)
!
        call gdstag(stoudy, kp, nno, ajacob, en, &
                    enprim, x0k, tetak, qim, qikm1, &
                    qik, x0pg, tetag, tetapg, rotm, &
                    rotkm1, rotk)
        call gdpetk(tetag, tetapg, petikm, petik)
!
!*** ECRITURE, DANS 'PVARIPR', DU VECTEUR-COURBURE ACTUALISE, CAR
!*** MAJSVT, UTILE POUR D'AUTRES ELEMENTS, COPIE VARIP DANS VARIM
        if (option .ne. 'RIGI_MECA_TANG') then
            call jevech('PVARIPR', 'E', ivarip)
            call gdecva(kp, petik, zr(ivarip))
        end if
!
        call gdsig('RIGI', kp, 1, x0pg, petik, &
                   rot0, rotk, granc, zi(imate), gn, &
                   gm, pn, pm)
        if (stoudy .gt. demi) then
!* ON TRAITE UN PROBLEME DYNAMIQUE
            call gddyng(kp, nno, en, x0sk, rmkm1, &
                        rmk, omkm1, ompkm1, omk, ompk, &
                        x0sec, rgmkm, rgmk, omgkm, ompgkm, &
                        omgk, ompgk)
        end if
!
        if (lMatr) then
            call gdmrig(kp, nno, ajacob, pjacob, en, &
                        enprim, x0pg, rot0, rotk, granc, &
                        pn, pm, rigi)
            if (stoudy .gt. demi) then
!* ON TRAITE UN PROBLEME DYNAMIQUE
                call gdmine(kp, nno, pjacob, en, grani, &
                            alfnmk, delnmk, pas, rot0, rotm, &
                            rotkm1, rotk, rmkm1, rmk, omgkm, &
                            ompgkm, omgk, ompgk, rigi)
            end if
        end if
!
        if (lVect) then
            ASSERT(lSigm)
            call gdfint(kp, nno, ajacob, pjacob, en, &
                        enprim, x0pg, pn, pm, fint)
            lsig = lsigma-1+(kp-1)*6
            do ks = 1, 3
                lsig = lsig+1
!*** ATTENTION : LE TORSEUR EST EXPRIME EN COORDONNEES LOCALES
                zr(lsig) = gn(ks)
            end do
            do ks = 1, 3
                lsig = lsig+1
                zr(lsig) = gm(ks)
            end do
            if (stoudy .gt. demi) then
                call gdfine(kp, nno, pjacob, en, grani, &
                            rot0, rotk, omgk, ompgk, fint)
            end if
        end if
!
!* FIN DE BOUCLE SUR LES POINTS DE GAUSS
!
    end do
!
    if (lMatr) then
        imat = imatuu-1
        do i = 1, nord
            do j = 1, nord
                imat = imat+1
                zr(imat) = rigi(i, j)
            end do
        end do
    end if
!
    if (lVect) then
        ifint = jefint-1
        do ne = 1, nno
            do kc = 1, 6
                ifint = ifint+1
                zr(ifint) = fint(kc, ne)
            end do
        end do
    end if
!
    if (lSigm) then
        zi(jcret) = codret
    end if
end subroutine
