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
! aslint: disable=W1501
!
subroutine te0409(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/coqrep.h"
#include "asterfort/coqgth.h"
#include "asterfort/dhrc_lc.h"
#include "asterfort/dhrc_recup_mate.h"
#include "asterfort/dkqbf.h"
#include "asterfort/dktbf.h"
#include "asterfort/dsqbfb.h"
#include "asterfort/dstbfb.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxqbm.h"
#include "asterfort/dxqloc.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtbm.h"
#include "asterfort/dxtloc.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/glrc_recup_mate.h"
#include "asterfort/glrc_lc.h"
#include "asterfort/glrcmm.h"
#include "asterfort/gquad4.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/maglrc.h"
#include "asterfort/kit_glrc_dm_vmis.h"
#include "asterfort/pmrvec.h"
#include "asterfort/q4gbc.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/t3gbc.h"
#include "asterfort/tecach.h"
#include "asterfort/utbtab.h"
#include "asterfort/utctab.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: Q4GG / DKTG
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

!
! TOUTES LES VARIABLES EN DEHORS
! DE LA BOUCLE SUR LES POINTS DE GAUSS PEUVENT ETRE
! INITIALISES AU MOMENT DE LA DECLARATION
! matloc et vecloc SONT INITIALISES EN FONCTION DE
! L OPTION FULL_MECA, RIGI_MECA_*
    integer(kind=8) :: codret
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: ul(6, 4), dul(6, 4), angmas(3)
    real(kind=8) :: matloc((6*4)*(6*4+1)/2), vecloc(6, 4)
!  ----------------
!  CMPS DE DEPLACEMENT :
!   - MEMBRANE : DX(N1), DY(N1), DX(N2), ..., DY(NNO)
!   - FLEXION  : DZ(N1), BETAX(N1), BETAY(N1), DZ(N2), ..., BETAY(NNO)
!  CMPS DE DEFORMATION COQUE :
!   - MEMBRANE : EPSIXX,EPSIYY,2*EPSIXY
!   - FLEXION  : KHIXX,KHIYY,2*KHIXY
!  CMPS D' EFFORTS COQUE :
!   - MEMBRANE : NXX,NYY,NXY
!   - FLEXION  : MXX,MYY,MXY
!   - CISAILLEMENT : QX , QY
! --------------------------------------------------------------------
!            NPG:    NOMBRE DE POINTS DE GAUSS PAR ELEMENT
!            NC :    NOMBRE DE COTES DE L'ELEMENT
!
    real(kind=8) :: poids
!            POIDS:  POIDS DE GAUSS (Y COMPRIS LE JACOBIEN)
!            AIRE:   SURFACE DE L'ELEMENT
    real(kind=8) :: um(2, 4), uf(3, 4), dum(2, 4), duf(3, 4)
!            UM:     DEPLACEMENT (MEMBRANE) "-"
!            UF:     DEPLACEMENT (FLEXION)  "-"
!           DUM:     INCREMENT DEPLACEMENT (MEMBRANE)
!           DUF:     INCREMENT DEPLACEMENT (FLEXION)
    real(kind=8) :: eps(3), khi(3), gam(2), deps(6), dkhi(3), dgam(2), n(3)
    real(kind=8) :: m(3), q(2)
!            EPS:    DEFORMATION DE MEMBRANE "-"
!            DEPS:   INCREMENT DE DEFORMATION DE MEMBRANE
!            KHI:    DEFORMATION DE FLEXION  "-" (COURBURE)
!            DKHI:   INCREMENT DE DEFORMATION DE FLEXION (COURBURE)
!            GAM :   DISTORSIONS TRANSVERSES
!            DGAM :  INCREMENT DES DISTORSIONS TRANSVERSES
!            N  :    EFFORT NORMAL "+"
!            M  :    MOMENT FLECHISSANT "+"
!            Q  :    CISAILLEMENT TRANSVERSE
!
    real(kind=8) :: sigmam(32), efform(32), effint(32), efforp(32)
!            SIGMAM : EFFORTS DANS LE REPERE UTILISATEUR  A T "-"
!            EFFORM : EFFORTS DANS LE REPERE DE L'ELEMENT A T "-"
!            EFFINT : EFFORTS DANS LE REPERE DE L'ELEMENT A T "+"
!            EFFORP : EFFORTS DANS LE REPERE UTILISATEUR  A T "+"
!
    real(kind=8) :: df(9), dm(9), dmf(9), dc(2, 2)
!            DF :    MATRICE DE RIGIDITE TANGENTE MATERIELLE  (FLEXION)
!            DM :    MATRICE DE RIGIDITE TANGENTE MATERIELLE  (MEMBRANE)
!            DMF:    MATRICE DE RIGIDITE TANGENTE MATERIELLE  (COUPLAGE)
!            DC:     MATRICE DE RIGIDITE ELASTIQUE MATERIELLE
!                                                         (CISAILLEMENT)
!
    real(kind=8) :: dff(3, 3), dmm(3, 3), dmff(3, 3)
    real(kind=8) :: dcc(2, 2), dci(2, 2), dmc(3, 2), dfc(3, 2)
!
    real(kind=8) :: bf(3, 3*4), bm(3, 2*4), bmq(2, 3), bc(2, 3*4)
!            BF :    MATRICE "B" (FLEXION)
!            BM :    MATRICE "B" (MEMBRANE)
!            BC :    MATRICE "B" (CISAILLEMENT)
    real(kind=8) :: flex(3*4, 3*4), memb(2*4, 2*4), flexi(3*4, 3*4)
    real(kind=8) :: mefl(2*4, 3*4), work(3, 3*4)
!           MEMB:    MATRICE DE RIGIDITE DE MEMBRANE
!           FLEX:    MATRICE DE RIGIDITE DE FLEXION
!           WORK:    TABLEAU DE TRAVAIL
!           MEFL:    MATRICE DE COUPLAGE MEMBRANE-FLEXION
!           LE MATERIAU EST SUPPOSE HOMOGENE
!
    real(kind=8) :: t2iu(4), t2ui(4), t1ve(9), c, s
!
    aster_logical :: t3g, q4g
    aster_logical :: leul, lLinear
    aster_logical :: lVect, lMatr, lVari, lSigm
    aster_logical :: q4gg
    aster_logical :: coupmf
    aster_logical :: is_param_opt(2)
!
    integer(kind=8) :: ndim, nno, npg, ipoids, icoopg
    integer(kind=8) :: imate, iret, icontm, ivarim, igeom, icarcr, ideplm, ideplp
    integer(kind=8) :: icacoq, icontp, ivarip, ino, nbcont, ivectu, jcret, imatuu
    integer(kind=8) :: nbvari, ipg
    integer(kind=8) :: i, i1, i2, j, k, l
    integer(kind=8) :: icpg, icpv
    integer(kind=8) :: jtab(7), nbsig
    integer(kind=8) :: multic
!
    real(kind=8) :: delas(6, 6), dsidep(6, 6)
    real(kind=8) :: lambda, deuxmu, deumuf, lamf, gt, gc, gf, seuil, alphaf
    real(kind=8) :: alpha, beta
    real(kind=8) :: excen

!   -- variables pour dhrc
    real(kind=8) :: a0(6, 6), c0(2, 2, 2)
    real(kind=8) :: aa_t(6, 6, 2), ab_(6, 2, 2), ac_(2, 2, 2), aa_c(6, 6, 2)
    real(kind=8) :: ga_t(6, 6, 2), gb_(6, 2, 2), gc_(2, 2, 2), ga_c(6, 6, 2)
    real(kind=8) :: cstseu(6)
!   -- attention la taille de ecp depend du nombre de variable interne
!   -- lors de l ajout de variable interne il faut incrementer ecr et ecrp
    integer(kind=8), parameter :: nbvarmax = 36
    real(kind=8) :: epst(6), ep, surfgp, sig(8), dsig(8), ecr(nbvarmax), ecrp(nbvarmax)
    real(kind=8) :: epsm(6), qsi, eta, ctor
    real(kind=8) :: carat3(21), jacob(5), caraq4(25)
    real(kind=8) :: matr(50), sigm(8), alfmc
    real(kind=8) :: epsi_c, epsi_els, epsi_lim, val_param_opt(10)
!
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: type_comp, mult_comp, rela_plas, rela_comp, defo_comp
!
    integer(kind=8)      ::  codret2(1)
    character(len=32) :: elasKeyword
!
! --------------------------------------------------------------------------------------------------
!
    pgl = 0.0
    xyzl = 0.0d0
    ul = 0.d0
    dul = 0.d0
    angmas = 0.0d0
    dff = 0.0d0
    dmm = 0.0d0
    dmff = 0.0d0
    dcc = 0.0d0
    dci = 0.0d0
    dmc = 0.0d0
    dfc = 0.0d0
    bf = 0.d0
    bm = 0.d0
    bmq = 0.d0
    bc = 0.d0
    flex = 0.d0
    memb = 0.d0
    flexi = 0.d0
    mefl = 0.d0
    work = 0.d0
    t2iu = 0.0d0
    t2ui = 0.0d0
    multic = -1
    codret = 0
    nbsig = 6
    q4gg = ASTER_FALSE
    t3g = ASTER_FALSE
    q4g = ASTER_FALSE
    leul = ASTER_FALSE
    coupmf = ASTER_FALSE
    mult_comp = ' '
!
    if (nomte .eq. 'MEDKTG3' .or. nomte .eq. 'MET3GG3') then
        t3g = ASTER_TRUE
    else if (nomte .eq. 'MEDKQG4' .or. nomte .eq. 'MEQ4GG4') then
        q4g = ASTER_TRUE
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (nomte .eq. 'MEQ4GG4' .or. nomte .eq. 'MET3GG3') then
        q4gg = ASTER_TRUE
        nbsig = 8
    end if
!
    lLinear = option .eq. 'RIGI_MECA'
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg)
    call jevech('PGEOMER', 'L', igeom)
!
    if (nno .eq. 3) then
        call dxtpgl(zr(igeom), pgl)
    else if (nno .eq. 4) then
        call dxqpgl(zr(igeom), pgl)
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call utpvgl(nno, 3, pgl, zr(igeom), xyzl)
!
    call jevech('PMATERC', 'L', imate)
! ---   COQUE HOMOGENEISEE ?
    if (.not. lLinear) then
        call rccoma(zi(imate), 'ELAS', 1, elasKeyword, codret2(1))
        if ((elasKeyword .eq. 'ELAS_COQUE') .or. (elasKeyword .eq. 'ELAS_COQMU') .or. &
            (elasKeyword .eq. 'ELAS_ORTH')) then
            call utmess('F', 'PLATE1_5')
        end if
    end if

    if (lLinear) then
        lMatr = ASTER_TRUE
        lVect = ASTER_FALSE
        lVari = ASTER_FALSE
        lSigm = ASTER_FALSE
        codret = 0
    else
! ----- Select objects to construct from option name
        call jevech('PCOMPOR', 'L', vk16=compor)
        call behaviourOption(option, compor, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)
    end if

    if (.not. lLinear) then
        call jevech('PCARCRI', 'L', icarcr)
        call jevech('PVARIMR', 'L', ivarim)

        rela_comp = compor(RELA_NAME)
        type_comp = compor(INCRELAS)
        defo_comp = compor(DEFO)
        leul = defo_comp .eq. 'GROT_GDEP'
        read (compor(NVAR), '(I16)') nbvari
!           -- on verifie que le nombre de varint tient dans ecr
        ASSERT(nbvari .le. nbvarmax)
!
        call tecach('OOO', 'PCONTMR', 'L', iret, nval=7, itab=jtab)
        icontm = jtab(1)
        ASSERT(npg .eq. jtab(3))
!
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)

        if (defo_comp(6:10) .eq. '_REAC') then
            call utmess('F', 'PLATE1_6')
        end if

        if (leul) then
            do i = 1, nno
                i1 = 3*(i-1)
                i2 = 6*(i-1)
                zr(igeom+i1) = zr(igeom+i1)+zr(ideplm+i2)+zr(ideplp+i2)
                zr(igeom+i1+1) = zr(igeom+i1+1)+zr(ideplm+i2+1)+zr(ideplp+i2+1)
                zr(igeom+i1+2) = zr(igeom+i1+2)+zr(ideplm+i2+2)+zr(ideplp+i2+2)
            end do
!
            if (nno .eq. 3) then
                call dxtpgl(zr(igeom), pgl)
            else if (nno .eq. 4) then
                call dxqpgl(zr(igeom), pgl)
            end if
!
            call utpvgl(nno, 3, pgl, zr(igeom), xyzl)
        end if
!
        call utpvgl(nno, 6, pgl, zr(ideplm), ul)
        call utpvgl(nno, 6, pgl, zr(ideplp), dul)
    else
        nbvari = 0
        call tecach('NNO', 'PMATERC', 'L', iret, iad=imate)
        call rccoma(zi(imate), 'ELAS', 1, rela_comp(1:10), iret)
        icarcr = 1
        ivarim = 1
        icontm = 1
    end if

    call r8inir(8, 0.d0, sig, 1)
    call r8inir(8, 0.d0, dsig, 1)
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PCACOQU', 'L', icacoq)
!       --  EPAISSEUR TOTALE
    ep = zr(icacoq)
    ctor = zr(icacoq+3)
    excen = zr(icacoq+4)
!
    icontp = icontm
    ivarip = ivarim
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
!
!   -- GRANDEURS GEOMETRIQUES :
    if (nno .eq. 3) then
        call gtria3(xyzl, carat3)
    else if (nno .eq. 4) then
        call gquad4(xyzl, caraq4)
    end if
!
!   -- MISES A ZERO :
    if (lMatr) then
        call r8inir((3*nno)*(3*nno), 0.d0, flexi, 1)
        call r8inir((3*nno)*(3*nno), 0.d0, flex, 1)
        call r8inir((2*nno)*(2*nno), 0.d0, memb, 1)
        call r8inir((2*nno)*(3*nno), 0.d0, mefl, 1)
        call r8inir((6*4)*(6*4+1)/2, 0.d0, matloc, 1)
    end if
!
    if (lVect) then
        call r8inir(6*nno, 0.d0, vecloc, 1)
    end if
    call r8inir(32, 0.d0, effint, 1)
    call r8inir(32, 0.d0, efforp, 1)
!
    call r8inir(36, 0.d0, delas, 1)
    call r8inir(32, 0.d0, sigmam, 1)
    call r8inir(32, 0.d0, efform, 1)
    call r8inir(2*4, 0.d0, um, 1)
    call r8inir(3*4, 0.d0, uf, 1)
    call r8inir(2*4, 0.d0, dum, 1)
    call r8inir(3*4, 0.d0, duf, 1)
!
!   -- PARTITION DU DEPLACEMENT EN MEMBRANE/FLEXION :
    do ino = 1, nno
        um(1, ino) = ul(1, ino)
        um(2, ino) = ul(2, ino)
        uf(1, ino) = ul(3, ino)
        uf(2, ino) = ul(5, ino)
        uf(3, ino) = -ul(4, ino)
        dum(1, ino) = dul(1, ino)
        dum(2, ino) = dul(2, ino)
        duf(1, ino) = dul(3, ino)
        duf(2, ino) = dul(5, ino)
        duf(3, ino) = -dul(4, ino)
    end do
!
!       -- CONTRAINTE 2D : NXX,NYY,NXY,MXX,MYY,MXY,QX,QY
    nbcont = 8
!
    if (.not. lLinear) then
!      -- SIGMAM : EFFORTS DANS REPERE UTILISATEUR
        do i = 1, nbcont*npg
            sigmam(i) = zr(icontm-1+i)
        end do
!
!      -- CALCUL DES MATRICES DE CHANGEMENT DE REPERE
!              T2IU : MATRICE DE PASSAGE (2x2) ; UTILISATEUR -> INTRINSEQUE
!              T2UI : MATRICE DE PASSAGE (2x2) ; INTRINSEQUE -> UTILISATEUR
!
        alpha = zr(icacoq+1)*r8dgrd()
        beta = zr(icacoq+2)*r8dgrd()
        call coqrep(pgl, alpha, beta, t2iu, t2ui, c, s)
!
!       -- PASSAGE DES EFFORTS GENERALISES AUX POINTS D'INTEGRATION
!              DU REPERE UTILISATEUR AU REPERE INTRINSEQUE
        call dxefro(npg, t2ui, sigmam, efform)
    end if
!
!     BOUCLE SUR LES POINTS DE GAUSS DE LA SURFACE:
!
    do ipg = 1, npg
        call r8inir(nbvarmax, 0.d0, ecrp, 1)
        call r8inir(3, 0.d0, n, 1)
        call r8inir(3, 0.d0, m, 1)
        call r8inir(2, 0.d0, q, 1)
        call r8inir(9, 0.d0, df, 1)
        call r8inir(9, 0.d0, dm, 1)
        call r8inir(9, 0.d0, dmf, 1)
        call r8inir(4, 0.d0, dc, 1)
!
        qsi = zr(icoopg-1+ndim*(ipg-1)+1)
        eta = zr(icoopg-1+ndim*(ipg-1)+2)
!
        icpg = (ipg-1)*nbcont
        icpv = (ipg-1)*nbvari
!
        if (nomte .eq. 'MEDKTG3') then
            call dxtbm(carat3(9), bm)
            call dktbf(qsi, eta, carat3, bf)
            poids = zr(ipoids+ipg-1)*carat3(7)
        else if (nomte .eq. 'MEDKQG4') then
            call jquad4(xyzl, qsi, eta, jacob)
            call dxqbm(qsi, eta, jacob(2), bm)
            call dkqbf(qsi, eta, jacob(2), caraq4, bf)
            poids = zr(ipoids+ipg-1)*jacob(1)
        else if (nomte .eq. 'MEQ4GG4') then
            call jquad4(xyzl, qsi, eta, jacob)
            call dxqbm(qsi, eta, jacob(2), bm)
            call dsqbfb(qsi, eta, jacob(2), bf)
            call q4gbc(qsi, eta, jacob(2), caraq4, bc)
            poids = zr(ipoids+ipg-1)*jacob(1)
        else if (nomte .eq. 'MET3GG3') then
            call dxtbm(carat3(9), bm)
            call dstbfb(carat3(9), bf)
            call t3gbc(xyzl, qsi, eta, bc)
            poids = carat3(8)
        end if
!

        call r8inir(3, 0.d0, eps, 1)
        call r8inir(6, 0.d0, deps, 1)
        call r8inir(6, 0.d0, epsm, 1)
        call r8inir(3, 0.d0, khi, 1)
        call r8inir(3, 0.d0, dkhi, 1)
        call pmrvec('ZERO', 3, 2*nno, bm, um, eps)
        call pmrvec('ZERO', 3, 2*nno, bm, dum, deps)
        call pmrvec('ZERO', 3, 3*nno, bf, uf, khi)
        call pmrvec('ZERO', 3, 3*nno, bf, duf, dkhi)
!
        if (q4gg) then
            call r8inir(2, 0.d0, gam, 1)
            call r8inir(2, 0.d0, dgam, 1)
            call pmrvec('ZERO', 2, 3*nno, bc, uf, gam)
            call pmrvec('ZERO', 2, 3*nno, bc, duf, dgam)
        end if

!           -- EULER_ALMANSI - TERMES QUADRATIQUES
        if (leul) then
            call r8inir(6, 0.d0, bmq, 1)
!
            do i = 1, 2
                do k = 1, nno
                    do j = 1, 2
                        bmq(i, j) = bmq(i, j)+bm(i, 2*(k-1)+i)*dum(j, k)
                    end do
                    bmq(i, 3) = bmq(i, 3)+bm(i, 2*(k-1)+i)*duf(1, k)
                end do
            end do
!
            do k = 1, 3
                do i = 1, 2
                    deps(i) = deps(i)-0.5d0*bmq(i, k)*bmq(i, k)
                end do
                deps(3) = deps(3)-bmq(1, k)*bmq(2, k)
            end do
        end if
!
        call r8inir(50, 0.d0, matr, 1)
!
        if (.not. lLinear) then
            do i = 1, 3
                deps(i+3) = dkhi(i)
                epsm(i) = eps(i)
                epsm(i+3) = khi(i)
            end do

            do i = 1, nbsig
                sig(i) = efform(icpg+i)
                sigm(i) = sig(i)
            end do
        end if
!
!EXCENTREMENT RAJOUTE UN COUPLAGE MEMBRANE_FLEXION : EPSI = EPSI+EXCENT*KHI
        if (rela_comp(1:4) .ne. 'ELAS' .and. excen .gt. 0.0d0) then
            do i = 1, 3
                epsm(i) = epsm(i)+excen*khi(i)
                deps(i) = deps(i)+excen*khi(i)
            end do
        end if

        if (rela_comp(1:4) .eq. 'ELAS') then

            call r8inir(3*3, 0.d0, dff, 1)
            call r8inir(3*3, 0.d0, dmm, 1)
            call r8inir(3*3, 0.d0, dmff, 1)
            call r8inir(2*2, 0.d0, dcc, 1)
            call r8inir(2*2, 0.d0, dci, 1)
            call r8inir(3*2, 0.d0, dmc, 1)
            call r8inir(3*2, 0.d0, dfc, 1)
            call dxmate('RIGI', dff, dmm, dmff, dcc, &
                        dci, dmc, dfc, nno, pgl, &
                        multic, coupmf, t2iu, t2ui, t1ve)
            call r8inir(36, 0.d0, dsidep, 1)
            do i = 1, 3
                do j = 1, 3
!                       -- MEMBRANE
                    dsidep(i, j) = dmm(i, j)
!                       -- FLEXION
                    dsidep(i+3, j+3) = dff(i, j)
!                       -- MEMBRANE-FLEXION
                    dsidep(i, j+3) = dmff(i, j)
                    dsidep(j+3, i) = dmff(i, j)
                end do
            end do
!
!EXCENTREMENT RAJOUTE UN COUPLAGE MEMBRANE_FLEXION : EPSI = EPSI+EXCENT*KHI
! ON NE FAIT RIEN CAR LE COUPLAGE MEMBRANE_FLEXION INDUIT EST PRIS EN COMPTE DANS DXMATE
!               --  prise en compte de la dilatation thermique
            if (.not. lLinear) then
                call coqgth(zi(imate), rela_comp, 'RIGI', ipg, ep, epsm, deps)
            end if
!
!               -- calcul de l'accroissement de contrainte
            do i = 1, 6
                dsig(i) = 0.d0
                do j = 1, 6
                    dsig(i) = dsig(i)+dsidep(i, j)*deps(j)
                end do
            end do
!               -- calcul de l'accoissement effort cisaillement
            if (q4gg) then
                dsig(7) = dcc(1, 1)*dgam(1)
                dsig(8) = dcc(2, 2)*dgam(2)
            end if
!
            do i = 1, nbsig
                sig(i) = sig(i)+dsig(i)
            end do
!
        else if (rela_comp(1:11) .eq. 'GLRC_DAMAGE') then
            do i = 1, nbvari
                ecr(i) = zr(ivarim-1+icpv+i)
            end do
!
            call maglrc(zi(imate), matr, delas, ecr)
            if (q4gg) then
                dcc(1, 1) = matr(14)
                dcc(2, 2) = matr(15)
                dcc(1, 2) = 0.d0
                dcc(2, 1) = 0.d0
            end if
!
!               -- aire de surface appartenant au point de g.
            surfgp = poids
!
!               --  prise en compte de la dilatation thermique
            call coqgth(zi(imate), rela_comp, 'RIGI', ipg, ep, epsm, deps)
            do i = 1, 6
                epst(i) = epsm(i)+deps(i)
            end do
!
            call r8inir(36, 0.d0, dsidep, 1)
            call glrcmm(zi(imate), matr, ep, surfgp, pgl, &
                        epst, deps, dsig, ecr, delas, &
                        dsidep, zr(icarcr), codret)
!
            do i = 1, 3
                dsig(i) = dsig(i)*ep
                dsig(3+i) = dsig(3+i)*ep*ep/6.d0
            end do
!
            do i = 1, nbvari
                ecrp(i) = ecr(i)
            end do
!
            if (q4gg) then
                dsig(7) = dcc(1, 1)*dgam(1)
                dsig(8) = dcc(2, 2)*dgam(2)
            end if
!
            do i = 1, nbsig
                sig(i) = sig(i)+dsig(i)
            end do
!
        else if (rela_comp(1:7) .eq. 'GLRC_DM') then
            if (.not. lLinear) then
                do i = 1, nbvari
                    ecr(i) = zr(ivarim-1+icpv+i)
                end do
            end if
!
            call glrc_recup_mate(zi(imate), rela_comp, lLinear, ep, lambda=lambda, &
                                 deuxmu=deuxmu, lamf=lamf, deumuf=deumuf, &
                                 gt=gt, gc=gc, gf=gf, seuil=seuil, &
                                 alpha=alphaf, alfmc=alfmc, epsic=epsi_c, &
                                 epsiels=epsi_els, epsilim=epsi_lim, &
                                 is_param_opt_=is_param_opt, val_param_opt_=val_param_opt)

!
!               --  prise en compte de la dilatation thermique
            call coqgth(zi(imate), rela_comp, 'RIGI', ipg, ep, epsm, deps)
!
!               -- endommagement seulement
            call r8inir(36, 0.d0, dsidep, 1)
            call glrc_lc(epsm, deps, ecr, option, sig, &
                         ecrp, dsidep, lambda, deuxmu, lamf, &
                         deumuf, gt, gc, gf, seuil, &
                         alphaf, alfmc, zr(icarcr), &
                         epsi_c, epsi_els, epsi_lim, codret, &
                         ep, is_param_opt, val_param_opt, t2iu)
!
        else if (rela_comp(1:4) .eq. 'DHRC') then
!
            if (.not. lLinear) then
                do i = 1, nbvari
                    ecr(i) = zr(ivarim-1+icpv+i)
                end do
            end if
!
            call dhrc_recup_mate(zi(imate), rela_comp, a0, c0, &
                                 aa_t, ga_t, ab_, gb_, ac_, &
                                 gc_, aa_c, ga_c, cstseu)

!
!               --  prise en compte de la dilatation thermique
            call coqgth(zi(imate), rela_comp, 'RIGI', ipg, ep, epsm, deps)
!
!               -- endommagement couple glissement acier beton
            call r8inir(36, 0.d0, dsidep, 1)
            call dhrc_lc(epsm, deps, ecr, pgl, option, &
                         sig, ecrp, a0, c0, aa_t, &
                         ga_t, ab_, gb_, ac_, gc_, &
                         aa_c, ga_c, cstseu, zr(icarcr), codret, &
                         dsidep, ASTER_FALSE)
!
        else if (rela_comp(1:7) .eq. 'KIT_DDI') then
!
            if (.not. lLinear) then
                do i = 1, nbvari
                    ecr(i) = zr(ivarim-1+icpv+i)
                end do
            end if

!
!               --  Prise en compte de la dilatation thermique
!                   En realite, on ne sait pas encore faire ...
!                   La routine coqgth nous arretera s'il y a de la dilatation
            call coqgth(zi(imate), rela_comp, 'RIGI', ipg, ep, epsm, deps)
!               -- endommagement plus plasticite
            call r8inir(3, r8vide(), angmas, 1)
            rela_plas = compor(PLAS_NAME)
            call kit_glrc_dm_vmis(zi(imate), rela_plas, epsm, deps, ecr, &
                                  option, sigm, sig, ecrp, dsidep, &
                                  zr(icarcr), codret, t2iu)
        else
            call utmess('F', 'PLATE1_7', sk=rela_comp)
        end if
!
        if (.not. lLinear) then
            do i = 1, nbvari
                zr(ivarip-1+icpv+i) = ecrp(i)
            end do
        end if
!
        if (lVect) then
!
!               -- EFFORTS RESULTANTS (N, M ET Q)

            do i = 1, 3
                n(i) = sig(i)
                m(i) = sig(i+3)
            end do
!               -- PRISE EN COMPTE DE L'EXCENTREMENT DANS LE TERME DE FLEXION (VOIR AUSSI DKTNLI)

            if (rela_comp(1:4) .ne. 'ELAS' .and. excen .gt. 0.d0) then
                do i = 1, 3
                    m(i) = m(i)+excen*n(i)
                end do
            end if
            do i = 1, 2
                q(i) = sig(i+6)
            end do
!
!               -- CALCUL DE DIV(SIGMA) ET RECOPIE DE N ET M DANS 'PCONTPR'
!                  BTSIG = BTSIG + BFT*M + BMT*N + BCT*Q
!
            do k = 1, 3
                effint((ipg-1)*8+k) = n(k)
                effint((ipg-1)*8+3+k) = m(k)
            end do
!
            if (q4gg) then
                do k = 1, 2
                    effint((ipg-1)*8+6+k) = q(k)
                end do
            end if
!
            do ino = 1, nno
                do k = 1, 3
                    vecloc(1, ino) = vecloc(1, ino)+bm(k, 2*(ino-1)+1)*n(k)*poids
                    vecloc(2, ino) = vecloc(2, ino)+bm(k, 2*(ino-1)+2)*n(k)*poids
                    vecloc(3, ino) = vecloc(3, ino)+bf(k, 3*(ino-1)+1)*m(k)*poids
                    vecloc(5, ino) = vecloc(5, ino)+bf(k, 3*(ino-1)+2)*m(k)*poids
                    vecloc(4, ino) = vecloc(4, ino)-bf(k, 3*(ino-1)+3)*m(k)*poids
                end do
!
!                   -- PRISE EN COMPTE DU CISAILLEMENT
                if (q4gg) then
                    do k = 1, 2
                        vecloc(3, ino) = vecloc(3, ino)+bc(k, 3*(ino-1)+1)*q(k)*poids
                        vecloc(5, ino) = vecloc(5, ino)+bc(k, 3*(ino-1)+2)*q(k)*poids
                        vecloc(4, ino) = vecloc(4, ino)-bc(k, 3*(ino-1)+3)*q(k)*poids
                    end do
                end if
            end do
        end if
!
        if (lMatr) then
!           -- CALCUL DES MATRICES TANGENTES MATERIELLES (DM,DF,DMF)
            l = 0
            do i = 1, 3
                do j = 1, 3
                    l = l+1
                    dm(l) = dm(l)+poids*dsidep(j, i)
                    dmf(l) = dmf(l)+poids*dsidep(j, i+3)
                    df(l) = df(l)+poids*dsidep(j+3, i+3)
                end do
            end do
!
            if (q4gg) then
                dc(1, 1) = poids*dcc(1, 1)
                dc(2, 2) = dc(1, 1)
                dc(1, 2) = 0.d0
                dc(2, 1) = 0.d0
            end if
!
!               CALCUL DE LA MATRICE TANGENTE :
!               KTANG = KTANG + BFT*DF*BF + BMT*DM*BM + BMT*DMF*BF
!                             + BCT*DC*BC
!           -- MEMBRANE
            call utbtab('CUMUL', 3, 2*nno, dm, bm, &
                        work, memb)
!
!           -- FLEXION
            call utbtab('CUMUL', 3, 3*nno, df, bf, &
                        work, flex)
!
!           -- CISAILLEMENT
            if (q4gg) then
                call utbtab('CUMUL', 2, 3*nno, dc, bc, &
                            work, flex)
            end if

!           -- COUPLAGE
            call utctab('CUMUL', 3, 3*nno, 2*nno, dmf, &
                        bf, bm, work, mefl)
        end if
    end do
!
    if (.not. lLinear) then
!           -- PASSAGE DES EFFORTS GENERALISES AUX POINTS D'INTEGRATION
!              DU REPERE INTRINSEQUE AU REPERE LOCAL
!              STOCKAGE DES EFFORTS GENERALISES
!
        call dxefro(npg, t2iu, effint, efforp)
        do i = 1, nbcont*npg
            zr(icontp-1+i) = efforp(i)
        end do
    end if
!
    if (lVect) then
        call utpvlg(nno, 6, pgl, vecloc, zr(ivectu))
    end if
!
    if (lSigm) then
        zi(jcret) = codret
    end if
!
    if (lMatr) then
        if (t3g) then
            call dxtloc(flex, memb, mefl, ctor, matloc)
        else if (q4g) then
            call dxqloc(flex, memb, mefl, ctor, matloc)
        end if
        call jevech('PMATUUR', 'E', imatuu)
        call utpslg(nno, 6, pgl, matloc, zr(imatuu))
    end if
!
end subroutine
