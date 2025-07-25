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

subroutine te0537(option, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!                   POU_D_EM : POUTRE MULTIFIBRE EULER BERNOULLI
!
!   Calcul de l'option EPSI_ELGA
!       - déformations dans les fibres (sous pts de gauss) à partir des déplacements
!   Calcul de l'option SIEF_ELGA
!       - contrainte dans les fibres (sous pts de gauss) comportement linéaire
!   Calcul de l'option STRX_ELGA
!
! --------------------------------------------------------------------------------------------------
!                   POU_D_TGM : POUTRE MULTIFIBRE TIMOSHENKO
!
!   Calcul de l'option EPSI_ELGA
!       - déformations dans les fibres (sous pts de gauss) à partir des déplacements
! --------------------------------------------------------------------------------------------------
!
    implicit none
    character(len=16) :: option, nomte
!
#include "jeveux.h"
#include "MultiFiber_type.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jsd1ff.h"
#include "asterfort/lonele.h"
#include "asterfort/matela.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/pmavec.h"
#include "asterfort/pmfdef.h"
#include "asterfort/pmfdge.h"
#include "asterfort/pmfinfo.h"
#include "asterfort/pmfitx.h"
#include "asterfort/pmfmats.h"
#include "asterfort/pmfpti.h"
#include "asterfort/pmfrig.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/vecma.h"
#include "asterfort/verifm.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jcont, lorien, jdepl, imate, nno, nc, i, iret
    integer(kind=8) :: ip, ipos, istrxr, ipos1, ipos2, nbfig, ig, icp, isdcom
    integer(kind=8) :: codres(2), ncomp
    integer(kind=8) :: npg, ndim, nnoel, nnos, ipoids, ivf
    integer(kind=8) :: jacf, jtab(7)
    parameter(nno=2)
    real(kind=8) :: ul(14), pgl(3, 3), dege(6), xl, e, nu
    real(kind=8) :: g
    real(kind=8) :: casect(6)
    real(kind=8) :: coa, cob, ex12, ex13
    real(kind=8) :: b(4), gg, xi, wi, valres(2), alpha
    real(kind=8) :: klv(78), klc(12, 12), effo(12)
    character(len=8) :: materi
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: ch16, nomres(2)
!
    integer(kind=8) :: nbfibr, nbgrfi, tygrfi, nbcarm, nug(10)
    real(kind=8) :: a, xiy, xiz, alfay, alfaz, phiy, phiz, ey, ez
    real(kind=8) :: epsthe, temp, d1b(7, 14), eps(7)
    integer(kind=8) :: lmater, itemp, j
    character(len=4) :: fami
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_cara = 9
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara), nomat
    data noms_cara/'A1', 'IY1', 'IZ1', 'AY1', 'AZ1', 'EY1', 'EZ1', 'EY2', 'EZ2'/
!
! --------------------------------------------------------------------------------------------------
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nnoel, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf)
    ASSERT(nno .eq. nnoel)
! --------------------------------------------------------------------------------------------------
!   Récupération des caractéristiques des fibres
    call pmfinfo(nbfibr, nbgrfi, tygrfi, nbcarm, nug, jacf=jacf)
!
    alpha = 0.d0
    call jevech('PCAORIE', 'L', lorien)
    call jevech('PDEPLAR', 'L', jdepl)
!   Recuperation des coordonnees des noeuds
    xl = lonele()
    nc = 6
    if (nomte .eq. 'MECA_POU_D_TGM') nc = 7
!
!   Recuperation des orientations
    call matrot(zr(lorien), pgl)
!   Passage des deplacements dans le repere local
    call utpvgl(nno, nc, pgl, zr(jdepl), ul)
!   Nombre de composantes des champs PSTRX par points de gauss
    ncomp = 18
!
    jcont = -1
    if (nomte .eq. 'MECA_POU_D_EM') then
        if (option .eq. 'EPSI_ELGA') then
            call tecach('OOO', 'PDEFOPG', 'E', iret, nval=7, itab=jtab)
            jcont = jtab(1)
        else if (option .eq. 'SIEF_ELGA') then
            call jevech('PMATERC', 'L', imate)
            call tecach('OOO', 'PCONTRR', 'E', iret, nval=7, itab=jtab)
            jcont = jtab(1)
        else if (option .eq. 'STRX_ELGA') then
            call jevech('PMATERC', 'L', imate)
            call jevech('PSTRXRR', 'E', istrxr)
!           Si excentrecement calcul de alpha appel intégration sur section et calcul G torsion
            call pmfitx(zi(imate), 1, casect, g)
            if ((casect(2) .ge. r8prem()) .or. (casect(3) .ge. r8prem())) then
                coa = 3.0d0/2.0d0/xl
                cob = 3.0d0/4.0d0
                ex13 = casect(2)/casect(1)
                ex12 = casect(3)/casect(1)
                alpha = coa*ex13*ul(2)-coa*ex12*ul(3) &
                        +cob*ex12*ul(5)+cob*ex13*ul(6) &
                        -coa*ex13*ul(8)+coa*ex12*ul(9) &
                        +cob*ex12*ul(11)+cob*ex13*ul(12)
            else
                alpha = 0.0d0
            end if
        else
            ch16 = option
            call utmess('F', 'ELEMENTS2_47', sk=ch16)
        end if
!
!       si option EPSI_ELGA ou SIEF_ELGA boucle sur les points de gauss
        if (option .ne. 'STRX_ELGA') then
!           alpha modes incompatibles
            call jevech('PSTRXRR', 'L', istrxr)
            alpha = zr(istrxr-1+15)
            do ip = 1, npg
!               Matrice B puis DEGE puis déformations sur les fibres
                call pmfpti(ip, zr(ipoids), zr(ivf), xl, xi, wi, b, gg)
                call pmfdge(b, gg, ul, alpha, dege)
                ipos = jcont+nbfibr*(ip-1)
                call pmfdef(tygrfi, nbfibr, nbcarm, zr(jacf), dege, zr(ipos))
            end do
        end if
!       Si EPSI_ELGA : jcont est l'adresse PDEFORR, on sort
        if (option .eq. 'SIEF_ELGA') then
!           Si option sief_elga on continue
!           Récupération des différents matériaux dans SDCOMP dans COMPOR
            call jevech('PCOMPOR', 'L', vk16=compor)
            call jeveuo(compor(MULTCOMP), 'L', isdcom)
!           boucle sur les groupes de fibre
            ipos1 = jcont-1
            ipos2 = ipos1+nbfibr
            do ig = 1, nbgrfi
                icp = isdcom-1+(nug(ig)-1)*MULTI_FIBER_SIZEK
!               nombre de fibres de ce groupe
                read (zk24(icp+MULTI_FIBER_NBFI), '(I24)') nbfig
                materi = zk24(icp+MULTI_FIBER_MATER) (1:8)
!               On multiplie par E (constant sur le groupe)
                nomres(1) = 'E'
                nomres(2) = 'NU'
                call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                            materi, 'ELAS', 0, ' ', [0.d0], &
                            2, nomres, valres, codres, 1)
                e = valres(1)
                nu = valres(2)
!               On multiplie les zr(jcont) (déformations) par E pour avoir des contraintes
                do i = 1, nbfig
                    zr(ipos1+i) = zr(ipos1+i)*e
                    zr(ipos2+i) = zr(ipos2+i)*e
                end do
                ipos1 = ipos1+nbfig
                ipos2 = ipos2+nbfig
            end do
        end if
!
        if (option .eq. 'STRX_ELGA') then
!           Calcul des efforts généralisés (cas élastique) on fait kele*ul
            call pmfrig(nomte, zi(imate), klv)
            call vecma(klv, 78, klc, 12)
            call pmavec('ZERO', 12, klc, ul, effo)
            do ip = 1, npg
                zr(istrxr-1+ncomp*(ip-1)+1) = effo(6*(ip-1)+1)
                zr(istrxr-1+ncomp*(ip-1)+2) = effo(6*(ip-1)+2)
                zr(istrxr-1+ncomp*(ip-1)+3) = effo(6*(ip-1)+3)
                zr(istrxr-1+ncomp*(ip-1)+4) = effo(6*(ip-1)+4)
                zr(istrxr-1+ncomp*(ip-1)+5) = effo(6*(ip-1)+5)
                zr(istrxr-1+ncomp*(ip-1)+6) = effo(6*(ip-1)+6)
                zr(istrxr-1+ncomp*(ip-1)+15) = alpha
            end do
        end if
!   MECA_POU_D_TGM : EPSI_ELGA
    else
!
        call tecach('OOO', 'PDEFOPG', 'E', iret, nval=7, itab=jtab)
        jcont = jtab(1)

!       RECUPERATION DES CARACTERISTIQUES GENERALES DES SECTIONS
        call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
        a = vale_cara(1)
        xiy = vale_cara(2)
        xiz = vale_cara(3)
        alfay = vale_cara(4)
        alfaz = vale_cara(5)
        ey = (vale_cara(6)+vale_cara(8))/2.d0
        ez = (vale_cara(7)+vale_cara(9))/2.d0

!       CARACTERISTIQUES MATERIAUX
        call jevech('PMATERC', 'L', lmater)
        call pmfmats(nomat)
!
        call verifm(fami, npg, 1, '+', zi(lmater), epsthe, iret)
        itemp = 0
        if (iret .eq. 0) itemp = 1
!
        call moytem(fami, npg, 1, '+', temp, iret)
        call matela(zi(lmater), nomat, itemp, temp, e, nu)
!
        g = e/(2.0d0*(1.0d0+nu))
        phiy = e*xiz*12.d0*alfay/(xl*xl*g*a)
        phiz = e*xiy*12.d0*alfaz/(xl*xl*g*a)
!
!       PASSAGE DE G (CENTRE DE GRAVITE) A C (CENTRE DE TORSION)
        do i = 1, 2
            ul(7*(i-1)+2) = ul(7*(i-1)+2)-ez*ul(7*(i-1)+4)
            ul(7*(i-1)+3) = ul(7*(i-1)+3)+ey*ul(7*(i-1)+4)
        end do
!
        do ip = 1, npg
!           calcul des déformations généralisées
            call jsd1ff(ip, xl, phiy, phiz, d1b)
            eps(1:nc) = 0.0d0
            do i = 1, nc
                do j = 1, 2*nc
                    eps(i) = eps(i)+d1b(i, j)*ul(j)
                end do
            end do
!           calcul des deformations sur les fibres
            ipos = jcont+nbfibr*(ip-1)
            call pmfdef(tygrfi, nbfibr, nbcarm, zr(jacf), eps, zr(ipos))
        end do
    end if
end subroutine
