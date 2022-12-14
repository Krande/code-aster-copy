! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine te0442(option, nomte)
    !
    implicit none
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/coqrep.h"
#include "asterfort/cylrep.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxsiro.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/normev.h"
#include "asterfort/tecach.h"
#include "asterfort/tpsivp.h"
#include "asterfort/utmess.h"
#include "blas/dgemm.h"
    !
    character(len=16), intent(in) :: option, nomte
!......................................................................
    !
!    - FONCTION REALISEE: CHANGEMENT DE REPERE POUR LES PLAQUES
    !
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                        'REPE_TENS'  :  TENSEURS
!                        'REPE_GENE'  :  QUANTITES GENERALISEES
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
    !
    !
    integer :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano, iret(4)
    integer :: jgeom, jin, jout, jang, np, itab(7), iret1, iret2, nbsp
    integer :: jcara, ncmp, vali(2)
    real(kind=8) :: alpha, beta, c, s
    real(kind=8), dimension(3, 3) :: pig, pgcyl, picyl
    real(kind=8), dimension(2, 2) :: t2iu1, t2ui1
    real(kind=8), dimension(2, 2) :: t2iu2, t2ui2
    integer, parameter :: nptmax=9, nspmax=162
    real(kind=8), dimension(3) :: axe_z, orig, x, xsp, xbary
    real(kind=8) :: rep
    character(len=4) :: fami
    character(len=8) :: pain, paout
    integer :: ipt, ino, joff, type_pt, ipaxe, ipaxe2
    real(kind=8) :: a, b, xnorm, epais, excen, zic, hicou
    integer, parameter :: pt_gauss = 1, pt_noeud = 2
    integer :: jnbspi, nbcou, icou, isp
    !
    if (option .ne. 'REPE_TENS' .and. option .ne. 'REPE_GENE') then
! OPTION DE CALCUL INVALIDE
        ASSERT(.false.)
    endif
    !
    if (option .eq. 'REPE_TENS') then
        ncmp=6
        call tecach('ONO', 'PCOGAIN', 'L', iret(1), nval=7,&
                    itab=itab)
        call tecach('ONO', 'PCONOIN', 'L', iret(2), nval=7,&
                    itab=itab)
        call tecach('ONO', 'PDEGAIN', 'L', iret(3), nval=7,&
                    itab=itab)
        call tecach('ONO', 'PDENOIN', 'L', iret(4), nval=7,&
                    itab=itab)
        iret1 = iret(1) + iret(2) + iret(3) + iret(4)
        ASSERT(iret1.eq.6)
        !
        if (iret(1) .eq. 0) then
            pain = 'PCOGAIN'
            paout = 'PCOGAOUT'
        else if (iret(2).eq.0) then
            pain = 'PCONOIN'
            paout = 'PCONOOUT'
        else if (iret(3).eq.0) then
            pain = 'PDEGAIN'
            paout = 'PDEGAOUT'
        else if (iret(4).eq.0) then
            pain = 'PDENOIN'
            paout = 'PDENOOUT'
        endif
        !
    else if (option.eq.'REPE_GENE') then
        ncmp=8
        call tecach('ONO', 'PEFGAIN', 'L', iret(1), nval=7,&
                    itab=itab)
        call tecach('ONO', 'PEFNOIN', 'L', iret(2), nval=7,&
                    itab=itab)
        call tecach('ONO', 'PDGGAIN', 'L', iret(3), nval=7,&
                    itab=itab)
        call tecach('ONO', 'PDGNOIN', 'L', iret(4), nval=7,&
                    itab=itab)
        iret1 = iret(1) + iret(2) + iret(3) + iret(4)
        ASSERT(iret1.eq.6)
        !
        if (iret(1) .eq. 0) then
            pain = 'PEFGAIN'
            paout = 'PEFGAOUT'
        else if (iret(2).eq.0) then
            pain = 'PEFNOIN'
            paout = 'PEFNOOUT'
        else if (iret(3).eq.0) then
            pain = 'PDGGAIN'
            paout = 'PDGGAOUT'
        else if (iret(4).eq.0) then
            pain = 'PDGNOIN'
            paout = 'PDGNOOUT'
        endif
    endif
! Localisation des champs: noeuds/pts de Gauss
    if (pain(4:5) .eq. 'NO') then
        fami = 'NOEU'
        type_pt = pt_noeud
    else if (pain(4:5).eq.'GA') then
        fami = 'RIGI'
        type_pt = pt_gauss
    endif
! Infos sur les noeuds et points de Gauss de l'??l??ment
    call elrefe_info(fami=fami,ndim=ndim,nno=nno,nnos=nnos,&
  npg=npg,jpoids=ipoids,jvf=ivf,jdfde=idfdx,jgano=jgano)
! Nombre de points en fonction de la localisation des champs
    if (pain(4:5) .eq. 'NO') then
        np = nno
    else if (pain(4:5).eq.'GA') then
        np = npg
    endif
    ASSERT(np.le.nptmax)
    !
! Adresses des champs locaux
    call jevech('PGEOMER', 'L', jgeom)
    call jevech('PCACOQU', 'L', jcara)
    call jevech('PANGREP', 'L', jang)
    call jevech(pain, 'L', jin)
    call jevech(paout, 'E', jout)
    call tecach('OOO', pain, 'L', iret2, nval=7,&
                itab=itab)
! Nombre de sous-points dans le champ local
    nbsp = itab(7)
    if ((nbsp.ne.1) .and. (mod(nbsp,3).ne.0)) then
        call utmess('F', 'ELEMENTS5_54', si=nbsp)
    endif
    vali(1)=nspmax
    vali(2)=nbsp
    if (nbsp .gt. nspmax) then
        call utmess('F', 'ELEMENTS5_4', ni=2, vali=vali)
    endif
!
    !
! Pig: matrice de passage du rep??re intrins??que au rep??re global
    if (nno .eq. 3) then
! pour une maille triangle
        call dxtpgl(zr(jgeom), pig)
    else if (nno.eq.4) then
! pour une maille quadrangle
        call dxqpgl(zr(jgeom), pig, 'S', iret1)
    endif
    !
! Param??tres de la coque :
    epais = zr(jcara-1+1)
    alpha = zr(jcara-1+2) * r8dgrd()
    beta = zr(jcara-1+3) * r8dgrd()
    excen = zr(jcara-1+5)
!
! T2iu1 : matrice (2x2) de passage du rep??re intrins??que au rep??re utilisateur 1 de
! la vari??t?? (d??fini par l'utilisateur dans caraelem)
! T2ui1 = tiu1^T
    call coqrep(pig, alpha, beta, t2iu1, t2ui1,&
                c, s)
!  Type de changement de rep??re
    rep = zr(jang-1+3)
    !
    if (nint(rep) == 0) then
        !
! --- PASSAGE DES CONTRAINTES DU REPERE LOCAL 1
! --- A L'ELEMENT AU REPERE INTRINSEQUE DE LA COQUE
!     ---------------------------------------
        if (option .eq. 'REPE_TENS') then
            call dxsiro(np*nbsp, t2ui1, zr(jin), zr(jout))
        else if (option.eq.'REPE_GENE') then
            call dxefro(np, t2ui1, zr(jin), zr(jout))
        endif
        !
! --- CALCUL DES MATRICES DE PASSAGE DU CHGT DE REPERE
        alpha = zr(jang-1+1) * r8dgrd()
        beta = zr(jang-1+2) * r8dgrd()
        call coqrep(pig, alpha, beta, t2iu2, t2ui2,&
                    c, s)
        !
! ---   PASSAGE DES QUANTITES DU REPERE INTRINSEQUE
! ---   A L'ELEMENT AU REPERE LOCAL 2 DE LA COQUE
        if (option .eq. 'REPE_TENS') then
            call dxsiro(np*nbsp, t2iu2, zr(jout), zr(jout))
        else if (option.eq.'REPE_GENE') then
            call dxefro(np, t2iu2, zr(jout), zr(jout))
        endif
        !
! --- PASSAGE DES CONTRAINTES DU REPERE INTRINSEQUE
! --- A L'ELEMENT AU REPERE LOCAL 1 DE LA COQUE
!     REPERE = 'COQUE_INTR_UTIL'
!     ---------------------------------------
    else if (nint(rep)== 1) then
        if (option .eq. 'REPE_TENS') then
            call dxsiro(np*nbsp, t2iu1, zr(jin), zr(jout))
        else if (option.eq.'REPE_GENE') then
            call dxefro(np, t2iu1, zr(jin), zr(jout))
        endif
        !
! --- PASSAGE DES CONTRAINTES DU REPERE LOCAL 1
! --- A L'ELEMENT AU REPERE INTRINSEQUE DE LA COQUE
!     REPERE = 'COQUE_UTIL_INTR'
!     ---------------------------------------
    else if (nint(rep) == 2) then
        if (option .eq. 'REPE_TENS') then
            call dxsiro(np*nbsp, t2ui1, zr(jin), zr(jout))
        else if (option.eq.'REPE_GENE') then
            call dxefro(np, t2ui1, zr(jin), zr(jout))
        endif
        !
    else if (nint(rep) == 3) then
!   PASSAGE DES CONTRAINTES DU REPERE LOCAL 1
!   A L'ELEMENT AU REPERE CYLINDRIQUE
!   REPERE ='COQUE_UTIL_CYL'
        if (option .eq. 'REPE_TENS') then
! Nombre de couches d??clar??es dans le caraelem.
            call jevech('PNBSP_I', 'L', jnbspi)
            nbcou = zi(jnbspi)
! V??rification de coh??rence :
! On v??rifie que le nb de sous-points dans le champ local
! est bien ??gal au nombre de sous-points (3 * nbcou) pr??vus
! pour la coque (i.e. qu'on n'a pas extrait un champ sur un
! sous point, ce qui rend le calcul de la cote des sous-points
! impossible)
            ASSERT ( nbsp /= np * 3 * nbcou )
            hicou = epais/nbcou
! D??finition du rep??re
            axe_z(:)= zr(jang-1+4: jang-1+6)
            call normev(axe_z, xnorm)
! Si l'axe z n'est pas initialis?? dans la carte, on ne veut pas changer de rep??re
! => on ne fait rien
            if (xnorm > r8prem()) then
                orig(:) = zr(jang-1+7: jang-1+9)
! Passage des contraintes du rep??re local 1 (rep??re utilisateur d??fini par
! cara_elem) au rep??re intrins??que pour tous les points et sous-points
! matrice de passage t2ui1
                call dxsiro(np*nbsp, t2ui1, zr(jin), zr(jout))
! Passage des contraintes du rep??re intrins??que au rep??re local cylindrique
                joff = 0
                do ipt = 1, np
! -- Calcul des coordonn??es du point
                    x(:) = 0.d0
                    if (type_pt == pt_gauss) then
!  Le point est un point de Gauss
                        do ino = 1, nno
                            x(1) = x(1) + zr(jgeom+3* (ino-1)-1+1)*zr(ivf+ (ipt-1)*nno+ ino-1)
                            x(2) = x(2) + zr(jgeom+3* (ino-1)-1+2)*zr(ivf+ (ipt-1)*nno+ ino-1)
                            x(3) = x(3) + zr(jgeom+3* (ino-1)-1+3)*zr(ivf+ (ipt-1)*nno+ ino-1)
                        enddo
                    else if (type_pt == pt_noeud) then
!  Le point est un noeud de l'??l??ment
                        x(1) = zr(jgeom+3* (ipt-1)-1+1)
                        x(2) = zr(jgeom+3* (ipt-1)-1+2)
                        x(3) = zr(jgeom+3* (ipt-1)-1+3)
                    endif
                    !
! Boucle sur les couches
!
                    do icou = 1, nbcou
                        !
!  Boucle sur les sous-points dans l'??paisseur de la couche
                        !
                        do isp = 1, 3
! Calcul de la cote du sous-point par rapport
! ?? la surface moyenne de la coque
                            zic = excen - epais/2.d0 + (icou-1)*hicou
                            if (isp .eq. 1) then
                                zic = zic
                            else if (isp.eq.2) then
                                zic = zic + hicou/2.d0
                            else
                                zic = zic + hicou
                            endif
!-- Coordonn??es du sous-point (pig(3,:) est la normale ?? l'??l??ment)
                            xsp(:) = x(:) + zic * pig(3,:)
!-- Calcul de la matrice de passage du rep??re global au rep??re local cylindrique
                            ipaxe = 0
                            call cylrep(ndim, xsp, axe_z, orig, pgcyl,&
                                        ipaxe)
                            if (ipaxe > 0) then
! le point est sur l'axe du repere cylindrique, on essaie de se placer au
! centre de gravit?? pour calculer le rep??re
                                xbary(:) = 0
                                do ino = 1, nno
                                    xbary(1)= xbary(1) + zr(jgeom+3* (ino-1)-1+1)
                                    xbary(2)= xbary(2) + zr(jgeom+3* (ino-1)-1+2)
                                    xbary(3)= xbary(3) + zr(jgeom+3* (ino-1)-1+3)
                                enddo
                                xbary(:)=xbary(:)/nno
                                ipaxe2 = 0
                                call cylrep(ndim, xbary, axe_z, orig, pgcyl,&
                                            ipaxe2)
                                if (ipaxe2 > 0) then
                                    call utmess('A', 'ALGORITH2_13')
                                endif
                            endif
! picyl: matrice de passage du rep??re intrins??que au rep??re cylindrique
! picyl = pig*pgcyl
                            a=1.d0
                            b=0.d0
                            call dgemm('N', 'N', 3, 3, 3,&
                                       a, pig(1, 1), 3, pgcyl(1, 1), 3,&
                                       b, picyl(1, 1), 3)
! Appliquer le changement de base
                            call tpsivp(picyl, zr(jout+joff-1+1:jout+joff-1+ncmp))
! Mettre ?? jour l'offset (un tenseur 3D sym??trique a 6 composantes)
                            joff = joff + ncmp
                        enddo
                    enddo
                enddo
            else
! on ne veut pas changer de rep??re, on recopie le champ
                zr(jout:jout+np*nbsp*ncmp) = zr(jin:jin+np*nbsp*ncmp)
            endif
        else
            ASSERT(.false.)
        endif
    endif
end subroutine te0442
