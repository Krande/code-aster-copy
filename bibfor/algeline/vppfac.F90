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
subroutine vppfac(lmasse, masgen, vect, neq, nbvect, &
                  mxvect, masmod, facpar, masmoduni, inemodz)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterfort/dismoi.h"
#include "asterfort/get_equa_info.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jemarq.h"
#include "asterfort/mrmult.h"
#include "asterfort/pteddl.h"
#include "asterfort/wkvect.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsvpar.h"
#include "blas/ddot.h"
#include "asterfort/gettco.h"
!
    integer(kind=8) :: lmasse, neq, nbvect, mxvect
    real(kind=8) :: masgen(*), vect(neq, *)
    real(kind=8) :: masmod(mxvect, *), facpar(mxvect, *), masmoduni(mxvect, *)
    real(kind=8), optional ::  inemodz(mxvect, *)
!     CALCUL DES PARAMETRES MODAUX :
!            FACTEUR DE PARTICIPATION ET MASSE MODALE UNITAIRE
!     ------------------------------------------------------------------
! IN  LMASSE : IS : DESCRIPTEUR NORMALISE DE LA MATRICE DE MASSE
!     ------------------------------------------------------------------
!     PRESUME L'EXECUTION PREALABLE DE VPPGEN : CALCUL DES PARAMETRES
!     MODAUX.
!     ------------------------------------------------------------------
!
!
    integer(kind=8) :: lddl, laux1, laux2, iddl, ia, ieq, ivect, mxddl, iadpar(1), l1, ibid
    integer(kind=8) :: iddl2, iddl3, ia2, ia3, ncdg, nbmrig, nb_nodes_mesh, kcoor, node
    parameter(mxddl=6)
    character(len=8) :: nomddl(mxddl), basemo, k8b, typeq, mesh
    character(len=14) :: nume
    character(len=16) :: nompar(3), typmas, typbas
    character(len=19) :: masse
    character(len=24) :: posddl, vecau1, vecau2
    real(kind=8) :: rmin, rmax, raux, rval, cdg(3), coorno(3), massTotDirUnit
    aster_logical :: gene
    character(len=24), pointer :: refn(:) => null()
    real(kind=8) :: rundef, epsi
    blas_int :: b_incx, b_incy, b_n
!     ------------------------------------------------------------------
    data nomddl/'DX      ', 'DY      ', 'DZ      ',&
     &              'DRX     ', 'DRY     ', 'DRZ     '/
    data nompar/'FACT_PARTICI_DX', 'FACT_PARTICI_DY', 'FACT_PARTICI_DZ'/
!
!     ------------------------------------------------------------------
    data posddl/'&&VPPFAC.POSITION.DDL'/
    data vecau1/'&&VPPFAC.VECTEUR.AUX1'/
    data vecau2/'&&VPPFAC.VECTEUR.AUX2'/
!
    rundef = r8vide()
    epsi = r8prem()
!
!     ------------------------------------------------------------------
!     ----------------- CREATION DE VECTEURS DE TRAVAIL ----------------
!     ------------------------------------------------------------------
!
    call wkvect(vecau1, 'V V R', neq, laux1)
    call wkvect(vecau2, 'V V R', neq, laux2)
    call wkvect(posddl, 'V V I', neq*mxddl, lddl)
!
!
!   recuperation des coordonnées du centre de gravité pour le calcul
!   des inerties effectives
    nbmrig = 3
    call getvr8('', 'CENTRE', nbval=3, vect=cdg, nbret=ncdg)
    if (ncdg .ne. 0 .and. present(inemodz)) then
        nbmrig = 6
    end if
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ------ RECUPERATION DES POSITIONS DU VECTEUR EXCITATION : --------
!     ---- VECTEUR DE 1 DANS LA DIRECTION iddl DANS LE CAS PHYSIQUE ----
!     --- FACTEURS DE PARTICIPATION DES MODES DANS LE CAS GENERALISE ---
!     ------------------------------------------------------------------
!
    call jemarq()
    masse = zk24(zi(lmasse+1)) (1:19)
    call dismoi('NOM_NUME_DDL', masse, 'MATR_ASSE', repk=nume)
    call gettco(masse, typmas)
!
    rmin = 100.d0*r8miem()
    rmax = sqrt(r8maem())
!
! DETERMINATION DU CAS : BASE PHYSIQUE OU BASE GENERALISEE
    gene = .false.
    if (typmas(1:14) .eq. 'MATR_ASSE_GENE') then
! SI MATR_ASSE_GENE : BASE GENERALISEE
        call jeveuo(nume(1:14)//'.NUME.REFN', 'L', vk24=refn)
        basemo = refn(1) (1:8)
! SAUF SI NUME.REFN POINTE VERS UN MODELE_GENE ET NON VERS UNE BASE
! ALORS CAS DE LA SSD, TRAITE COMME UN MODE_MECA CLASSIQUE
        call gettco(basemo, typbas)
        if (typbas(1:14) .eq. 'MODE_MECA') then
            gene = .true.
        end if
    elseif (nbmrig .eq. 6) then
        call dismoi('NOM_MAILLA', nume, 'NUME_DDL', repk=mesh)
        call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_nodes_mesh)
        call jeveuo(mesh//'.COORDO    .VALE', 'L', kcoor)
    end if

    if (.not. gene) then
        call pteddl('NUME_DDL', nume, mxddl, nomddl, neq, &
                    tabl_equa=zi(lddl))
    end if

    do iddl = 1, nbmrig
        if (gene) then
            if (iddl .gt. 3) exit
            do ieq = 1, neq
                call rsvpar(basemo, 1, nompar(iddl), ibid, rundef, &
                            k8b, l1)
                if (l1 .eq. 100) then
                    zr(laux1+ieq-1) = 0.d0
                else
                    call rsadpa(basemo, 'L', 1, nompar(iddl), ieq, &
                                0, tjv=iadpar)
                    zr(laux1+ieq-1) = zr(iadpar(1))
                end if
! SECURITE SI ON EST PASSE PAR DES MODES HETERODOXES AVEC FACTEURS DE PARTICIPATIONS HERETIQUES
            end do
        else
            if (iddl .eq. 4) then
                iddl2 = 2
                iddl3 = 3
            elseif (iddl .eq. 5) then
                iddl2 = 3
                iddl3 = 1
            elseif (iddl .eq. 6) then
                iddl2 = 1
                iddl3 = 2
            else
                iddl2 = 0
                iddl3 = 0
            end if
            ia = (iddl-1)*neq
            ia2 = (iddl2-1)*neq
            ia3 = (iddl3-1)*neq
!
!           calcul des vecteurs de déplacement unitaire U_d
!           voir R5.01.30(Vecteur déplacement unitaire)
            if (iddl .le. 3) then
!               en translation
                do ieq = 1, neq
                    zr(laux1+ieq-1) = zi(lddl+ia+ieq-1)
                end do
            else
!               en rotation
                do ieq = 1, neq
                    zr(laux1+ieq-1) = 0.d0
                    call get_equa_info(nume, ieq, typeq, nume_nodez=node)
                    if (typeq(1:1) .ne. 'A') cycle
                    coorno(1) = zr(kcoor-1+3*(node-1)+1)
                    coorno(2) = zr(kcoor-1+3*(node-1)+2)
                    coorno(3) = zr(kcoor-1+3*(node-1)+3)
                    if (iddl2 .gt. 0) then
!                       il y a au plus un zi(lddl+iaX non nul
!                       car une équation ne correspond qu'à une composante
                        zr(laux1+ieq-1) = zi(lddl+ia+ieq-1) &
                                          -zi(lddl+ia2+ieq-1) &
                                          *(coorno(iddl3)-cdg(iddl3)) &
                                          +zi(lddl+ia3+ieq-1) &
                                          *(coorno(iddl2)-cdg(iddl2))
                    end if
                end do
            end if
        end if
!
!     ------------------------------------------------------------------
!     ----------- CALCUL DE  FREQ * MASSE * UNITAIRE_DIRECTION (Ud) ----
!     ------------------------------------------------------------------
!       M*Ud
        call mrmult('ZERO', lmasse, zr(laux1), zr(laux2), 1, &
                    .false._1)
!       masse ou inertie de rotation globale selon la direction
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
!       Ud^t*M*Ud
        massTotDirUnit = ddot(b_n, zr(laux1), b_incx, zr(laux2), b_incy)
!
        do ivect = 1, nbvect
!           Vectp(freq)^t*M*Ud
            rval = ddot(b_n, vect(1, ivect), b_incx, zr(laux2), b_incy)
            raux = masgen(ivect)
            if ((abs(raux) .lt. rmin) .or. (abs(rval) .gt. rmax)) then
                if (iddl .le. 3) then
                    masmod(ivect, iddl) = rmax
                    facpar(ivect, iddl) = rmax
                    if (massTotDirUnit .gt. epsi) then
                        masmoduni(ivect, iddl) = rmax/massTotDirUnit
                    else
                        masmoduni(ivect, iddl) = 0.d0
                    end if
                else
                    iddl2 = iddl-3
                    ! INER_EFFE_D*
                    inemodz(ivect, iddl2) = rmax
                    ! INER_EFFE_UN_D*
                    if (massTotDirUnit .gt. epsi) then
                        inemodz(ivect, iddl) = rmax/massTotDirUnit
                    else
                        inemodz(ivect, iddl) = 0.d0
                    end if
                end if

            else
                raux = rval/raux
                if (iddl .le. 3) then
                    masmod(ivect, iddl) = rval*raux
                    facpar(ivect, iddl) = raux
                    if (massTotDirUnit .gt. epsi) then
                        masmoduni(ivect, iddl) = rval*raux/massTotDirUnit
                    else
                        masmoduni(ivect, iddl) = 0.d0
                    end if
                else
                    iddl2 = iddl-3
                    ! INER_EFFE_D*
                    inemodz(ivect, iddl2) = rval*raux
                    ! INER_EFFE_UN_D*
                    inemodz(ivect, iddl) = rval*raux/massTotDirUnit
                    if (massTotDirUnit .gt. epsi) then
                        inemodz(ivect, iddl) = rval*raux/massTotDirUnit
                    else
                        inemodz(ivect, iddl) = 0.d0
                    end if
                end if
            end if
        end do
    end do
!
!     ------------------------------------------------------------------
!     ----------------- DESTRUCTION DES VECTEURS DE TRAVAIL ------------
!     ------------------------------------------------------------------
!
    call jedetr(posddl)
    call jedetr(vecau1)
    call jedetr(vecau2)
!
    call jedema()
end subroutine
