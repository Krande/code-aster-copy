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
subroutine libint(imped, nume91, nbint, lisint, nbeq1)
    implicit none
!
!---------------------------------------------------------------C
!--       ROUTINE XXXXX3        M. CORUS - AOUT 2011          --C
!--       CONSTRUCTION DE LA MATRICE LIBRE A L'INTERFACE      --C
!--                                                           --C
!--   REALISE APRES TEST DE LA METHODE DECRITE DANS LA DOC R3 --C
!--   (QUI MARCHE PAS) CONSSITANT A JOUER SUR LES COEFF DES   --C
!--   MULTIPLICATEURS DE LAGRANGE                             --C
!--                                                           --C
!--                                                           --C
!--       METHODE PERSO : ON VIRE LES LAGRANGES :             --C
!--          K(IND_LAG,IND_LAG)=IDENTITE*COEFF_LAGRANGE       --C
!--          K(IND_INTERF,IND_LAG)=0                          --C
!--  ET DONC K(IND_LAG,IND_INTERF)=0                          --C
!--                                                           --C
!--  AVANTAGES : ON CONSERVE LA NUMROTATION                   --C
!--              ON PEUT "MIXER" DES RESULTATS                --C
!--                  INTERFACE LIBRE / INTERFACE FIXE         --C
!--                                                           --C
!--  INCONVENIENT : ON RISQUE DE FAIRE N'IMPORTE QUOI SI      --C
!--                 ON NE FAIT PAS UN PEU ATTENTION           --C
!--                                                           --C
!-- RECOMMANDATION : N'UTILISER CETTE ROUTINE QUE SUR UNE     --C
!--                  COPIE DE LA MATICE DE RAIDEUR INITIALE   --C
!--                                                           --C
!--                                                           --C
!-- APRES TEST DE LA METHODE DECRITE DANS LA DOC R3           --C
!-- (QUI MARCHE PAS) CONSSITANT A JOUER SUR LES COEFF DES     --C
!-- MULTIPLICATEURS DE LAGRANGE                               --C
!--                                                           --C
!---------------------------------------------------------------C
!--   VARIABLES E/S  :
!--   IMPED    /IN/  : NOM K19 DE LA MATRICE DE RAIDEUR
!--   NUME91   /IN/  : NOM DU NUME_DDL ASSOCIE
!--   NBINT    /IN/  : NOMBRE D'INTERFACE DONT IL FAUT LIBERER
!--                       LES LAGRANGES
!--   LISINT   /IN/  : LISTE DES NOMS D'INTERFACES PERMETTANT DE
!--                       RECUPERER LES DDL CONCERNES PAR LA LIBERAtION
!--                       VOIR DANS LE CODE ET DANS OP0091 POUR UNE
!--                       UTILISATION DANS UN AUTRE CADRE
!--   NBEQ1    /IN/  : NB DE DDL DE LA MATRICE
!
!
!
!
!
#include "jeveux.h"
#include "asterfort/ddllag.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
!
    character(len=19) :: imped, nume91
    character(len=24) :: indin1
    integer(kind=8) :: j1, k1, l1, m1, n1, nbeq1, llint1, nbddl1, lintf, nbint, lklibr
    integer(kind=8) :: lag1, lag2, ind, lsmhc
    real(kind=8) :: abs
    character(len=24) :: lisint
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=8), pointer :: delg(:) => null()
!
    call jeveuo(lisint, 'L', lintf)
!-- RECUPERATION DE LA MATRICE DE RAIDEUR
    call jeveuo(jexnum(imped(1:19)//'.VALM', 1), 'E', lklibr)
!
!-- RECUPERATION DES INFOS DU NUME_DDL
    call jeveuo(nume91(1:14)//'.SMOS.SMDI', 'L', vi=smdi)
    call jeveuo(nume91(1:14)//'.SMOS.SMHC', 'L', lsmhc)
    call jeveuo(nume91(1:14)//'.NUME.DELG', 'L', vi=delg)
!
!
    do k1 = 1, nbint
!
        indin1 = '&&VEC_DDL_INTF_'//zk8(lintf+k1-1)
        call jeveuo(indin1, 'L', llint1)
        call jelira(indin1, 'LONMAX', nbddl1)
!
        do m1 = 1, nbddl1
            if (zi(llint1+m1-1) .gt. 0) then
                call ddllag(nume91, zi(llint1+m1-1), nbeq1, lag1, lag2)
!-- SUPRESSION DES COUPLAGES L1 / L2
                if (lag1 .gt. 1) then
                    l1 = smdi(lag1)-smdi(1+lag1-2)-1
                    ind = smdi(1+lag1-2)
                    do n1 = 1, l1
                        zr(lklibr+ind+n1-1) = 0.d0
                    end do
                end if
                if (lag2 .gt. 1) then
                    l1 = smdi(lag2)-smdi(1+lag2-2)-1
                    ind = smdi(1+lag2-2)
                    do n1 = 1, l1
                        zr(lklibr+ind+n1-1) = 0.d0
                    end do
                end if
!
!-- SUPPRESSION DES COUPLAGES EQ / L1
                if (zi(llint1+m1-1) .gt. 1) then
                    l1 = smdi(1+zi(llint1+m1-1)-1)-smdi(1+zi( &
                                                        llint1+m1-1)-2)-1
                    ind = smdi(1+zi(llint1+m1-1)-2)
                    do j1 = 1, l1
!-- ON TESTE DANS LE NUME.DELG SI LA VALEUR EST NEGATIVE
                        if (delg(1+zi4(lsmhc+ind+j1-1)-1) .lt. 0) then
                            zr(lklibr+ind+j1-1) = 0.d0
                        end if
                    end do
                end if
!-- ON REND LA DIAGONALE POSITIVE
                zr(lklibr+smdi(lag1)-1) = abs(zr(lklibr+smdi(1+ &
                                                             lag1-1)-1))
                zr(lklibr+smdi(lag2)-1) = abs(zr(lklibr+smdi(1+ &
                                                             lag2-1)-1))
            end if
!
        end do
!
    end do
!
end subroutine
