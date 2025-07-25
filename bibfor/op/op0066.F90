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

subroutine op0066()
!
    implicit none
!
! ----------------------------------------------------------------------
!
!  OPERATEUR CALC_ERC_DYN
!
!  CALCUL DE L'ERREUR EN RELATION DE COMPORTEMENT EN DYNAMIQUE
!  SOUS UNE FORMULATION FREQUENTIELLE
!---------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/etausr.h"
#include "asterfort/blkobs.h"
#include "asterc/getres.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/materc.h"
#include "asterfort/titre.h"
#include "asterfort/jeveuo.h"
#include "asterfort/crea_nume_erc.h"
#include "asterfort/sigusr.h"
#include "asterfort/uttcpu.h"
#include "asterfort/uttcpr.h"
#include "asterfort/utmess.h"
#include "asterc/r8depi.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdscr.h"
#include "asterfort/asse_matr_erc.h"
#include "asterfort/asse_vect_erc.h"
#include "asterfort/preres.h"
#include "asterfort/resoud.h"
#include "asterfort/utcrre.h"
#include "asterfort/rsadpa.h"
#include "asterfort/refdaj.h"
#include "asterfort/dismoi.h"
#include "asterfort/eval_erc.h"
#include "asterfort/archi_erc.h"
#include "asterfort/jedetr.h"
#include "asterfort/getvid.h"
!
    character(len=1) :: typcst(2)
    character(len=8) :: baseno, result, matmas, matrig, matamo, nommes, numnu, mymat(2), nomo
    character(len=8) :: k8bid
    character(len=14) :: nom_nume_erc
    character(len=16) :: typcon, nomcmd
    character(len=19) :: dynam1, dynam2, nom_matr_erc, nom_vect_erc, solveu, maprec
    character(len=24) :: matobs(3), lfreqs, matprod(4), mate, carele
    real(kind=8) ::  tps1(4), rtab(2), deuxpi, freq, omega, gamma, alpha, coef_gamma, m_omega2
    real(kind=8) ::  mycoef(2), cout_fon, cout_uv
    complex(kind=8) :: cbid
    integer(kind=8) :: obsdim(3), nbfreq, ifreq, lfreq, iret, ibid, ivecterc, ii, ladpa, valei(8)
    aster_logical :: amor, eval
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call titre()
!
    baseno = '&&OP0066'
    deuxpi = r8depi()
    maprec = '&&OP0066.MAPREC'
!
    call getres(result, typcon, nomcmd)
!
! --- RECUPERATION DES MATRICES ET INFOS POUR LA RESOLUTION DE L'ERC
    call materc(matmas, matrig, matamo, numnu, amor, nommes, &
                lfreqs, nbfreq, matobs, obsdim, gamma, alpha, eval)
!
! --- CREATION DU SOUS-BLOC RELATIF AU PRODUIT H^T*Gr*H
    call blkobs(matobs, obsdim, alpha, matprod)
!
! --- ININTIALISATION DES STRUCTURES DES BLOCS LIES AUX IMPEDANCES
    dynam1 = baseno//'.DYNAM1_MTX'
    dynam2 = baseno//'.DYNAM2_MTX'
    call mtdefs(dynam1, matrig, 'V', ' ')
    call mtdscr(dynam1)
    call mtdefs(dynam2, matrig, 'V', ' ')
    call mtdscr(dynam2)
! --- ALLOCATION DES SD POUR STOCKAGE DES RESULTATS

    call utcrre(result, 2*nbfreq)
    call refdaj('F', result, 2*nbfreq, numnu, 'DYNAMIQUE', [matrig, matmas, ' '], iret)

! --- ASSEMBLAGE INITIAL DE LA MATRICE ERC
!
! --- --- CREATION DU NUME_DDL_GENE ET INITIALISATION DU MATR/VECT_ASSE_GENE

    call crea_nume_erc(baseno, numnu, matprod, nom_nume_erc, nom_matr_erc, &
                       nom_vect_erc, solveu, valei)

! --- PRINT DES DIMENSIONS DU PROBLEME D'ERC
    valei(1) = nbfreq
    valei(2) = obsdim(2)
    valei(4) = obsdim(1)
    valei(5) = obsdim(3)
    call utmess('I', 'ALGELINE7_22', ni=8, vali=valei)
!

!
! --- BOUCLE SUR LES FREQUENCES ---
!====
    call uttcpu('CPU.OP0066', 'INIT', ' ')
!
    do ifreq = 1, nbfreq
        call uttcpu('CPU.OP0066', 'DEBUT', ' ')

! --- --- CALCUL DES COEFF POUR LES MATRICES
!
        call jeveuo(lfreqs, 'L', lfreq)
        freq = zr(lfreq-1+ifreq)
        omega = deuxpi*freq
        m_omega2 = -omega*omega
        coef_gamma = gamma/(1.0d0-gamma)
        typcst(1) = 'R'
        typcst(2) = 'R'
        mymat(1) = matrig
        mymat(2) = matmas
        mycoef(1) = gamma
        mycoef(2) = -coef_gamma*gamma*m_omega2

! --- --- CALCUL DE LA PREMIERE MATRICE DYNAMIQUE DU BLOC DIAGONAL
        call mtcmbl(2, typcst, mycoef, mymat, dynam1, ' ', ' ', 'ELIM=')
! --- --- CALCUL DE LA DEUXIEME MATRICE DYNAMIQUE DU BLOC HORS-DIAGONAL
        mycoef(1) = -gamma
        mycoef(2) = -gamma*m_omega2
        call mtcmbl(2, typcst, mycoef, mymat, dynam2, ' ', ' ', 'ELIM=')
!
! --- --- REMPLISSAGE DU .VALE DE LA MATRICE ASSEMBLEE DE L'ERC
!
        call asse_matr_erc(nom_matr_erc, nom_nume_erc, dynam1, dynam2, matprod)
!
! --- --- REMPLISSAGE DU .VALE DU SECOND MEMBRE DE L'ERC
!
        call asse_vect_erc(baseno, nom_vect_erc, nommes, matobs, obsdim, alpha, ifreq, omega)

! --- --- FACTORISATION DE LA MATRICE ERC

        call preres(solveu, 'V', iret, maprec, nom_matr_erc, ibid, -9999)
!
! --- --- RESOLUTION DU PROBLEME D'ERC SOUS FORME A*x=b
        call jeveuo(nom_vect_erc//'.VALE', 'E', ivecterc)
        cbid = dcmplx(0.d0, 0.d0)
        call resoud(nom_matr_erc, maprec, solveu, '', 1, '', '', 'v', zr(ivecterc), &
                    [cbid], '', .false._1, 0, iret)
!
! --- --- ON EVALUE LA FONCTION COUT SI DEMANDE
        if (eval) then
            call eval_erc(baseno, dynam1, zr(ivecterc), nommes, matobs, obsdim, ifreq, &
                          omega, alpha, cout_fon, cout_uv)
        end if
! --- --- ARCHIVAGE DES RESULTATS
        call archi_erc(result, ifreq, matmas, obsdim, zr(ivecterc), freq, eval, cout_fon, cout_uv)
! --- --- VERIFICATION SI INTERRUPTION DEMANDEE PAR SIGNAL USR1
!
        if (etausr() .eq. 1) call sigusr()
!
! --- --- MESURE CPU
!
        call uttcpu('CPU.OP0066', 'FIN', ' ')
        call uttcpr('CPU.OP0066', 4, tps1)
        if (tps1(4) .gt. .90d0*tps1(1) .and. ifreq .ne. nbfreq) then
            rtab(1) = tps1(4)
            rtab(2) = tps1(1)
            call utmess('Z', 'DYNAMIQUE_13', si=ifreq, nr=2, valr=rtab, &
                        num_except=ASTER_TIMELIMIT_ERROR)
        end if
!
    end do
    ! boucle sur les frequences
!
!     NETOYAGE DES OBJETS JEVEUX TEMPORAIRES ET STOCKAGES FINAUX
!
!
    call dismoi('NOM_MODELE', matrig, 'MATR_ASSE', repk=nomo)
    call getvid(' ', 'CHAM_MATER', scal=mate, nbret=iret)
    if (iret .eq. 0) mate = '        '
    call getvid(' ', 'CARA_ELEM', scal=carele, nbret=iret)
    if (iret .eq. 0) carele = '        '
!
    do ii = 1, 2*nbfreq
        call rsadpa(result, 'E', 1, 'MODELE', ii, &
                    0, sjv=ladpa, styp=k8bid)
        zk8(ladpa) = nomo
        call rsadpa(result, 'E', 1, 'CHAMPMAT', ii, &
                    0, sjv=ladpa, styp=k8bid)
        zk8(ladpa) = mate(1:8)
        call rsadpa(result, 'E', 1, 'CARAELEM', ii, &
                    0, sjv=ladpa, styp=k8bid)
        zk8(ladpa) = carele(1:8)
    end do

    call jedetr(baseno//'.COUNT.ERC.COL')
!
    call jedema()
end subroutine
