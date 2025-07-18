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

subroutine prjlis(moda, maa, modb, mab, nbnoa, &
                  nbnob, motcle, linta, lintb, intfa, &
                  intfb, fpliao, fplibo, iada, iadb, &
                  numlis, matprj, modgen, ssta, sstb)
    implicit none
!  O. NICOLAS     DATE 01/08/04
!-----------------------------------------------------------------------
!  BUT : < CALCUL DE LA MATRICE DE PROJECTION >
!
!  CALCULER LA MATRICE DE PROJECTION D'UNE INTERFACE ESCLAVE SUR UNE
!  MAITRE POUR LA SOUS STRUCTURATION DYNAMIQUE
!
!-----------------------------------------------------------------------
!
! MODA /I/ : NOM DU MODELE MAITRE
! MODB /I/ : NOM DU MODELE ESCLAVE
! MAA /I/ : NOM DU MAILLAGE MAITRE
! MAB /I/ : NOM DU MAILLAGE ESCLAVE
! NBNOA /I/ : NOMBRE DE NOEUDS D'INTERFACE DU MAILLAGE MAITRE
! NBNOB /I/ : NOMBRE DE NOEUDS D'INTERFACE DU MAILLAGE ESCLAVE
! MOTCLE /I/ : MOT CLE RENSEIGNANT LE GROUPE OU LA MAILLE MAITRE
! LINTA /I/ : NOM DE L'INTERFACE AMONT MAITRE
! LINTB /I/ : NOM DE L'INTERFACE AMONT ESCLAVE
! INTFA /I/ : NOM DE L'INTERFACE MAITRE
! INTFB /I/ : NOM DE L'INTERFACE ESCLAVE
! MATPRJ /O/ : NOM DE LA MATRICE D'OBSERVATION
! NUMLIS /I/ : NUMERO INTERFACE COURANTE
! FPLIAO /I/ : FAMILLE DES PROFNO MATRICES DE LIAISON ORIENTEES SSTA
! FPLIBO /I/ : FAMILLE DES PROFNO MATRICES DE LIAISON ORIENTEES SSTB
! IADA   /I/ : VECTEUR DES CARACTERISTIQUES LIAISON SSTA
! IADB   /I/ : VECTEUR DES CARACTERISTIQUES LIAISON SSTB
! MODGEN  /I/ : NOM K8 DU MODELE GENERALISE
! SSTA    /I/ : NOM K8 DE LA SOUS-STRUCTURE MAITRE
! SSTB    /I/ : NOM K8 DE LA SOUS-STRUCTURE ESCLAVE
!
!
!
!
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/geolis.h"
#include "asterfort/infniv.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/pj2dco.h"
#include "asterfort/pj3dco.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
!
    character(len=4) :: zcst
    character(len=8) :: linta, lintb, moda, modb, maa, mab, intfa, intfb
    character(len=8) :: matprj, nonob, nonoa, nomg, modgen, ssta, sstb
    character(len=16) :: tymocl(2), motcle(2), motfac, corres
    character(len=24) :: inta, intb, fpliao, fplibo, toto, geoma, geomb
    integer(kind=8) :: ibid, nbnoa, nbnob, llinta, llintb, nbmaa, ndim
    integer(kind=8) ::    nbnob2, idecal, inoa, inob, nbterm
    integer(kind=8) :: itemcm, itemtm, nnoa, nunoa, nunob, nunoa2, nunoai, nunobi, i, j
    integer(kind=8) :: iinob, llplia, llplib, icompa, icompb, iadoa, iadob, nbec
    integer(kind=8) :: numlis, iada(3), iadb(3), nbcmpm, ifm, niv
    parameter(nbcmpm=10)
    integer(kind=8) :: idecoa(nbcmpm), idecob(nbcmpm)
    real(kind=8) :: rbid, beta, coefa
    integer(kind=8), pointer :: pjef_nb(:) => null()
    real(kind=8), pointer :: pjef_cf(:) => null()
    integer(kind=8), pointer :: limanua(:) => null()
    integer(kind=8), pointer :: nldesca(:) => null()
    integer(kind=8), pointer :: nldescb(:) => null()
    integer(kind=8), pointer :: pjef_nu(:) => null()
!
!-----------------------------------------------------------------------
!
    call jemarq()
!
!     -- IMPRESSION DE LA RELATION SI INFO:2:
!     ---------------------------------------
    call infniv(ifm, niv)
!
    motfac = 'LIAISON'
    tymocl(1) = 'MAILLE'
    tymocl(2) = 'GROUP_MA'
!
!--------------LES NOMBRES DES NOEUDS DES INTERFACES
    inta = linta//'.IDC_LINO'
    call jenonu(jexnom(inta(1:13)//'NOMS', intfa), ibid)
    call jelira(jexnum(inta, ibid), 'LONMAX', nbnoa)
!
    intb = lintb//'.IDC_LINO'
    call jenonu(jexnom(intb(1:13)//'NOMS', intfb), ibid)
    call jelira(jexnum(intb, ibid), 'LONMAX', nbnob)
!
!--------------LES LISTES DES NUMEROS DES NOEUDS DES INTERFACES
    call jenonu(jexnom(linta//'.IDC_NOMS', intfa), ibid)
    call jeveuo(jexnum(linta//'.IDC_LINO', ibid), 'L', llinta)
    call jeveuo(linta//'.IDC_DEFO', 'L', vi=nldesca)
!
    call jenonu(jexnom(lintb//'.IDC_NOMS', intfb), ibid)
    call jeveuo(jexnum(lintb//'.IDC_LINO', ibid), 'L', llintb)
    call jeveuo(lintb//'.IDC_DEFO', 'L', vi=nldescb)
!
!
!--------------LE NOMBRE DES MAILLES DE L'INTERFACE MAITRE
    call reliem(moda, maa, 'NU_MAILLE', motfac, 1, &
                2, motcle, tymocl, '&&PRJLIS.LIMANUA', nbmaa)
!--------------LA LISTE DES NUMEROS DES MAILLES DE L'INTERFACE MAITRE
    call jeveuo('&&PRJLIS.LIMANUA', 'L', vi=limanua)
!
!---ON FAIT LA PROJECTION
    ndim = 3
    call dismoi('Z_CST', maa, 'MAILLAGE', repk=zcst)
    if (zcst .eq. 'OUI') ndim = 2
!
    corres = '&&PRJLIS.CORRES'
!
!--------------TRANSFORMATION DE LA GEOMETRIE POUR LA PROJECTION
    geoma = '&&PRJLIS.GEOM_TRANSA'
    geomb = '&&PRJLIS.GEOM_TRANSB'
    call geolis(modgen, ssta, sstb, intfa, intfb, &
                geoma, geomb, '&&PRJLIS.LIMANUA', nbmaa)
!
!--------------CALCUL DE CORRES
!
    if (ndim .eq. 2) then
        call pj2dco('PARTIE', moda, modb, nbmaa, limanua, &
                    nbnob, nldescb, geoma, geomb, corres, &
                    .false._1, rbid, 0.d0)
    else if (ndim .eq. 3) then
        call pj3dco('PARTIE', moda, modb, nbmaa, limanua, &
                    nbnob, nldescb, geoma, geomb, corres, &
                    .false._1, rbid, 0.d0)
    end if
!
    call jeveuo(corres//'.PJEF_NB', 'L', vi=pjef_nb)
    call jeveuo(corres//'.PJEF_NU', 'L', vi=pjef_nu)
    call jeveuo(corres//'.PJEF_CF', 'L', vr=pjef_cf)
    call jelira(corres//'.PJEF_NB', 'LONMAX', nbnob2)
!
!      CALL UTIMSD(6,2,.TRUE.,.TRUE.,'&&PRJLIS.CORRES',1,' ')
!
!
!-------------ON RECUPERE LES COEFFICIENTS DE LA PROJECTION
    toto = 'TUTU'
    call wkvect(toto, 'V V R', nbnob*nbnoa, itemtm)
!
! Initialisation de la matrice d'observation
    do inob = 1, nbnob
        do inoa = 1, nbnoa
            zr(itemtm+(inob-1)*nbnoa+inoa-1) = 0.d0
        end do
    end do
!
! Remplissage et impression de la matrice d'observation
    idecal = 0
    beta = 0.0d0
    nbterm = 0
    iinob = 1
! boucle sur l'ensemble des noeuds du modele maitre
    do inob = 1, nbnob2
! on recupere le nombre de noeuds maitre lie au noeud esclave courant
        nnoa = pjef_nb(inob)
        nbterm = nnoa+1
! si le nbre de noeud maitre lie au noeud esclave courant est > 0
        if (nnoa .gt. 0) then
            nunob = zi(llintb+iinob-1)
            nunobi = nldescb(nunob)
            nonob = int_to_char8(nunobi)
            if (niv .eq. 2) then
                write (ifm, *) ' '
                write (ifm, *) '_RELA IMPRESSION D''UNE RELATION'//&
     &'            LINEAIRE ENTRE '&
     &            , nbterm, ' DDLS. (AVANT NORMALISATION DE LA RELATION)'
                write (ifm, 1001)-1.d0, nonob
            end if
! boucle sur le nombre de noeud maitre lie au noeud esclave courant
            do inoa = 1, nnoa
                nunoa = pjef_nu(1+idecal-1+inoa)
                coefa = pjef_cf(1+idecal-1+inoa)
                nonoa = int_to_char8(nunoa)
! boucle sur le nombre de noeud maitre present dans l'interface
                do j = 1, nbnoa
! si le noeud maitre courant est present dans la liste des noeuds
! maitres d'interface on stocke la valeur du coefficient
                    nunoa2 = zi(llinta+j-1)
                    nunoai = nldesca(nunoa2)
                    if (nunoa .eq. nunoai) then
! On stocke la valeur du coefficient dans la matrice d'observation
! le stockage est donc C(Nbre Noeud esclave,Nbre Noeud maitre)
! l'ordre est donc celui de l'interface esclave pour
! les lignes de la matrice et celui de l'interface maitre pour les
! colonnes
                        zr(itemtm+(iinob-1)*nbnoa+j-1) = coefa
                        if (niv .eq. 2) then
                            write (ifm, 1001) coefa, nonoa
                        end if
                    end if
                end do
            end do
            if (niv .eq. 2) then
                write (ifm, *) '_RELA = ', beta
            end if
            idecal = idecal+nnoa
            iinob = iinob+1
        end if
    end do
!
! ************************************************************
! Recuperation des donnees par composantes
    nomg = 'DEPL_R'
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    if (nbec .gt. 10) then
        call utmess('F', 'MODELISA_94')
    end if
!
    call jeveuo(jexnum(fpliao, numlis), 'L', llplia)
    call jeveuo(jexnum(fplibo, numlis), 'L', llplib)
!
    call wkvect(matprj, 'G V R', iada(1)*iadb(1), itemcm)
! Initialisation de la matrice d'observation
    do inob = 1, iadb(1)
        do inoa = 1, iada(1)
            zr(itemcm+(inob-1)*iada(1)+inoa-1) = 0.d0
        end do
    end do
!
    do inob = 1, nbnob
        iadob = zi(llplib+(inob-1)*(1+nbec))
        call isdeco(zi(llplib+(inob-1)*(1+nbec)+1), idecob, nbcmpm)
        icompb = iadob-1
        do i = 1, nbcmpm
            if (idecob(i) .gt. 0) then
                icompb = icompb+1
                do inoa = 1, nbnoa
                    iadoa = zi(llplia+(inoa-1)*(1+nbec))
                    call isdeco(zi(llplia+(inoa-1)*(1+nbec)+1), idecoa, nbcmpm)
                    icompa = iadoa-1
                    do j = 1, nbcmpm
                        if ((idecoa(j) .gt. 0) .and. (i .eq. j)) then
! On se limite au repere globaux
                            icompa = icompa+i
                            zr(itemcm+(icompb-1)*iada(1)+icompa-1) = &
                                zr(itemtm+(inob-1)*nbnoa+inoa-1)
                        end if
                    end do
                end do
            end if
        end do
    end do
!
!
! ************************************************************
!
!-------------FORMAT D'IMPRESSION
1001 format(' _RELA ', e14.7, a10, a10)
!
    call jedetr(toto)
    call jedetr(geoma)
    call jedetr(geomb)
    call jedetr(corres)
    call jedetr('&&PRJLIS')
!
    call jedema()
end subroutine
