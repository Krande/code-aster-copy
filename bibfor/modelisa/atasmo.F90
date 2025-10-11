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

subroutine atasmo(neq, az, apddl, apptr, numedz, &
                  ataz, basez, nblia, nmul, numatz)
    implicit none
!
!     ATASMO  --  LE BUT DE CETTE ROUTINE EST DE CREER LA MATR_ASSE
!                 DE NOM ATA.
!                 LE .VALM DE CETTE MATR_ASSE VA CONTENIR LES TERMES
!                 DU PRODUIT DE MATRICES AT*A.
!                 A EST 'CONCEPTUELLEMENT' UNE MATRICE RECTANGLE
!                 DE 'HAUTEUR' NBLIG ET DE LARGEUR NEQ.
!                 A EST 'INFORMATIQUEMENT' UNE COLLECTION NUMEROTEE
!                 COMPORTANT NBLIG OBJETS QUI SONT DES VECTEURS DE REELS
!                 DE LONGUEUR NEQ.
!                 CHACUN DE CES VECTEURS EST UNE LIGNE DE LA MATRICE A.
!                 LA MATR_ASSE ATA VA DONC ETRE SYMETRIQUE ET A
!                 VALMURS REELLES. D'AUTRE PART ON VA LUI AFFECTER
!                 UN PROFIL MORSE.
!
!   ARGUMENT        E/S  TYPE         ROLE
!    AZ              IN    K*     NOM DE LA COLLECTION DES VECTEURS
!                                 LIGNES (I.E. AZ EST LA MATRICE
!                                 RECTANGULAIRE POUR LAQUELLE ON VA
!                                 CALCULER LE PRODUIT AZ_T*AZ).
!    NUMEDZ         IN    K*      NOM DU NUME_DDL DECRIVANT LES
!                                 LIGNES DE LA MATRICE AZ
!    BASEZ           IN    K*     NOM DE LA BASE SUR LAQUELLE ON
!                                 CREE LA MATR_ASSE.
!    ATAZ           OUT    K*     NOM DE LA MATR_ASSE SYMETRIQUE
!                                 A VALMURS REELLES DONT LE .VALM
!                                 CONTIENT LE PRODUIT AT*A.
!                                 LE PROFIL DE CETTE MATRICE EST
!                                 EN LIGNE DE CIEL.
!                                 CE PROFIL EST DETERMINE DANS LA
!                                 ROUTINE.
!    NUMATZ         OUT    K*     NOM DU NUME_DDL A CREER POUR ATAZ
!                                 ON LE DETRUIT S'Il EXISTE DEJA
!.========================= DEBUT DES DECLARATIONS ====================
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedup1.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/trir_i4.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterc/vector_vector_add_values.h"
#include "asterc/vector_vector_delete.h"
#include "asterc/vector_vector_get_data.h"
#include "asterc/vector_vector_get_data_sizes.h"
#include "asterc/vector_vector_new.h"
#include "asterc/vector_vector_resize.h"
#include "asterc/vector_vector_sort_unique.h"
!
    character(len=8) :: ma
! -----  ARGUMENTS
    character(len=*) :: az, numedz, ataz, basez, numatz
    integer(kind=8) :: neq, nblia, nmul
    integer(kind=8) :: apddl(*), apptr(*)
! -----  VARIABLES LOCALES
    integer(kind=8) :: j, k, iimax, jhbid, idsuiv, dimacv, jconl, nblig, nblig2
    integer(kind=8) ::  ilig, idligm, nddlt, jacv, jaci, iilib, idlm, iddl
    integer(kind=8) :: ieq, ncoef, jvalm, decal, jrefa, vec_ptr, vecSize, position
    integer(kind=8) :: i, jsmde, ii1, ii2, iii, ii, jj, jdecal, nddltm, kdeb
    character(len=1) :: base
    character(len=14) :: numddl, numedd
    character(len=19) :: ata
    character(len=24) :: ksmhc, ksmdi, krefa, kconl, kvalm
    character(len=24) :: a, ksuiv, khbid
    real(kind=8) :: un, vi, vj, vij
    integer(kind=8), pointer :: acompac_1er(:) => null()
    integer(kind=8), pointer :: acompac_nbt(:) => null()
    integer(kind=8), pointer :: v_smdi(:) => null()
    integer(kind=4), pointer :: v_smhc(:) => null()
!
! ========================= DEBUT DU CODE EXECUTABLE ==================
    call jemarq()
!
!
! --- 1. INITIALISATIONS :
!     ---------------------
    base = basez
    a = az
    ata = ataz
    numedd = numedz
    nblig = nblia*nmul
    numddl = numatz
    call detrsd('MATR_ASSE', ata)
    call detrsd('NUME_DDL', numddl)
    call copisd('NUME_EQUA', base, numedd//'.NUME', numddl//'.NUME')
!     UNE DROLE DE GLUTE A RESORBER :
    call jedup1(numedd//'.MLTF.RENU', base, numddl//'.MLTF.RENU')
!
    un = 1.0d0
!
!
    call jelira(a, 'NMAXOC', nblig2)
    ASSERT(nblig .gt. 0)
    ASSERT(nblig .le. nblig2)
!
    ksmdi = numddl//'.SMOS.SMDI'
!
!
!     IIMAX   LONGUEUR MAXIMUM ADMISSIBLE DE KHBID ET ISUIV
!             DANS LA ROUTINE MOINSR IIMAX EST AUGMENTE SI NECESS.
!     KHBID   TABLE DES NUMEROS DE LIGNE
!     ISUIV   TABLE DES CHAINAGES DE LA STRUCTURE CHAINEE
!            (SMDI,SMHC,ISUIV) QUI EST CONTRUITE AVANT D'OBTENIR LA
!               STRUCTURE COMPACTE (SMDI,SMHC) DE LA MATRICE .
    khbid = '&&ATASMO.SMOS.SMHC'
    ksuiv = '&&ATASMO.ANCIEN.ISUIV'
!     ON COMMENCE AVEC IIMAX=100
    iimax = 100
    call wkvect(khbid, 'V V I', iimax, jhbid)
    call wkvect(ksuiv, 'V V I', iimax, idsuiv)
!
!
!     2. CALCUL DE DIMACV
!        DIMACV= NOMBRE TOTAL DE TERMES NON NULS DANS A
!        -- ALLOCATION DE &&ATASMO.ACOMPAC_NBT
!        -- ALLOCATION DE &&ATASMO.ACOMPAC_1ER
!     ----------------------------------------------------------------
    dimacv = 0
    AS_ALLOCATE(vi=acompac_nbt, size=nblig)
    AS_ALLOCATE(vi=acompac_1er, size=nblig)
    do j = 1, nmul
        do ilig = 1, nblia
            call jeveuo(jexnum(a, ilig+(j-1)*nblia), 'L', idligm)
            nddltm = apptr(ilig+1)-apptr(ilig)
            nddlt = 0
            do i = 1, nddltm
                if (zr(idligm-1+i) .ne. 0.d0) nddlt = nddlt+1
            end do
!           ASSERT(NDDLT.GT.0)
            acompac_nbt(ilig+(j-1)*nblia) = nddlt
            acompac_1er(ilig+(j-1)*nblia) = dimacv+1
            dimacv = dimacv+nddlt
            call jelibe(jexnum(a, ilig+(j-1)*nblia))
        end do
    end do
!       ASSERT(DIMACV.GT.0)
    dimacv = max(dimacv, 1)
!
!
!     3. COMPACTAGE DE A : ON NE CONSERVE QUE LES TERMES /= 0 AINSI
!        QUE LES INDICES CORRESPONDANT :
!     ----------------------------------------------------------------
    call wkvect('&&ATASMO.ACOMPAC_I', 'V V S', dimacv, jaci)
    call wkvect('&&ATASMO.ACOMPAC_V', 'V V R', dimacv, jacv)
    k = 0
    do j = 1, nmul
        do ilig = 1, nblia
            call jeveuo(jexnum(a, ilig+(j-1)*nblia), 'L', idligm)
            nddlt = apptr(ilig+1)-apptr(ilig)
            jdecal = apptr(ilig)
            kdeb = k
            do i = 1, nddlt
                if (zr(idligm-1+i) .ne. 0.d0) then
                    k = k+1
                    zi4(jaci-1+k) = apddl(jdecal+i)
                    zr(jacv-1+k) = zr(idligm-1+i)
                end if
            end do
!         ON DOIT TRIER LE TABLEAU DES NUMEROS D'EQUATIONS
!         CAR MOINSR S'ATTEND A UN TABLEAU ORDONNE
!         LE TABLEAU DES VALEURS EST AUSSI PERMUTE POUR
!         RESPECTER CE TRI
            call trir_i4(zi4(jaci+kdeb), zr(jacv+kdeb), 1, k-kdeb)
            call jelibe(jexnum(a, ilig+(j-1)*nblia))
        end do
    end do
!
!
!     4. OBJETS DU NUME_DDL : .SMHC ET .SMDI :
!     ----------------------------------------
!     IILIB  : 1-ERE PLACE LIBRE
    iilib = 1
    call wkvect('&&ATASMO.LMBID', 'V V S', 1, idlm)
!
    call vector_vector_new(vec_ptr)
    call vector_vector_resize(vec_ptr, neq)
!
!     4.1 : ON FORCE LA PRESENCE DES TERMES DIAGONAUX:
    do ieq = 1, neq
        zi4(idlm) = ieq
        call vector_vector_add_values(vec_ptr, ieq-1, 1, zi4(idlm))
    end do
!
!
!     4.2 : ON INSERE LES VRAIS TERMES :
    do ilig = 1, nblig
!       NDDLT : NOMBRE DE TERMES NON NULS POUR ILIG
        nddlt = acompac_nbt(ilig)
        idlm = jaci-1+acompac_1er(ilig)
!
!       -- INSERTION DES COLONNES DE L'ELEMENT DANS
!           LA STRUCTURE CHAINEE
        do iddl = 0, nddlt-1
            position = zi4(idlm+iddl)-1
            call vector_vector_add_values(vec_ptr, position, iddl+1, zi4(idlm))
        end do
    end do
!
!
!     DESIMBRIQUATION DE CHAINES POUR OBTENIR LA STRUCTURE COMPACTE
!     (ZI(JSMDI),SMHC) DE LA MATRICE
    ksmhc = numddl//'.SMOS.SMHC'
    call vector_vector_sort_unique(vec_ptr)
    call vector_vector_get_data_sizes(vec_ptr, vecSize, ncoef)
    call wkvect(ksmdi, base//' V I', neq, vi=v_smdi)
    call wkvect(ksmhc, base//' V S', ncoef, vi4=v_smhc)
    call vector_vector_get_data(vec_ptr, v_smdi, v_smhc)
    call vector_vector_delete(vec_ptr)
    v_smdi(1) = 1
!   Le contenu du SMDI ne correspond pas à la position du premier terme
!   d'une ligne donnée dans le SMHC mais à la position du terme diagonal
!   de la ligne dans le SMHC, c'est pourquoi on doit réaliser la boucle
!   suivante
    do iddl = 2, neq-1
        v_smdi(iddl) = v_smdi(iddl+1)-1
    end do
    v_smdi(neq) = ncoef
!
!
!     5. OBJET DU NUME_DDL :  .SMDE :
!     ----------------------------------------------------
    call wkvect(numddl//'.SMOS.SMDE', base//' V I', 6, jsmde)
    zi(jsmde+1-1) = neq
    zi(jsmde+2-1) = ncoef
    zi(jsmde+3-1) = 1
!
!
!
!     6. OBJETS: MATR_ASSE.REFA ET MATR_ASSE.CONL:
!     --------------------------------------------------
    krefa = ata//'.REFA'
    call wkvect(krefa, base//' V K24', 20, jrefa)
    call dismoi('NOM_MAILLA', numedd, 'NUME_DDL', repk=ma)
    zk24(jrefa-1+1) = ma
    zk24(jrefa-1+2) = numddl
    zk24(jrefa-1+9) = 'MS'
    zk24(jrefa-1+10) = 'NOEU'
    zk24(jrefa-1+11) = 'MPI_COMPLET'
    kconl = ata//'.CONL'
    call wkvect(kconl, base//' V R', neq, jconl)
    do i = 1, neq
        zr(jconl+i-1) = un
    end do
!
!
!     7. OBJET: MATR_ASSE.VALM :
!     --------------------------------------------------
    kvalm = ata//'.VALM'
    call jecrec(kvalm, base//' V R', 'NU', 'DISPERSE', 'CONSTANT', &
                1)
    call jeecra(kvalm, 'LONMAX', ncoef)
    call jecroc(jexnum(kvalm, 1))
    call jeveuo(jexnum(kvalm, 1), 'E', jvalm)
    do ilig = 1, nblig
!       NDDLT : NOMBRE DE TERMES NON NULS POUR ILIG
        nddlt = acompac_nbt(ilig)
        decal = acompac_1er(ilig)
!
!       -- CALCUL DE .VALM(II,JJ) :
        do j = 1, nddlt
            vj = zr(jacv-1+decal-1+j)
            jj = zi4(jaci-1+decal-1+j)
            ASSERT(jj .le. neq)
            ii2 = v_smdi(jj)
            if (jj .eq. 1) then
                ii1 = 1
            else
                ii1 = v_smdi(jj-1)+1
            end if
            ASSERT(ii2 .ge. ii1)
            do i = 1, j
                vi = zr(jacv-1+decal-1+i)
                ii = zi4(jaci-1+decal-1+i)
                vij = vi*vj
!           -- CUMUL DE VIJ DANS .VALM :
                do iii = ii1, ii2
                    if (v_smhc(iii) .eq. ii) goto 110
                end do
                ASSERT(.false.)
110             continue
                zr(jvalm-1+iii) = zr(jvalm-1+iii)+vij
            end do
        end do
    end do
!
!
!     9. MENAGE :
!     ------------
    call jedetr('&&ATASMO.SMOS.SMHC')
    call jedetr('&&ATASMO.ANCIEN.ISUIV')
    AS_DEALLOCATE(vi=acompac_nbt)
    AS_DEALLOCATE(vi=acompac_1er)
    call jedetr('&&ATASMO.ACOMPAC_V')
    call jedetr('&&ATASMO.ACOMPAC_I')
    call jedetr('&&ATASMO.LMBID')
!
    call jedema()
end subroutine
