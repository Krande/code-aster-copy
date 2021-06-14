# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------


import argparse
import glob
import logging
import os
import re
import sys

import medcoupling as mc

logger = logging.getLogger()

class ColoredFormatter(logging.Formatter):
    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(30,38)
    COLORS = { 'WARNING': YELLOW, 'INFO': WHITE, 'DEBUG': BLUE, 'CRITICAL': YELLOW, 'ERROR': RED }
    def __init__(self, *args, **kwargs):
        logging.Formatter.__init__(self, *args, **kwargs)
    def format(self, record):
        RESET_SEQ = "\033[0m"
        COLOR_SEQ = "\033[1;%dm"
        record.levelname = COLOR_SEQ % ColoredFormatter.COLORS[record.levelname] + record.levelname + RESET_SEQ
        return logging.Formatter.format(self, record)


def setVerbose(verbose=1, code_aster=False):
    global logger

    if not code_aster:
        logger = logging.getLogger()
        formatter = ColoredFormatter('%(levelname)s : %(asctime)s : %(message)s ',style='%')
        formatter.default_time_format = '%H:%M:%S'
        formatter.default_msec_format = "%s.%03d"
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    verbose_map = { 0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}
    logger.setLevel(verbose_map[verbose])

def operate(isTopCall,effective_dest_file_name,nb_dest_par,origin_seq):
    pat = re.compile("([\d]+)([\s\S]+)$")

    check_remain = set()

    dest_fns = []

    for fn in effective_dest_file_name:
        bn = os.path.basename(fn)
        base,ext = os.path.splitext(bn)
        if ext not in [".med",".rmed"]:
            logger.error("file {} in pattern does not look MED file".format(fn))
            sys.exit(1)
        base_rev = base[::-1]
        m = pat.match(base_rev)
        if not m:
            logger.warn("Le fichier MED {} du pattern de destination : impossible d extraire le proc Id !".format(fn))
            sys.exit(1)
        proc_id = int(m.group(1)[::-1])
        remain_fn = m.group(2)[::-1]
        check_remain.add( remain_fn )
        dest_fns.append((fn,proc_id))

    if len(check_remain) != 1:
        logger.error("Il y a un probleme dans la detection des procIds")
        sys.exit(1)

    if set([nb for _,nb in dest_fns]) != set(range(nb_dest_par)):
        logger.error("Il y a un probleme dans la detection des procIds (sur les entiers dans les noms de fichiers)")
        sys.exit(1)

    logger.info("Lecture du fichier MED original {}".format(origin_seq))
    mm_orginal = mc.MEDFileMesh.New(origin_seq)
    meshDim = mm_orginal.getMeshDimension()
    orig_pt1 = mc.MEDCoupling1SGTUMesh(mm_orginal[-meshDim]).getNodalConnectivity() # -mm_orginal.getMeshDimension() pour dire que l'on prend les MED_POINT1
    # orig_pt1 contient les node ids qui doivent porter un MED_POINT
    logger.info("On compte {} MED_POINT1 dans l'orginal".format(len(orig_pt1)))
    meshDim = mm_orginal.getMeshDimension()
    num_orig_pt1 = None
    fam_orig_pt1 = None
    orig_pt1 = mc.MEDCoupling1SGTUMesh(mm_orginal[-meshDim]).getNodalConnectivity()
    if mm_orginal.getNumberFieldAtLevel(-meshDim):
        num_orig_pt1 = mm_orginal.getNumberFieldAtLevel(-meshDim)

    if mm_orginal.getFamilyFieldAtLevel(-meshDim):
        fam_orig_pt1 = mm_orginal.getFamilyFieldAtLevel(-meshDim)

    check_all_pts1 = []

    for dest_fn,dest_proc in dest_fns:
        dest_mm = mc.MEDFileMesh.New(dest_fn)
        pt1_in_current_dest = dest_mm.getGlobalNumFieldAtLevel(1).buildIntersection(orig_pt1)
        check_all_pts1.append(pt1_in_current_dest)
        pt1_to_be_present = dest_mm.getGlobalNumFieldAtLevel(1).findIdForEach(pt1_in_current_dest)

        pt1_cur_mesh = None
        pt1_already_present = mc.DataArrayInt([])
        pt1_fam_before = None
        pt1_num_before = None
        if -meshDim in dest_mm.getNonEmptyLevels():
            pt1_cur_mesh = mc.MEDCoupling1SGTUMesh(dest_mm[-meshDim])
            pt1_fam_before = dest_mm.getFamilyFieldAtLevel(-meshDim)
            pt1_already_present = pt1_cur_mesh.getNodalConnectivity()
            if not pt1_to_be_present.buildIntersection(pt1_already_present).isEqual(pt1_already_present):
                raise RuntimeError("Les POINT1 dans {} ne sont pas exactement dans ceux attendus de l'origine !".format(fn))
            if dest_mm.getNumberFieldAtLevel(-meshDim):
                pt1_num_before = dest_mm.getNumberFieldAtLevel(-meshDim)

        pt1_to_add = pt1_to_be_present.buildSubstraction(pt1_already_present)
        if pt1_to_add.empty():
            logger.info("Pour le proc {} il y a {} MED_POINT1 et rien à faire !".format(dest_proc,len(pt1_already_present)))
            continue

        pt1_entry_to_be_present = orig_pt1.findIdForEach(dest_mm.getGlobalNumFieldAtLevel(1)[pt1_to_add])
        logger.info("Pour le proc {} on rajoute {} MED_POINT1 (avant {} apres {})".format(dest_proc,len(pt1_to_add),len(pt1_already_present),len(pt1_to_be_present)))
        # on rajoute le mesh
        pt1_to_add_mesh = mc.MEDCoupling1SGTUMesh(dest_mm.getName(),mc.NORM_POINT1)
        pt1_to_add_mesh.setCoords(dest_mm.getCoords())
        pt1_to_add_mesh.setNodalConnectivity(pt1_to_add)
        pt1_after = mc.MEDCoupling1SGTUMesh.Merge1SGTUMeshesOnSameCoords([elt for elt in [pt1_cur_mesh,pt1_to_add_mesh] if elt is not None])
        pt1_after.setName(dest_mm.getName())
        pt1_fam = mc.DataArrayInt.Aggregate([elt for elt in [pt1_fam_before,fam_orig_pt1[pt1_entry_to_be_present]] if elt is not None])
        pt1_num = None
        if pt1_num_before:
            pt1_num = mc.DataArrayInt.Aggregate([pt1_num_before,num_orig_pt1[pt1_entry_to_be_present]])

        dest_mm[-meshDim] = pt1_after
        if fam_orig_pt1:
            dest_mm.setFamilyFieldArr(-meshDim,pt1_fam)

        if pt1_num:
            dest_mm.setRenumFieldArr(-meshDim,pt1_num)
        logger.warn("On reecrit le proc {} avec les {} MED_POINT1 dans {}".format(dest_proc,len(pt1_to_be_present),dest_fn))
        dest_mm.write(dest_fn,2)
    check = mc.DataArrayInt.Aggregate(check_all_pts1)
    check.sort()
    check = check.buildUnique()
    isok = check.isEqual(orig_pt1)
    if isok:
        return True
    if not isTopCall:
        return False
    pt1_simply_ignored = orig_pt1.buildSubstraction(check)
    # par defaut on prend le premier fichier pour lui rajouter les points manquants
    filename_with_pts_to_add = dest_fns[0][0]
    logger.warn("Des noeuds (au nombre de {}) ont été oubliés ainsi que les cellules MED_POINT1 dessus ! On va les rajouter dans le fichier {}".format(len(pt1_simply_ignored),filename_with_pts_to_add))
    first_part_mm = mc.MEDFileMesh.New(filename_with_pts_to_add)
    coo = first_part_mm.getCoords()
    node_num = first_part_mm.getNumberFieldAtLevel(1)
    node_fam = first_part_mm.getFamilyFieldAtLevel(1)
    node_glob = first_part_mm.getGlobalNumFieldAtLevel(1)
    coo = mc.DataArrayDouble.Aggregate([coo,mm_orginal.getCoords()[pt1_simply_ignored]])
    if node_num:
        node_num = mc.DataArrayInt.Aggregate([node_num,mm_orginal.getNumberFieldAtLevel(1)[pt1_simply_ignored]])
    if node_fam:
        node_fam = mc.DataArrayInt.Aggregate([node_fam,mm_orginal.getFamilyFieldAtLevel(1)[pt1_simply_ignored]])
    if node_glob:
        node_glob = mc.DataArrayInt.Aggregate([node_glob,pt1_simply_ignored])
    first_part_mm.setCoords(coo)
    first_part_mm.setFamilyFieldArr(1,node_fam)
    first_part_mm.setRenumFieldArr(1,node_num)
    first_part_mm.setGlobalNumFieldAtLevel(1,node_glob)
    first_part_mm.write(filename_with_pts_to_add,2)
    logger.warn("Les {} noeuds ont ete rajouter dans ! On va lui rajouter les MED_POINT1 dessus".format(len(pt1_simply_ignored),filename_with_pts_to_add))
    return operate(False,effective_dest_file_name,nb_dest_par,origin_seq)

def add_pt1(args):
    setVerbose(args["verbosity"])

    if args["nb_dest_par"] < 2:
        logger.warn("Le nb de proc de destination {} est < 2 : Rien à faire".format(args["nb_dest_par"]))
        sys.exit(0)

    logger.info("Le fichier MED original depuis lequel on prend les MED_POINT1 : \"{}\"".format(args["origin_seq"]))

    effective_dest_file_name = glob.glob(args["dest_par"])

    if len(effective_dest_file_name) != args["nb_dest_par"]:
        logger.error("Le nombre de fichiers MED dans le pattern {} ({}) est different de celui attendu {}".format(args["dest_par"],len(effective_dest_file_name),args["nb_dest_par"]))
        sys.exit(1)

    logger.info("Les fichiers MED destinations dans lesquels on va rajouter les MED_POINT1 : \"{}\"".format(effective_dest_file_name))

    if not args["force"]:
        res = input("Confirmez vous l'ecriture dans les fichiers {}. Taper 0 pour arreter".format(effective_dest_file_name))
        if res == "0":
            logger.warn("Operation arretée. Les fichiers destinations ne seront pas modifiés")
            sys.exit(1)

    isok = operate(True,effective_dest_file_name,args["nb_dest_par"],args["origin_seq"])

    if not isok:
        raise RuntimeError("Il y a des MED_POINT1 manquant :(")
    logger.info("Tout semble OK.")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-orig','--origin_seq', dest="origin_seq", required=True, help="origin sequential med file path")
    parser.add_argument('-dest','--dest_par', dest="dest_par", required=True, help="destination parallel med file pattern")
    parser.add_argument('-nb_dest','--nb_dest_par', dest="nb_dest_par", type=int, required=True, help="nb of proc destination parallel med file pattern. Just for safety")
    parser.add_argument('-v','--verbose', dest = "verbosity", type=int , default=1, help="verbosity. Default 1 (print info). 0 only errors reported. 2 and higher all debug messages")
    parser.add_argument('-f','--force', dest = "force", type=int , default=0, help="No question, file in destination parallel will be overwritten.")

    args = parser.parse_args()

    setVerbose(args.verbosity)

    if args.nb_dest_par < 2:
        logger.warn("Le nb de proc de destination {} est < 2 : Rien à faire".format(args.nb_dest_par))
        sys.exit(0)

    logger.info("Le fichier MED original depuis lequel on prend les MED_POINT1 : \"{}\"".format(args.origin_seq))

    effective_dest_file_name = glob.glob(args.dest_par)

    if len(effective_dest_file_name) != args.nb_dest_par:
        logger.error("Le nombre de fichiers MED dans le pattern {} ({}) est different de celui attendu {}".format(args.dest_par,len(effective_dest_file_name),args.nb_dest_par))
        sys.exit(1)

    logger.info("Les fichiers MED destinations dans lesquels on va rajouter les MED_POINT1 : \"{}\"".format(effective_dest_file_name))

    if not args.force:
        res = input("Confirmez vous l'ecriture dans les fichiers {}. Taper 0 pour arreter".format(effective_dest_file_name))
        if res == "0":
            logger.warn("Operation arretée. Les fichiers destinations ne seront pas modifiés")
            sys.exit(1)

    isok = operate(True,effective_dest_file_name,args.nb_dest_par,args.origin_seq)

    if not isok:
        raise RuntimeError("Il y a des MED_POINT1 manquant :(")
    logger.info("Tout semble OK.")
