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
# person_in_charge: nicolas.tardieu@edf.fr

"""
:py:mod:`mesh_builder` --- Simple mesh builder
**********************************************

The :py:mod:`mesh_builder` helps generate simple meshes (disk, square,
cylinder, cube). The user can indicate the geometrical sizes of the shape
and the mesh fineness using the nrefine parameter.

The functions are exposed as class-method in the :py:class:`Mesh` class.

"""


from ..Commands import MACR_ADAP_MAIL
from ..MacroCommands.macr_adap_mail_ops import HOMARD_INFOS
from ..Objects import Mesh
from ..Supervis import CO
import tempfile
from math import sqrt, cos, sin, pi


def buildSquare(cls, lx=1, ly=1, nrefine=3):
    """Build the quadrilateral mesh of a square.
    Arguments:
        lx [float] : length of the square along the x axis (default 1.).
        ly [float] : length of the square along the y axis (default 1.).
        nrefine [int] : number of mesh refinement iterations (default 3).
    """
    # Mesh creation
    mesh = cls()
    with tempfile.NamedTemporaryFile(mode="w+") as f:
        txtMesh = ("""TITRE\n SQUARE\nFINSF\nCOOR_3D\nN1 0. {1:f} 0.\nN2 {0:f} {1:f} """
        """0.\nN3 0. 0. 0.\nN4 {0:f} 0. 0.\nFINSF\nQUAD4\nSURFACE N2 N4 N3 N1\nFINSF\n"""
        """SEG2\nS13 N1 N3\nS34 N3 N4\nS42 N4 N2\nS21 N2 N1\nFINSF\nGROUP_MA\nSURFACE SURFACE"""
        """\nFINSF\nGROUP_MA\nRIGHT S42\nFINSF\nGROUP_MA\nLEFT S13\nFINSF\nGROUP_MA\nTOP S21\n"""
        """FINSF\nGROUP_MA\nBOTTOM S34\nFINSF\nGROUP_NO\nN1 N1\nFINSF\nGROUP_NO\nN2 N2\nFINSF"""
        """\nGROUP_NO\nN3 N3\nFINSF\nGROUP_NO\nN4 N4\nFINSF\nFIN\n""")
        f.file.write(txtMesh.format(lx, ly))
        f.file.seek(0)  # file ready to be read
        mesh.readAsterFile(f.name)

    # Mesh refinement
    newMesh = mesh
    for _ in range(nrefine):
        HOMARD_INFOS.new()
        resu = MACR_ADAP_MAIL(__use_namedtuple__=True,
                                ADAPTATION='RAFFINEMENT_UNIFORME',
                                MAILLAGE_N=newMesh,
                                MAILLAGE_NP1=CO('newMesh'))
        HOMARD_INFOS.pop()
        newMesh = resu.newMesh
    return newMesh


def buildDisk(cls, radius=1, nrefine=3):
    """Build the quadrilateral mesh of a disk.
    Arguments:
        radius [float] : radius of the disk (default 1).
        nrefine [int] : number of mesh refinement iterations (default 3).
    """
    # Mesh creation
    mesh = cls()
    with tempfile.NamedTemporaryFile(mode="w+") as f:
        txtMesh = ("""TITRE\nDISK\nFINSF\nCOOR_3D\nN1 0 -{1:f} 0\nN4 -{1:f} 0 0\n"""
        """N5 {1:f} 0 0\nN8 0 {1:f} 0\nN10 {0:f} 0 0\nN11 0 -{0:f} 0\nN14 0 {0:f} 0\n"""
        """N15 -{0:f} 0 0\nN20 -{2:f} -{2:f} 0\nN24 {2:f} {2:f} 0\nN25 {2:f} -{2:f} 0\n"""
        """N27 -{2:f} {2:f} 0\nN30 {3:f} -{3:f} 0\nN35 {3:f} {3:f} 0\nN36 -{3:f} -{3:f} 0\n"""
        """N40 -{3:f} {3:f} 0\nN45 0 0 0\nFINSF\nQUAD4\nM73 N1 N20 N45 N25\nM74 N25 N45 N24 N5"""
        """\nM75 N20 N4 N27 N45\nM76 N45 N27 N8 N24\nM87 N10 N30 N25 N5\nM88 N30 N11 N1 N25\n"""
        """M101 N10 N5 N24 N35\nM102 N35 N24 N8 N14\nM107 N11 N36 N20 N1\nM108 N36 N15 N4 N20\n"""
        """M119 N14 N8 N27 N40\nM120 N40 N27 N4 N15\nFINSF\nSEG2\nS1 N14 N35\nS2 N35 N10"""
        """\nS3 N10 N30\nS4 N30 N11\nS5 N11 N36\nS6 N36 N15\nS7 N15 N40\nS8 N40 N14\nFINSF\n"""
        """GROUP_MA\nCIRCLE\nS1 S2 S3 S4 S5 S6 S7 S8\nFINSF\nGROUP_MA\nSURFACE\nM73  M74  """
        """M75  M76  M87  M88  M101 M102 M107 M108 M119 M120 \nFINSF\nFIN\n""")
        f.file.write(txtMesh.format(radius,radius*sqrt(2)/2.5, radius/sqrt(2)/2.5, radius/sqrt(2)))

        f.file.seek(0)  # file ready to be read
        mesh.readAsterFile(f.name)

    # Definition of the boundary mesh
    nSeg2 = 400
    sNodes = u"TITRE\nBORD\nFINSF\nCOOR_3D\n"
    sNodes += "".join(["N{} {} {} 0\n".format(i, radius*cos(2*pi/(nSeg2+1)*i),
                                 radius*sin(2*pi/(nSeg2+1)*i)) for i in range(nSeg2+1)])
    sNodes += "FINSF\nSEG2\n"
    sNodes += "".join(["M{} N{} N{}\n".format(i, i, i + 1) for i in range(nSeg2)])
    sNodes += "M{} N{} N{}\n".format(nSeg2, nSeg2, 0)
    sNodes += "FINSF\nGROUP_MA\nCIRCLE\n"
    sNodes += "".join(["M{} \n".format(i) for i in range(nSeg2+1)])
    sNodes += "\nFINSF\nFIN\n"

    border = Mesh()
    with tempfile.NamedTemporaryFile(mode="w+") as f:
        f.file.write(sNodes)
        f.file.seek(0)
        border.readAsterFile(f.name)

    # Mesh refinement
    newMesh = mesh
    for _ in range(nrefine):
        HOMARD_INFOS.new()
        resu = MACR_ADAP_MAIL(__use_namedtuple__=True,
                                ADAPTATION='RAFFINEMENT_UNIFORME',
                                MAILLAGE_N=newMesh,
                                MAILLAGE_NP1=CO('newMesh'),
                                MAILLAGE_FRONTIERE=border,
                                GROUP_MA_FRONT='CIRCLE',)
        HOMARD_INFOS.pop()
        newMesh = resu.newMesh
    return newMesh


def buildCube(cls, lx=1, ly=1, lz=1, nrefine=3):
    """Build the hexaedral mesh of a cube.
    Arguments:
        lx [float] : length of the cube along the x axis (default 1.).
        ly [float] : length of the cube along the y axis (default 1.).
        lz [float] : length of the cube along the z axis (default 1.).
        nrefine [int] : number of mesh refinement iterations (default 3).
    """
    # Mesh creation
    mesh = cls()
    with tempfile.NamedTemporaryFile(mode="w+") as f:
        txtMesh = ("""TITRE\n CUBE\nFINSF\nCOOR_3D\nN1 0. {1:f} 0.\nN2 {0:f} {1:f} 0.\n"""
        """N3 0. 0. 0.\nN4 {0:f} 0. 0.\nN5 0. {1:f} {2:f}\nN6 {0:f} {1:f} {2:f}\nN7 0. 0. {2:f}"""
        """\nN8 {0:f} 0. {2:f}\nFINSF\nHEXA8\nCUBE N2 N4 N8 N6 N1 N3 N7\n N5\nFINSF\nQUAD4\n"""
        """BACK N1 N3 N7 N5\nFRONT N2 N6 N8 N4\nRIGHT N5 N6 N2 N1\nLEFT N3 N4 N8 N7\nTOP N7 """
        """N8 N6 N5\nBOTTOM N2 N4 N3 N1\nFINSF\nSEG2\nS13 N1 N3\nS37 N3 N7\nS75 N7 N5\nS26 """
        """N2 N6\nS68 N6 N8\nS84 N8 N4\nS56 N5 N6\nS21 N2 N1\nS34 N3 N4\nS24 N2 N4\nS78 N7 """
        """N8\nS51 N5 N1\nFINSF\nGROUP_MA\nVOLUME CUBE\nFINSF\nGROUP_MA\nBACK BACK\nFINSF\n"""
        """GROUP_MA\nFRONT FRONT\nFINSF\nGROUP_MA\nRIGHT RIGHT\nFINSF\nGROUP_MA\nLEFT LEFT\n"""
        """FINSF\nGROUP_MA\nTOP TOP\nFINSF\nGROUP_MA\nBOTTOM BOTTOM\nFINSF\nGROUP_MA\nS13 S13"""
        """\nFINSF\nGROUP_MA\nS37 S37\nFINSF\nGROUP_MA\nS75 S75\nFINSF\nGROUP_MA\nS26 S26\n"""
        """FINSF\nGROUP_MA\nS68 S68\nFINSF\nGROUP_MA\nS84 S84\nFINSF\nGROUP_MA\nS56 S56\nFINSF"""
        """\nGROUP_MA\nS21 S21\nFINSF\nGROUP_MA\nS34 S34\nFINSF\nGROUP_MA\nS24 S24\nFINSF\n"""
        """GROUP_MA\nS78 S78\nFINSF\nGROUP_MA\nS51 S51\nFINSF\nGROUP_NO\nN1 N1\nFINSF\nGROUP_NO"""
        """\nN2 N2\nFINSF\nGROUP_NO\nN3 N3\nFINSF\nGROUP_NO\nN4 N4\nFINSF\nGROUP_NO\nN5 N5\n"""
        """FINSF\nGROUP_NO\nN6 N6\nFINSF\nGROUP_NO\nN7 N7\nFINSF\nGROUP_NO\nN8 N8"""
        """\nFINSF\nFIN\n""")
        f.file.write(txtMesh.format(lx, ly, lz))
        f.file.seek(0)  # file ready to be read
        mesh.readAsterFile(f.name)

    # Mesh refinement
    newMesh = mesh
    for _ in range(nrefine):
        HOMARD_INFOS.new()
        resu = MACR_ADAP_MAIL(__use_namedtuple__=True,
                                ADAPTATION='RAFFINEMENT_UNIFORME',
                                MAILLAGE_N=newMesh,
                                MAILLAGE_NP1=CO('newMesh'))
        HOMARD_INFOS.pop()
        newMesh = resu.newMesh
    return newMesh


def buildCylinder(cls, height=3, radius=1, nrefine=3):
    """Build the hexaedral mesh of a cylinder.
    Arguments:
        height [float] : height of the cylinder along the z axis (default 3).
        radius [float] : radius of the cylinder (default 1).
        nrefine [int] : number of mesh refinement iterations (default 3).
    """
    # Mesh creation
    mesh = cls()
    with tempfile.NamedTemporaryFile(mode="w+") as f:
        txtMesh = ("""TITRE\nCYLINDRE\nFINSF\nCOOR_3D\nN1 0 -{3:f} 0\nN2 0 -{3:f} {1:f}\n"""
        """N3 -{3:f} 0 {1:f}\nN4 -{3:f} 0 0\nN5 {3:f} 0 0\nN6 {3:f} 0 {1:f}\nN7 0 {3:f} {1:f}"""
        """\nN8 0 {3:f} 0\nN9 {0:f} 0 {1:f}\nN10 {0:f} 0 0\nN11 0 -{0:f} 0\nN12 0 -{0:f} {1:f}"""
        """\nN13 0 {0:f} {1:f}\nN14 0 {0:f} 0\nN15 -{0:f} 0 0\nN16 -{0:f} 0 {1:f}\nN17 0 -{3:f}"""
        """ {2:f}\nN18 -{4:f} -{4:f} {1:f}\nN19 -{3:f} 0 {2:f}\nN20 -{4:f} -{4:f} 0\nN21 {3:f}"""
        """ 0 {2:f}\nN22 {4:f} {4:f} {1:f}\nN23 0 {3:f} {2:f}\nN24 {4:f} {4:f} 0\nN25 {4:f} """
        """-{4:f} 0\nN26 {4:f} -{4:f} {1:f}\nN27 -{4:f} {4:f} 0\nN28 -{4:f} {4:f} {1:f}\nN29 """
        """{0:f} 0 {2:f}\nN30 {5:f} -{5:f} 0\nN31 0 -{0:f} {2:f}\nN32 {5:f} -{5:f} {1:f}\nN33 """
        """{5:f} {5:f} {1:f}\nN34 0 {0:f} {2:f}\nN35 {5:f} {5:f} 0\nN36 -{5:f} -{5:f} 0\nN37 """
        """-{0:f} 0 {2:f}\nN38 -{5:f} -{5:f} {1:f}\nN39 -{5:f} {5:f} {1:f}\nN40 -{5:f} {5:f} 0"""
        """\nN41 -{4:f} -{4:f} {2:f}\nN42 {4:f} {4:f} {2:f}\nN43 {4:f} -{4:f} {2:f}\nN44 -{4:f}"""
        """ {4:f} {2:f}\nN45 0 0 0\nN46 0 0 {1:f}\nN47 {5:f} -{5:f} {2:f}\nN48 {5:f} {5:f} """
        """{2:f}\nN49 -{5:f} -{5:f} {2:f}\nN50 -{5:f} {5:f} {2:f}\nN51 0 0 {2:f}\nFINSF\nQUAD4"""
        """\nM57 N1 N17 N41 N20\nM58 N20 N41 N19 N4\nM59 N17 N2 N18 N41\nM60 N41 N18 N3 N19\n"""
        """M61 N6 N21 N42 N22\nM62 N22 N42 N23 N7\nM63 N21 N5 N24 N42\nM64 N42 N24 N8 N23\n"""
        """M65 N1 N25 N43 N17\nM66 N17 N43 N26 N2\nM67 N25 N5 N21 N43\nM68 N43 N21 N6 N26\nM69"""
        """ N8 N27 N44 N23\nM70 N23 N44 N28 N7\nM71 N27 N4 N19 N44\nM72 N44 N19 N3 N28\nM73 N1 """
        """N20 N45 N25\nM74 N25 N45 N24 N5\nM75 N20 N4 N27 N45\nM76 N45 N27 N8 N24\nM77 N3 N18 """
        """N46 N28\nM78 N28 N46 N22 N7\nM79 N18 N2 N26 N46\nM80 N46 N26 N6 N22\nM81 N10 N29 N47"""
        """ N30\nM82 N30 N47 N31 N11\nM83 N29 N9 N32 N47\nM84 N47 N32 N12 N31\nM85 N6 N9 N29 """
        """N21\nM86 N21 N29 N10 N5\nM87 N10 N30 N25 N5\nM88 N30 N11 N1 N25\nM89 N2 N17 N31 N12"""
        """\nM90 N17 N1 N11 N31\nM91 N12 N32 N26 N2\nM92 N32 N9 N6 N26\nM93 N13 N33 N48 N34\nM94"""
        """ N34 N48 N35 N14\nM95 N33 N9 N29 N48\nM96 N48 N29 N10 N35\nM97 N6 N9 N33 N22\nM98 N22"""
        """ N33 N13 N7\nM99 N7 N13 N34 N23\nM100 N23 N34 N14 N8\nM101 N10 N5 N24 N35\nM102 N35"""
        """ N24 N8 N14\nM103 N11 N31 N49 N36\nM104 N36 N49 N37 N15\nM105 N31 N12 N38 N49\nM106"""
        """ N49 N38 N16 N37\nM107 N11 N36 N20 N1\nM108 N36 N15 N4 N20\nM109 N16 N3 N19 N37\n"""
        """M110 N37 N19 N4 N15\nM111 N16 N38 N18 N3\nM112 N38 N12 N2 N18\nM113 N16 N39 N50 N37"""
        """\nM114 N37 N50 N40 N15\nM115 N39 N13 N34 N50\nM116 N50 N34 N14 N40\nM117 N7 N13 N39 """
        """N28\nM118 N28 N39 N16 N3\nM119 N14 N8 N27 N40\nM120 N40 N27 N4 N15\nFINSF\nHEXA8\n"""
        """M121 N2 N17 N41 N18 N26 N43 N51\n N46\nM122 N26 N43 N51 N46 N6 N21 N42\n N22\nM123 """
        """N18 N41 N19 N3 N46 N51 N44\n N28\nM124 N46 N51 N44 N28 N22 N42 N23\n N7\nM125 N17 N1"""
        """ N20 N41 N43 N25 N45\n N51\nM126 N43 N25 N45 N51 N21 N5 N24\n N42\nM127 N41 N20 N4 """
        """N19 N51 N45 N27\n N44\nM128 N51 N45 N27 N44 N42 N24 N8\n N23\nM129 N9 N29 N47 N32 """
        """N6 N21 N43\n N26\nM130 N32 N47 N31 N12 N26 N43 N17\n N2\nM131 N29 N10 N30 N47 N21 """
        """N5 N25\n N43\nM132 N47 N30 N11 N31 N43 N25 N1\n N17\nM133 N9 N33 N48 N29 N6 N22 N42"""
        """\n N21\nM134 N29 N48 N35 N10 N21 N42 N24\n N5\nM135 N33 N13 N34 N48 N22 N7 N23\n N42"""
        """\nM136 N48 N34 N14 N35 N42 N23 N8\n N24\nM137 N12 N31 N49 N38 N2 N17 N41\n N18\nM138"""
        """ N38 N49 N37 N16 N18 N41 N19\n N3\nM139 N31 N11 N36 N49 N17 N1 N20\n N41\nM140 N49 """
        """N36 N15 N37 N41 N20 N4\n N19\nM141 N13 N39 N50 N34 N7 N28 N44\n N23\nM142 N34 N50 N40"""
        """ N14 N23 N44 N27\n N8\nM143 N39 N16 N37 N50 N28 N3 N19\n N44\nM144 N50 N37 N15 N40 """
        """N44 N19 N4\n N27\nFINSF\nGROUP_MA\nBOTTOM\nM73 M74 M75 M76 M87 M88 M101\nM102 M107 """
        """M108 M119 M120\nFINSF\nGROUP_MA\nTOP\nM77 M78 M79 M80 M91 M92 M97\nM98 M111 M112 M117"""
        """ M118\nFINSF\nGROUP_MA\nVOLUME\nM121 M122 M123 M124 M125 M126 M127 M128\nM129 M130 """
        """M131 M132 M133 M134 M135\nM136 M137 M138 M139 M140 M141 M142\nM143 M144\nFINSF\n"""
        """GROUP_MA\nSURFEXT\nM81 M82 M83 M84 M93 M94 M95 M96 M103 M104 M105 M106 M113 M114  """
        """M115 M116\nFINSF\nFIN\n""")
        f.file.write(txtMesh.format(radius, height, height/2., radius*sqrt(2)/2.5,
        radius/sqrt(2)/2.5, radius/sqrt(2)))
        f.file.seek(0)  # file ready to be read
        mesh.readAsterFile(f.name)

    # Mesh refinement
    newMesh = mesh
    for _ in range(nrefine):
        HOMARD_INFOS.new()
        resu = MACR_ADAP_MAIL(__use_namedtuple__=True,
                              ADAPTATION='RAFFINEMENT_UNIFORME',
                              MAILLAGE_N=newMesh,
                              MAILLAGE_NP1=CO('newMesh'),
                              FRONTIERE_ANALYTIQUE=_F(GROUP_MA='SURFEXT', NOM='SURFEXT',
                                                      TYPE='CYLINDRE', X_CENTRE=0.0, Y_CENTRE=0.0,
                                                      Z_CENTRE=0.0, X_AXE=0., Y_AXE=0., Z_AXE=1.,
                                                      RAYON=radius),
                              )
        HOMARD_INFOS.pop()
        newMesh = resu.newMesh
    return newMesh
