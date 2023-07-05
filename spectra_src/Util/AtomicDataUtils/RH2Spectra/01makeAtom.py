# usage : 
# $ conda activate <spectra-env>
# $ cd /path/to/spectra/data/atom/Ca_I-II-III
# $ python ../../../spectra_src/Util/AtomicDataUtils/RH2Spectra/01makeAtom.py CaI+II_45.atom.RH.20220802.txt CaI+II_45.RH.configuration-table.txt ./ -file-prefix rh.


#from glob import glob
import argparse
from pprint import pprint
from dataclasses import dataclass
import os

from spectra_src import Elements
from spectra_src import Constants as Cst

#-------------------------------------------------------------------
# global constants
#-------------------------------------------------------------------
#c_ = Cst.c_#2.9979246  * 1.E+10         #: speed of light, [:math:`cm \cdot s^{-1}`]
#h_ = Cst.h_#6.62606885 * 1.E-27         #: Planck constant, [:math:`erg \cdot s`]
#eV2erg_ = Cst.eV2erg_#1.60217662 * 1.E-12         #: unit conversion from eV to erg, [:math:`erg/eV`]
CM_INVERSE_TO_EV_ = Cst.c_ * Cst.h_ / Cst.eV2erg_


#-------------------------------------------------------------------
# struct
#-------------------------------------------------------------------

@dataclass
class LevelRecord:
    conf : str
    term : str
    j : str
    g : int
    stage : int
    ev : float

@dataclass
class LevelRH:
    ev : float

@dataclass
class LineRH:
    j: int
    i: int
    f: float
    type: str
    nlambda: int
    symmetr: str
    qcore: float
    qwing: float
    vdWapprx: str
    vdWaals_H_1: float 
    vdWaals_H_2: float 
    vdWaals_He_1: float 
    vdWaals_He_2: float 
    radiative: float
    stark: float

@dataclass
class AlphaRH:
    j: int
    i: int
    alpha0 : float
    nlambda : int
    wavedep : str
    lambmin : float
    waves  : list[float]
    alphas : list[float]

@dataclass 
class ColTempRH:
    ntemp : int
    temps : list[float]

@dataclass
class CollisionRH:
    j : int
    i : int
    kind : str
    coes : list[float]

@dataclass
class AtomFiles:
    level : str
    Aji   : str
    RadiativeLine : str
    Alpha : str
    CEe   : str
    CIe   : str
    Gro   : str 
    Conf  : str


#-------------------------------------------------------------------
# real rh level index --> configuration,term,j mapping table
#-------------------------------------------------------------------


def is_skip_(text: str):
    if text.startswith('#'): return True
    if text == '': return True
    return False

def read_table_(file: str):
    datas : dict[int,LevelRecord] = {}
    with open(file, 'r') as f:
        #count = 0
        prefix = None
        for line in f:
            line = line.strip()
            if is_skip_(line): continue

            words = [w.strip() for w in line.split()]
            if prefix is None and words[0].lower() == 'prefix': 
                prefix = words[1]
                continue
            if prefix is None: continue
            assert len(words)==7, "table must have 7 columns"
            
            idx = int(words[0])
            conf = words[1]
            term = words[2]
            J    = words[3]
            g    = int(words[4])
            stage= int(words[5])
            if '+' in words[6]:
                ev = sum([float(v) for v in words[6].split('+')])
            else:
                ev   = float( words[6] )
            # datas[idx] = {
            #     'conf': conf,
            #     'term': conf,
            #     'J'   : J,
            #     'g'   : g,
            #     'stage': stage,
            #     'ev' : ev,
            # }
            datas[idx] = LevelRecord(
                conf = conf, term=term, j=J, g=g, stage=stage, ev=ev
            )
            #count += 1
    return datas, prefix

#-------------------------------------------------------------------
# read RH level table and make .Level
#-------------------------------------------------------------------
def read_rh_level_(file: str):
    
    def is_start_(text: str):
        if not text.startswith('#'): return False
        if 'Energy levels' in text : return True
        return False

    def is_end_(text: str):
        if text.startswith('#') and \
            'qcore' in text: return True
        return False

    datas: dict[int,LevelRH] = {}
    with open(file, 'r') as f:
        isstart = False
        for line in f:
            line = line.strip()
            ##: check whether to start reading
            if not isstart:
                if is_start_(line):
                    isstart = True
                continue
            ##: check whether to end reading
            if is_end_(line): break
            ##: check whether to skip line
            if is_skip_(line): continue
            
            ##: read
            words = line.split("'")
            assert len(words)==3, "text splited by ' must have length==3"
            label = words[1]
            line = words[0] + words[2]
            words = [w.strip() for w in line.split()]

            cm_inv = float(words[0])
            #g      = int(float(words[1]))
            idx    = int(words[3])
            ev = 0.0 if cm_inv==0.0 else cm_inv * CM_INVERSE_TO_EV_   
            level = LevelRH(ev=ev)
            datas[idx] = level          

    return datas

COMMENT_ = '#'
END_ = 'END'
def text_bar_(n:int=98):
    return COMMENT_+'-'*n

def get_snum_(text: str):
    s = ''
    for c in text:
        if not c.isdigit(): break
        s += c
    return s

def conf2n_(conf: str):
    if conf == '-': return '-'
    s = conf.split('.')[-1]
    return get_snum_(s)

def term2s2p1_(term: str):
    if term == '-': return '-'
    return get_snum_(term)

def term2L_(term: str):
    if term == '-': return '-'
    sym = term[-1]
    return Cst.L_s2i_[sym]


def make_Level_file_( levels: dict[int, LevelRH], tables: dict[int,LevelRecord], 
                    outfile: str, title: str, sym: str, prefix: str ):

    with open(outfile, 'w') as f:
        f.write(text_bar_() + '\n')
        f.write(f"Title: {title}" + '\n')
        f.write(text_bar_() + '\n')
        
        for key, val in zip(('Z','Element','nLevel'),(Elements.ELEMENT_DICT[sym]['Z'],sym,len(tables))):
            f.write(f"{key:<20s}{val}\n")
        f.write(COMMENT_+'\n')
        f.write(END_ + '\n')
        f.write(text_bar_() + '\n')

        ## level table
        f.write(f"{'prefix':<16s}{prefix}\n")
        text = "#"+' '*3
        names = 'conf', 'term', 'J', 'n', 'L', '2S+1', 'g=2J+1', 'stage', 'E[eV]'
        for name in names:
            text += f'{name:<12s}'
        f.write(text + '\n')

        for k, rec in tables.items():
            n = conf2n_(rec.conf)
            L = term2L_(rec.term)
            s2p1 = term2s2p1_(rec.term)
            
            values = (
                rec.conf, rec.term, rec.j, n, L, s2p1, rec.g, rec.stage, rec.ev 
            )
            text = ' '*4
            for v in values:
                if isinstance(v, str): 
                    text += f'{v:<12s}'
                    continue
                if isinstance(v, int): 
                    text += f'{v:<12d}'
                    continue
                if isinstance(v, float): 
                    text += f'{v:+.12E}'
                    continue
                raise TypeError(f"unsupported type(v)={type(v)}")
            f.write(text + '\n')
        
        f.write(END_+'\n')
        f.write(text_bar_()+'\n')

    return 0

#-------------------------------------------------------------------
# read RH line table and make .Aji and .RadiativeLine
#-------------------------------------------------------------------
def read_rh_line_(file: str):
    
    def is_start_(text: str):
        if not text.startswith('#'): return False
        if 'qcore' in text : return True
        return False

    def is_end_(text: str):
        if text.startswith('#') and \
            'alpha [m^-2]' in text: return True
        return False

    nwords = 15
    datas: list[LineRH] = []
    with open(file, 'r') as f:
        isstart = False
        for line in f:
            line = line.strip()
            ##: check whether to start reading
            if not isstart:
                if is_start_(line):
                    isstart = True
                continue
            ##: check whether to end reading
            if is_end_(line): break
            ##: check whether to skip line
            if is_skip_(line): continue
            
            ##: read        
            words = [w.strip() for w in line.split()]
            assert len(words)==nwords, f"text splited by whitespace must have length=={nwords}"
            
            j = int(words[0])
            i = int(words[1])
            ff= float(words[2])
            t = words[3]
            nlambda = int(words[4])
            symmetr = words[5]
            qcore = float(words[6])
            qwing = float(words[7])
            vdWapprx = words[8]
            vdWaals_H_1 = float(words[9])
            vdWaals_H_2 = float(words[10])
            vdWaals_He_1= float(words[11])
            vdWaals_He_2= float(words[12])
            radiative   = float(words[13])
            stark       = float(words[14])

            linerh = LineRH(
                j=j,i=i,f=ff,type=t,nlambda=nlambda,symmetr=symmetr,
                qcore=qcore,qwing=qwing,vdWapprx=vdWapprx,
                vdWaals_H_1=vdWaals_H_1,vdWaals_H_2=vdWaals_H_2,
                vdWaals_He_1=vdWaals_He_1,vdWaals_He_2=vdWaals_He_2,
                radiative=radiative,stark=stark
            )
            datas.append( linerh )   

    return datas

def make_Aji_file_( lines: list[LineRH], 
                    tables: dict[int,LevelRecord], 
                    outfile: str, prefix: str):
    with open(outfile, 'w') as f:
        f.write(f"{'prefix':<16s}{prefix}\n")
        
        text = "#"+' '*3
        names = (
            'conf_i', 'term_i', 'J_i',
            'conf_j', 'term_j', 'J_j',
            'Aji[s^-1]', 'Wave[AA]'
        )
        for name in names:
            text += f'{name:<12s}'
        f.write(text + '\n')

        for rhline in lines:
            j = rhline.j
            i = rhline.i
            Aji = rhline.radiative

            confj = tables[j].conf
            termj = tables[j].term
            Jj    = tables[j].j
            confi = tables[i].conf
            termi = tables[i].term
            Ji    = tables[i].j
            
            evj   = tables[j].ev 
            evi   = tables[i].ev
            wave_AA = Cst.eV2AA_div_ / abs(evj-evi)
            text = ' '*4
            text += f"{confi:<12s}"
            text += f"{termi:<12s}"
            text += f"{Ji:<12s}"
            text += f"{confj:<12s}"
            text += f"{termj:<12s}"
            text += f"{Jj:<12s}"
            text += f"{Aji:.4E}  "
            text += f"{wave_AA:.5E} "
            f.write(text + '\n')

def make_RadiativeLine_file_( lines: list[LineRH], 
                    tables: dict[int,LevelRecord], 
                    outfile: str, prefix: str):
    with open(outfile, 'w') as f:
        f.write(text_bar_() + '\n')
        nRadiative = len( lines )
        f.write(f"nRadiative    {nRadiative}" + '\n')
        f.write(END_ + '\n')
        f.write(text_bar_() + '\n')
        f.write(f"{'prefix':<16s}{prefix}\n")

        text = "#"+' '*3
        names = (
            'conf_i','term_i','J_i',
            'conf_j','term_j','J_j',
            'Emission','Nlambda','BackWidth[cm/s]',
            'Fixed','qcore','qwing','filename'
        )
        for name in names:
            if name == 'BackWidth[cm/s]': text += f'{name:<18s}' 
            else : text += f'{name:<12s}'
        f.write(text + '\n')

        for rhline in lines:
            j = rhline.j
            i = rhline.i

            confj = tables[j].conf
            termj = tables[j].term
            Jj    = tables[j].j
            confi = tables[i].conf
            termi = tables[i].term
            Ji    = tables[i].j

            emission = rhline.type.capitalize()
            nlambda  = rhline.nlambda
            backwidth= '5.E+07'
            fixed    = '0'
            qcore    = rhline.qcore
            qwing    = rhline.qwing
            filename = 'tmp.dat'

            if nlambda%2 == 0: nlambda += 1

            text = ' '*4
            text += f"{confi:<12s}"
            text += f"{termi:<12s}"
            text += f"{Ji:<12s}"
            text += f"{confj:<12s}"
            text += f"{termj:<12s}"
            text += f"{Jj:<12s}"
            text += f"{emission:<12s}"
            text += f"{nlambda:<12d}"
            text += f"{backwidth:<18s}"
            text += f"{fixed:<12s}"
            s = f"{qcore:.2f}"
            text += f"{s:<12s}"
            s = f"{qwing:.2f}"
            text += f"{s:<12s}"
            text += f"{filename:<12s}"
            f.write(text + '\n')
        f.write(END_ + '\n')
    return 0
#-------------------------------------------------------------------
# read RH photoionization table and make .Alpha
#-------------------------------------------------------------------
def read_rh_alpha_(file: str):
    
    def is_start_(text: str):
        if not text.startswith('#'): return False
        if 'alpha [m^-2]' in text : return True
        return False

    def is_end_(text: str):
        if not text.startswith('#') and \
            'TEMP' in text: return True
        return False

    nwords = 6
    datas: list[AlphaRH] = []
    with open(file, 'r') as f:
        isstart = False
        isblock = False
        for line in f:
            line = line.strip()
            ##: check whether to start reading
            if not isstart:
                if is_start_(line):
                    isstart = True
                continue
            ##: check whether to end reading
            if is_end_(line): break
            ##: check whether to skip line
            if is_skip_(line): continue
            
            ##: read        
            words = [w.strip() for w in line.split()]
            if not isblock and len(words)==nwords:
                j = int(words[0])
                i = int(words[1])
                alpha0 = float(words[2]) * 1.E4  # [m^2] --> [cm^2]
                nlambda = int(words[3])
                wavedep = words[4]
                lambmin = float(words[5])# [nm]
                rhalpha = AlphaRH(
                    j=j,i=i,alpha0=alpha0,nlambda=nlambda,
                    wavedep=wavedep,lambmin=lambmin,
                    waves=[],alphas=[],
                )
                if wavedep != 'EXPLICIT': 
                    datas.append( rhalpha ) 
                    continue

                isblock = True
                count = 0
                continue
            if not isblock: continue
            assert len(words)==2, f"text splited by whitespace must have length==2"
            wave  = float(words[0])   # [nm]
            alpha = float(words[1]) * 1.E4  # [m^2] --> [cm^2]
            rhalpha.waves.append( wave )
            rhalpha.alphas.append( alpha )

            count += 1
            if count == rhalpha.nlambda:
                assert len(rhalpha.alphas)==rhalpha.nlambda
                isblock = False
                datas.append( rhalpha )   

    return datas

def make_Alpha_file_( alphas: list[AlphaRH], 
                    tables: dict[int,LevelRecord], 
                    outfile: str, prefix: str):
    header = """
#-----------------------------------------------------------------------------
# Note :
#
#       1. if Hydrogenic, the photoinization cross section at other frequency:
#
#               alpha(nu) = alpha0 * (nu_0 / nu)^3
#
#           so the bound-free Gount factor is constant for each line.
#           and the J effect on bound-free3p6.3d --> 3p6 Gaunt factor is ignored.
#
#       2. if data, the photoionization cross section is calculated by interpolation
#-----------------------------------------------------------------------------
#   i : lower level
#   j : upper level (continuum)
#   nLambda : # frequency mesh
#   alpha0[cm^2] : Photoionization cross section at frequency edge
#   Wavelength Dependence :
#       - Hydrogenic    : 0     -->     nu^{-3} dependency
#       - data          : 1     -->     reads data from experimental result
#-----------------------------------------------------------------------------
    """
    depenmap = {
        "EXPLICIT" : "data",
        "HYDROGENIC" : "Hydrogenic"
    }
    with open(outfile, 'w') as f:
        f.write(header + '\n')
        nCont = len( alphas )
        f.write(f"nCont    {nCont}"+'\n')
        f.write(END_ + '\n')
        f.write(text_bar_() + '\n')
        f.write(f"{'prefix':<16s}{prefix}\n")
        


        text = "#"+' '*3
        names = (
            'conf_i', 'term_i', 'J_i',
            'conf_j', 'term_j', 'J_j',
            'nLambda', 'Wave Depen', 'alpha0[cm^2]'
        )
        for name in names:
            text += f'{name:<12s}'
        f.write(text + '\n')
        f.write(COMMENT_ + '\n')
        text = "#"+' '*3
        text += f"wavelength[nm]    cross section [cm^2]"
        f.write(text+'\n')
        f.write(COMMENT_ + '\n')

        for rhalpha in alphas:
            j = rhalpha.j
            i = rhalpha.i
            confj = tables[j].conf
            termj = tables[j].term
            Jj    = tables[j].j
            confi = tables[i].conf
            termi = tables[i].term
            Ji    = tables[i].j
            nlambda = rhalpha.nlambda
            wavedep = depenmap[rhalpha.wavedep]
            alpha0  = rhalpha.alpha0

            text = ' '*4
            text += f"{confi:<12s}"
            text += f"{termi:<12s}"
            text += f"{Ji:<12s}"
            text += f"{confj:<12s}"
            text += f"{termj:<12s}"
            text += f"{Jj:<12s}"
            text += f"{nlambda:<12d}"
            text += f"{wavedep:<12s}"
            text += f"{alpha0:.4E}"
            f.write(text + '\n')
            f.write(COMMENT_ + '\n')
            for w, a in zip( rhalpha.waves, rhalpha.alphas ):
                s1 = f"{w:.1f}"
                s2 = f"{a:.4E}"
                f.write(f"{s1:>18s}    {s2:>12s}" + '\n')
            f.write(COMMENT_ + '\n')
        f.write(END_ + '\n')

    return 0

#-------------------------------------------------------------------
# read RH collisional excitation/ionization coefficient table 
# and make .CE.Electron and .CI.Electron 
#-------------------------------------------------------------------
def read_rh_collision_(file: str):
    
    def is_start_(text: str):
        if text.startswith('#'): return False
        if 'TEMP' in text : return True
        return False

    def is_end_(text: str):
        if not text.startswith('#') and \
            'END' in text: return True
        return False

    nwords = 6
    datas: list[AlphaRH] = []
    with open(file, 'r') as f:
        isstart = False

        for line in f:
            line = line.strip()
            ##: check whether to start reading
            if not isstart:
                if is_start_(line):
                    isstart = True
                    words = [w.strip() for w in line.split()]
                    ntemp = int(words[1])
                    rhtemp = ColTempRH(
                        ntemp=ntemp,
                        temps=[float(w) for w in words[2:2+ntemp]])
                continue
            ##: check whether to end reading
            if is_end_(line): break
            ##: check whether to skip line
            if is_skip_(line): continue
            


            ##: read        
            words = [w.strip() for w in line.split()]
            
            kind = words[0]
            j = int(words[1])
            i = int(words[2])
            fac = 1.E+6 if kind in ('CE','CI') else 1.0
            rhcol = CollisionRH(
                j=j,i=i,kind=kind,
                coes=[float(w)*fac for w in words[3:3+ntemp]] 
            )
            ##: CE & CI [m^3 K^-1/2] --> *1.E-6 [cm^3 K^-1/2]
            datas.append( rhcol )   

    return datas, rhtemp

def make_CE_file_( temp: ColTempRH, cols: list[CollisionRH], 
                    tables: dict[int,LevelRecord], 
                    outfile: str, prefix: str):
    header = """
#--------------------------------------------------------------------------------------------------
# Collisional Excitation Rate
# type :
#       ECS -> Effective Collision Strength
#       CRC -> Collision Rate Coefficient
# Effective Collision Strengths (Omega_ij[-]) : RH/Atoms_example/CaII.atom
# Collision Rate : n_e * C_ij = n_e * ( 8.63e-6 * (Omega_ij*f1/f2) ) / (g_i * T_e^0.5) * e^(-dE_ji / (k*T_e))
#                             = n_e *        CE * T_e^0.5 * f1/f2                      * e^(-dE_ji / (k*T_e)) 
#  --> Omega_ij = CE * T_e * g_i / 8.63E-6
"""
    with open(outfile, 'w') as f:
        f.write(header + '\n')
        f.write(f"Type   ECS"+'\n')

        text = f"Temperature    "
        for t in temp.temps:
            s = f'{t:.2E}'
            text += f"{s:<12s}"
        f.write(text + '\n')

        f.write(END_ + '\n')
        f.write(text_bar_() + '\n')
        f.write(f"{'prefix':<16s}{prefix}\n")
        


        text = "#"+' '*3
        names = (
            'conf_i', 'term_i', 'J_i',
            'conf_j', 'term_j', 'J_j',
        )
        for name in names:
            text += f'{name:<12s}'
        for t in temp.temps:
            s = f'{t:.2E}'
            text += f"{s:<12s}"
        for name in ('f1','f2'):
            text += f"{name:<8s}"
        f.write(text + '\n')

        for rhcol in cols:
            j = rhcol.j
            i = rhcol.i
            confj = tables[j].conf
            termj = tables[j].term
            Jj    = tables[j].j
            confi = tables[i].conf
            termi = tables[i].term
            Ji    = tables[i].j
            kind  = rhcol.kind
            coes0  = rhcol.coes
            
            gi    = tables[i].g
            coes = []
            if not kind in ('CE','OMEGA'): continue
            if kind == 'CE':
                for t, c in zip( temp.temps, coes0 ):
                    coes.append( c * gi * t / 8.63E-6 )
            else:
                coes = coes0
            
            text = ' '*4
            text += f"{confi:<12s}"
            text += f"{termi:<12s}"
            text += f"{Ji:<12s}"
            text += f"{confj:<12s}"
            text += f"{termj:<12s}"
            text += f"{Jj:<12s}"

            for c in coes:
                s = f"{c:.3E}"
                text += f"{s:<12s}"
            for v in (1,1):
                text += f"{v:<8d}"
            f.write(text + '\n')

    return 0

def make_CI_file_( temp: ColTempRH, cols: list[CollisionRH], 
                    tables: dict[int,LevelRecord], 
                    outfile: str, prefix: str, nCont: int):
    header = """
#--------------------------------------------------------------------------------------------------
# Collisional Ionization Rate
# CI [s^-1 K^-1/2 cm^3] : RH/Atoms_example/CaII.atom
# Collision Rate : n_e * C_ik = n_e * CI * exp(-dE/kT) * sqrt(T) / f2
"""
    with open(outfile, 'w') as f:
        f.write(header + '\n')
        f.write(f"nCont   {nCont}" + '\n')

        text = f"Temperature    "
        for t in temp.temps:
            s = f'{t:.2E}'
            text += f"{s:<12s}"
        f.write(text + '\n')

        f.write(END_ + '\n')
        f.write(text_bar_() + '\n')
        f.write(f"{'prefix':<16s}{prefix}\n")
        


        text = "#"+' '*3
        names = (
            'conf_i', 'term_i', 'J_i',
            'conf_j', 'term_j', 'J_j',
        )
        for name in names:
            text += f'{name:<12s}'
        for t in temp.temps:
            s = f'{t:.2E}'
            text += f"{s:<12s}"
        for name in ('f2',):
            text += f"{name:<8s}"
        f.write(text + '\n')

        for rhcol in cols:
            j = rhcol.j
            i = rhcol.i
            confj = tables[j].conf
            termj = tables[j].term
            Jj    = tables[j].j
            confi = tables[i].conf
            termi = tables[i].term
            Ji    = tables[i].j
            kind  = rhcol.kind
            coes  = rhcol.coes
            
            if kind != 'CI': continue
            
            text = ' '*4
            text += f"{confi:<12s}"
            text += f"{termi:<12s}"
            text += f"{Ji:<12s}"
            text += f"{confj:<12s}"
            text += f"{termj:<12s}"
            text += f"{Jj:<12s}"

            for c in coes:
                s = f"{c:.3E}"
                text += f"{s:<12s}"
            for v in (1,):
                text += f"{v:<8d}"
            f.write(text + '\n')

    return 0
#-------------------------------------------------------------------
# make .Grotrian file 
#-------------------------------------------------------------------
def make_Gro_file_( outfile: str, prefix: str):
    
    with open(outfile, 'w') as f:
        
        f.write(text_bar_() + '\n')
        f.write(f"{'prefix':<16s}{prefix}\n")
        f.write(text_bar_() + '\n')
        f.write("position")
        f.write(text_bar_() + '\n')
        f.write("lineplot")

        text = "#"+' '*3
        names = (
            'conf_i', 'term_i', 'J_i',
            'conf_j', 'term_j', 'J_j',
            'text',
        )
        for name in names:
            text += f'{name:<12s}'
        for name in ('r1','r2'):
            text += f"{name:<8s}"
        f.write(text + '\n')

        f.write(END_ + '\n')
        f.write(text_bar_() + '\n')
    return 0

#-------------------------------------------------------------------
# make .conf file
#-------------------------------------------------------------------
def make_conf_file_( afiles:AtomFiles, outfile: str, outdir: str):

    names = (
        "folder",
        "Level",
        "Aji",
        "CEe",
        "CIe",
        "PI",
        "RadiativeLine",
        "Grotrian"
    )
    fnames = (
        os.path.abspath(outdir),
        afiles.level,
        afiles.Aji,
        afiles.CEe,
        afiles.CIe,
        afiles.Alpha,
        afiles.RadiativeLine,
        afiles.Gro,
    )
    with open(outfile, 'w') as f:
        count = 0
        for na, fna in zip( names, fnames ):
            if count ==0: fna = fna.replace('\\', '\\\\')
            if count > 0: fna = os.path.basename(fna)
            text = f"{na:<20s}{fna:<50s}"
            f.write(text + '\n')
            count += 1
    return 0
#-------------------------------------------------------------------
# main function
#-------------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(description='convert RH *.atom to *.Level configuration file', add_help=True)
    parser.add_argument('rhpath', metavar='rhpath', type=str,
                    help='path of rh *.atom file')
    parser.add_argument('conftable', metavar='conftable', type=str,
                    help='path of RH level index --> *.configuration-table.txt')
    parser.add_argument('symbol', metavar='symbol', type=str,
                    help='atomic symbol of the element')
    parser.add_argument('--outdir', metavar='DIR', 
                    type=str, default='./',
                    help='directory to output configuration files')
    parser.add_argument('--title', metavar='title',
                    type=str, default='please modify this title line',
                    help='title information')
    parser.add_argument('--use-rh-energy', 
                    action='store_true', default=False,
                    help='whether to use RH level energy value (wavelength in air?), defualt: False')
    parser.add_argument('--file-prefix', metavar='fprefix',
                    type=str, default='tmp.',
                    help='prefix of generated configuration files')
    parser.add_argument('--indir', metavar="indir",
                    type=str, default='',
                    help='if directory is defined, then use directory/{basename of rhpath/conftable}')
    parser.add_argument('--degenerate', metavar='degeneratefile',
                    type=str, default='',
                    help='json file to define what levels should be degenerated')
    
    args = parser.parse_args()


    ##: step 0. check atomic symbol
    sym = args.symbol
    if sym not in Elements.ELEMENT_DICT.keys():
        raise ValueError(f"symbol [{sym}] not found in :\n{Elements.ELEMENT_DICT.keys()}")

    ##: step 1 create output folder if not exist
    outdir = args.outdir
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir): os.makedirs(outdir, exist_ok=True)

    ##: step 2. create struct for output filenames
    
    fn = f'{args.file_prefix}{sym}'
    afile = AtomFiles(
        level  = os.path.join(outdir, f'{fn}.Level'),
        Aji    = os.path.join(outdir, f'{fn}.Aji'),
        RadiativeLine = os.path.join(outdir, f'{fn}.RadiativeLine'),
        Alpha = os.path.join(outdir, f'{fn}.Alpha'),
        CEe   = os.path.join(outdir, f'{fn}.CE.electron'),
        CIe   = os.path.join(outdir, f'{fn}.CI.electron'),
        Gro   = os.path.join(outdir, f'{fn}.Grotrian'),
        Conf  = os.path.join(outdir, f'{fn}.conf'),
    )

    if len(args.indir) > 0:
        args.indir = os.path.abspath(args.indir)
        args.rhpath = os.path.join(args.indir, os.path.basename(args.rhpath))
        args.conftable = os.path.join(args.indir, os.path.basename(args.conftable))

    ##: step 3. read level index --> configuration, term, j table
    tables, prefix = read_table_(args.conftable)
    #pprint(tables)

    ##: step 4. read RH levels
    rhlevels = read_rh_level_(args.rhpath)

    ##: step 5. use rh level energy or NIST level energy
    if args.use_rh_energy:
        for k in tables.keys(): tables[k].ev = rhlevels[k].ev

    ##: step 6. make .Level file
    make_Level_file_(rhlevels, tables, afile.level ,args.title,sym, prefix)
    
    ##: step 7. read RH lines
    rhlines = read_rh_line_(args.rhpath)
    #pprint(rhlines)

    ##: step 8. make .Aji file
    make_Aji_file_(rhlines, tables, afile.Aji, prefix)

    ##: step 9. make .RadiativeLine file
    make_RadiativeLine_file_(rhlines, tables, afile.RadiativeLine, prefix)

    ##: step 10. read RH photoionization data
    rhalphas = read_rh_alpha_(args.rhpath)

    ##: step 11. make .Alpha file
    make_Alpha_file_(rhalphas, tables, afile.Alpha, prefix)
    
    ##: step 12. read RH collisional Excitaion/Ionization data
    rhcols, rhcoltemp = read_rh_collision_(args.rhpath)

    ##: step 13. make .CE.Electron file
    make_CE_file_(rhcoltemp, rhcols, tables, afile.CEe, prefix)

    ##: step 14. make .CI.Electron file
    nCont = len(rhalphas)
    make_CI_file_(rhcoltemp, rhcols, tables, afile.CIe, prefix, nCont)

    ##: step 15. make .Grotrian file
    make_Gro_file_( afile.Gro, prefix)

    ##: step 16. make .conf file
    make_conf_file_( afile, afile.Conf, outdir)

    ##: step 17. check level combine
    if len(args.degenerate) > 0:
        print(f"to degenerate levels with conf file : {args.degenerate}")    
    

    return 0


if __name__ == '__main__':
    main()