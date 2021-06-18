

#-------------------------------------------------------------------------------
# function/class definition of Grotrian Plotting
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0 
#    2021/06/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

import matplotlib.pyplot as plt
from collections import Counter
import math
#from matplotlib.cm import ScalarMappable

from ..Util.AtomUtils import AtomIO, AtomInfo

from . import Plotting


def _prepare_dict_(_atom , _conf_prefix, _scaleFunc , _exclude ):
    r"""
    separate singlet and multiplet

        - singlet   : { (conf_clean,term) : { J : (E[eV], L, conf_origin) } }
        - multiplet : { (conf_clean,term) : { J : (E[eV], L, conf_origin) } }

    Parameters
    ----------
    _atom : AtomCls.Atom
        object of the atomic model

    _conf_prefix : str
        common configuration string of the inner shell

    Returns
    -------

    _singlet : { (conf_clean,term) : { J : (E[eV], L, conf_origin) } }

    _multiplet : { (conf_clean,term) : { J : (E[eV], L, conf_origin) } }

    _Lset : { "singlet" : {int,}, "multiplet" : {int,} }
        set of Quantum number L in integer for singlet and multiplet, respectively.
    """    
    _L_prefix = len(_conf_prefix)

    _Level_info = AtomInfo.Level_ctj_table_to_Level_info_( _atom._ctj_table.Level )
    #-------------------------------------------------------------------------
    # create and count list of (conf, term)
    #-------------------------------------------------------------------------
    _conf_term = list( zip( _Level_info["configuration"], _Level_info["term"] ) )
    _count = Counter( _conf_term )
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    # separate singlet and multiplet
    #-------------------------------------------------------------------------
    _L_prefix = len(_conf_prefix)

    _singlet = {}
    _multiplet = {}
    _Lset = {"singlet" : set(), "multiplet" : set()}

    for k in range(_atom.nLevel):


        #---------------------------------------------------------------------
        # - remove inner shell electron configuration
        # - remove all '.' in conf
        #---------------------------------------------------------------------
        _conf = _Level_info["configuration"][k]
        _conf_clean = _conf[_L_prefix:].replace('.','')
        #---------------------------------------------------------------------
        _term = _Level_info["term"][k]
        _J = _Level_info["J"][k]

        try:
            _L = CST.L_s2i_[ _term[-1] ]
        except:
            _L = 0

        if _count[ (_conf, _term) ] == 1:
            _d = _singlet
            _Lset["singlet"].add(_L)
        else:
            _d = _multiplet
            _Lset["multiplet"].add(_L)

        _key = (_conf_clean ,_term)

        if _scaleFunc is not None:
            _y = _atom.Level["erg"][k] / _atom.Level["erg"][:].max()
            _y = _scaleFunc(_y)
        else:
            _y = _atom.Level["erg"][k] / CST.eV2erg_

        if _key not in _d.keys():
            _d[_key] = {}

        #---------------------------------------------------------------------
        # assign y_ctj to each excluded ctj
        #---------------------------------------------------------------------
        _ctj = ( _conf, _term, _J )
        if _ctj in _exclude.keys():
            _exclude[_ctj].append( _y )
        else:
            ## add to normal plot level only when the ctj is not excluded
            _d[_key][_J] = ( _y, _L, _conf )

    #-------------------------------------------------------------------------
    # exclude[ctj] = [text1, text2, x, y_ctj]
    # -->group_dict[ (_text_list[0], _text_list[1]) ] = { "count":float,"y_group":float,"x":float,ctj_list=[...]}
    #-------------------------------------------------------------------------

    _group_dict = {}
    for _ctj, _text_list in _exclude.items():
        _group = ( _text_list[0], _text_list[1] )
        try:
            _group_dict[ _group ]['count'] += 1
            _group_dict[ _group ]['y_group'] += _text_list[-1]
            _group_dict[ _group ]['ctj_list'].append( _ctj )
        except KeyError:
            _group_dict[ _group ] = {}
            _group_dict[ _group ]['count'] = 1
            _group_dict[ _group ]['y_group'] = _text_list[-1]
            _group_dict[ _group ]['x_group'] = _text_list[-2]
            _group_dict[ _group ]['ctj_list'] = [_ctj,]
    #-------------------------------------------------------------------------

    return _singlet, _multiplet, _Lset, _group_dict

def line_with_text_(_ax, _lp, _rp, _text, _tsize, _r, _tangle=0, _lcolor="black",
                   _lwidth=0.7, _lstyle='-', _tcolor="black"):
    r"""
    Parameters
    ----------

    _ax : matlotlib.pyplot.Axe
        the axe to plot a line

    _lp : tuple of float/int, (x,y)
        left point of the line

    _rp : tuple of float/int, (x,y)
        right point of the line

    _text : str
        string of the text

    _tsize : int
        texture size

    _r : float, in the range of [0 : 1]
        relative position along the line starting from left point

    _tangle : int or float
        angle to rotate the text, default : 0

    _lcolor : str
        color of the line, default : "black"

    _lwidth : float
        line width, default : 0.7

    _lstyle : str
        line style, default : '-'

    _tcolor : str
        color of the text, default : "black"

    Returns
    -------

    _line_obj : matplotlib.lines.Line2D
        object of the line we plot

    _text_obj : matplotlib.text.Text
        object of the text we plot
    """
    # {'color': 'black', 'fontsize': 6, 'ha': 'center', 'va': 'center',
    #    'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)}
    _line_obj, = _ax.plot([_lp[0], _rp[0]], [_lp[1], _rp[1]], linestyle=_lstyle, linewidth=_lwidth, color=_lcolor)

    if _text is not None:
        _tx = _lp[0] + (_rp[0]-_lp[0]) * _r
        _ty = _lp[1] + (_rp[1]-_lp[1]) * _r
        _text_obj = _ax.text( _tx, _ty, _text, {'color': _tcolor, 'fontsize': _tsize,
                            'ha': 'center', 'va': 'center','rotation':_tangle,
                            'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)} )
    else:
        _text_obj = None

    return _line_obj, _text_obj

def arrow_without_text_(_ax, _pi, _pj, _direction, _cmap, _norm, _v, _abserr, _asize=10, _lwidth=2):
    r"""
    Parameters
    ----------

    _ax : matlotlib.pyplot.Axe
        the axe to plot a arrow

    _pi : tuple of float/int, (x,y)
        xy of point of level i

    _pj : tuple of float/int, (x,y)
        xy of point of level j

    _direction :
        "Left->Right", "Right->Left"

    _lwidth : float
        line width, default : 2

    _asize : float/int
        arrow size, default : 10

    _cmap :
        colormap to be used

    _norm :
        normalizaiton of colormap to to used

    _v : float
        value to decide the color

    _abserr : float
        absolute error, if 0 <= _v < _abserr, skip plotting arrow

    Returns
    -------

    _line_obj : matplotlib.text.Annotation
        object of the arrow we plot
    """
    # {'color': 'black', 'fontsize': 6, 'ha': 'center', 'va': 'center',
    #    'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)}

    assert _direction in ("i->j", "j->i")

    if _direction == "i->j":
        _xytext, _xy = _pi, _pj
    else:
        _xytext, _xy = _pj, _pi

    if 0 <= abs(_v) < _abserr:
        return None
    elif _v < 0:
        _xytext, _xy = _xy, _xytext
        _v = -_v

    _annotation_obj = _ax.annotate('', xy=_xy, xycoords='data',
                                   xytext=_xytext,textcoords='data',
                                   arrowprops=dict(color=_cmap(_norm(_v)),width=_lwidth, headwidth=_asize))

    return _annotation_obj

def read_Grotrian_(_lns):
    r"""
    read default line connection setup for Grotrian diagram
    """
    _line_plot  = {}
    _position   = {}
    _exclude    = {}
    #_prefix = ''
    _isPosition, _isLineplot = False, False
    _isGroup = False
    for _i, _ln in enumerate(_lns[:]):

        if AtomIO.skip_line_(_ln):
            continue
        elif AtomIO.check_end_(_ln):
            break

        _words = _ln.split()
        _words = [_v.strip() for _v in _words]

        if _words[0] == "prefix":
            _prefix = _words[1] if _words[1] != '-' else ''
            continue
        elif _words[0] == 'position':
            #_position = {}
            _isPosition, _isLineplot = True, False
            _isGroup = False
            continue
        elif _words[0] == 'lineplot':
            _isPosition, _isLineplot = False, True
            _isGroup = False
            continue
        elif _words[0] == 'group':
            _isPosition, _isLineplot = False, False
            _isGroup = True
            _group_str = _ln
            #_group_text = [ _words[1], _words[2], float(_words[3]) ]
            continue

        if _isPosition or _isGroup:

            if _words[0] == '-':
                if _prefix == '':
                    _ctj = (_words[0],_words[1],_words[2])
                elif _prefix[-1] == '.':
                    _ctj = (_prefix[:-1],_words[1],_words[2])
                else:
                    assert False
            else:
                _ctj = (_prefix+_words[0],_words[1],_words[2])

            if _isPosition:
                _position[_ctj] = float(_words[3])
            elif _isGroup:
                _exclude[_ctj] = [_v.strip() for _v in _group_str.split()][1:]


        if _isLineplot:
            _params = []
            # get ctj pair
            if _words[3] == '-':
                if _prefix == '':
                    _ctj_ij = ( (_prefix+_words[0],_words[1],_words[2]), (_words[3],_words[4],_words[5]) )
                elif _prefix[-1] == '.':
                    _ctj_ij = ( (_prefix+_words[0],_words[1],_words[2]), (_prefix[:-1],_words[4],_words[5]) )
                else:
                    assert False
            else:
                _ctj_ij = ( (_prefix+_words[0],_words[1],_words[2]), (_prefix+_words[3],_words[4],_words[5]) )

            _params.append( _ctj_ij[0] )
            _params.append( _ctj_ij[1] )
            _params.append( _words[6] )
            _params.append( float(_words[7]) )
            _params.append( float(_words[8]) )

            #_line_plot.append(_params)
            _line_plot[_ctj_ij] = _params[2:]

    return _line_plot, _prefix, _position, _exclude

def _filter_term(_term):

    _term = '' if _term == '-' else _term
    return _term


class Grotrian:

    def __init__(self, _atom, _path=None, _conf_prefix=None, _scaleFunc=None, _scaleFunc_inv=None):
        r"""

        Parameters
        -----------

        _atom : AtomCls.Atom
            object of the atomic model

        _path : str
            a path to the *.Grotrian configuration file

        _conf_prefix : strfilepath_dict["Grotrian"]
            common configuration string of the inner shell

        """

        self.set_atom( _atom, _path=_path,
             _conf_prefix=_conf_prefix, _scaleFunc=_scaleFunc, _scaleFunc_inv=_scaleFunc_inv )

    def set_atom(self, _atom, _path, _conf_prefix=None, _scaleFunc=None, _scaleFunc_inv=None):
        r""" """

        self.atom = _atom
        #_path = _path if _path is not None else _atom.filepath_dict["Grotrian"]
        self._path = _path
        self._scaleFunc = _scaleFunc
        self._scaleFunc_inv = _scaleFunc_inv

        if _path is None:
            line_plot = {}
            self.prefix = ''
            self.position = {}
            _exclude = {}
        else:
            with open(_path, 'r') as file:
                _fLines = file.readlines()
            line_plot, self.prefix, self.position, _exclude  = read_Grotrian_(_fLines)

        if _conf_prefix is not None:
            self.prefix = _conf_prefix



        #---------------------------------------------------------------------
        # prepare structures for plotting
        #---------------------------------------------------------------------
        singlet, multiplet, Lset, self.group_dict = _prepare_dict_(_atom=_atom, _conf_prefix=self.prefix,
                                                  _scaleFunc=_scaleFunc, _exclude=_exclude)


        self.singlet = singlet
        self.multiplet = multiplet
        self.Lset = Lset
        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        # a dictionary to store level position of levels in the plot.
        #---------------------------------------------------------------------
        self.pos_level = {}
        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        #---------------------------------------------------------------------
        self.member_to_head = self._init_group_member2head()
        self.line_plot_defined = line_plot


        self.fig = None

    def make_fig(self, _fig=None, _axe=None, _figsize=(6,8), _dpi=120, _f=200, _removeSpline=[], _resetFig=True):
        r"""

        Parameters
        ----------

        _figsize : tuple of int
            size of the figure, default : (6,8)

        _dpi : int
            "dot per inch" (resolution) of the figure, default : 120

        _f : int
            enlarge factor to show fine structure explicitly, default : 200
        """

        Lset = self.Lset
        singlet = self.singlet
        multiplet = self.multiplet
        pos_level = self.pos_level

        #---------------------------------------------------------------------
        # config
        #---------------------------------------------------------------------
        _hw = 0.5                       # half width of a level in the plot
        _b = len(Lset["singlet"])       # bias in the x axis of multiplet
        _fontsize = 14                  # fontsize of labels
        _textsize = 10                  # fontsize of text of term
        _Jsize = 7                      # fontsize of text of J
        _st = 0.1                       # horizontal spacing between term and text
        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        # set figure and axe
        #---------------------------------------------------------------------

        if _fig is not None:
            fig = _fig
            self.fig = fig
        else:
            if _resetFig or self.fig is None:
                fig = plt.figure(figsize=_figsize, dpi=_dpi)
                self.fig = fig
            else:
                fig = self.fig

        if _axe is not None:
            plt.sca( _axe )

        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        # remove splines
        #---------------------------------------------------------------------
        if _axe is not None:
            Plotting.remove_spline_(_axe, pos=_removeSpline)
        else:
            Plotting.remove_spline_(*self.fig.get_axes(), pos=_removeSpline)


        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        # excluded levels
        #---------------------------------------------------------------------
        for _text_tuple, _val in self.group_dict.items():

            _x_mid = float(_val['x_group'])
            xs_level = _x_mid-_hw, _x_mid+_hw

            ys_level = _val['y_group']/_val['count'], _val['y_group']/_val['count']

            plt.plot(xs_level, ys_level, "-k", linewidth=1)

            # store level posiiton, (conf_origin, term, J)
            for _ctj in _val["ctj_list"]:
                pos_level[ _ctj ] = {"xs" : xs_level,'ys' : ys_level}

            # plot text
            x_text = xs_level[1] + _st
            y_text = ys_level[0]
            plt.text(x_text, y_text, "{} {}".format(_text_tuple[0],_filter_term(_text_tuple[1])), fontsize=_textsize, color="k")


        #---------------------------------------------------------------------
        # singlet
        #---------------------------------------------------------------------
        for k0, v0 in singlet.items():
            for k1, v1 in v0.items():
                _ctj_ = (v1[2],k0[1],k1)
                _idx = self.atom._ctj_table.Level.index(_ctj_)

                #-------------------------------------------------------------
                try:
                    _x_mid = self.position[ _ctj_ ]
                except KeyError:
                    _x_mid = v1[1]
                #-------------------------------------------------------------

                if self.atom.Level["isGround"][_idx]:
                    xs_level = _x_mid-_hw, _x_mid+_hw+2*len(Lset["singlet"])+len(Lset["multiplet"])
                else:
                    xs_level = _x_mid-_hw, _x_mid+_hw
                ys_level = v1[0], v1[0]
                plt.plot(xs_level, ys_level, "-k", linewidth=1)
                # store level posiiton, (conf_origin, term, J)
                pos_level[ (v1[2], k0[1], k1) ] = {}
                pos_level[ (v1[2], k0[1], k1) ]["xs"] = xs_level
                pos_level[ (v1[2], k0[1], k1) ]["ys"] = ys_level
                # plot text
                x_text = xs_level[1] + _st
                y_text = ys_level[0]
                plt.text(x_text, y_text, "{} {}".format(k0[0],_filter_term(k0[1])), fontsize=_textsize, color="k")

        # a vertical line separates singlet panel and multiplet panel
        #plt.axvline(x=_b+1, linestyle="--", linewidth=0.5, color="k")
        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        # plotting multiplet
        #---------------------------------------------------------------------
        _hw = 0.6                    # multiplet need wider
        _sf = 0.1                    # space for connecting fine structure
        for k0, v0 in multiplet.items():
            # compute mean term energy
            y_count = 0.
            y_mean = 0
            for k1, v1 in v0.items():
                y_mean += v1[0]
                y_count += 1
            y_mean /= y_count

            for k1, v1 in v0.items():
                # plot level
                y_pos = y_mean + _f * (v1[0]-y_mean)            # enlarge space
                ys_level = y_pos, y_pos
                #-------------------------------------------------------------
                try:
                    _x_mid = self.position[ _ctj_ ]
                except KeyError:
                    _x_mid = v1[1]
                #-------------------------------------------------------------
                xs_level = _x_mid+1-_hw+_b, _x_mid+1+_hw+_b-_sf
                plt.plot(xs_level, ys_level, "-k", linewidth=1)
                # store level posiiton, (conf_origin, term, J)
                pos_level[ (v1[2], k0[1], k1) ] = {}
                pos_level[ (v1[2], k0[1], k1) ]["xs"] = xs_level
                pos_level[ (v1[2], k0[1], k1) ]["ys"] = ys_level
                # connect fine structure
                plt.plot([xs_level[1],xs_level[1]+_sf], [y_pos,y_mean], "--k", linewidth=0.5)
                # plot text of J
                plt.text(xs_level[0]-2*_st, y_pos, k1, fontsize=_Jsize, color="k")
            # plot text of term
            x_pos_text = xs_level[1]+_sf + _st
            plt.text(x_pos_text, y_mean, "{} {}".format(k0[0],_filter_term(k0[1])), fontsize=_textsize, color="k")
        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        # config of the plot
        #---------------------------------------------------------------------

        #--- tick of singlet
        xtick1 = list(Lset["singlet"])
        xticklabel1 = [CST.L_i2s_[xt] for xt in xtick1]
        #--- tick of multiplet
        xtick2 = [ x+_b+1 for x in list(Lset["multiplet"]) ]
        xticklabel2 = [CST.L_i2s_[xt-_b-1] for xt in xtick2]
        #--- customize xticks and xticklabels
        plt.xticks( xtick1 + xtick2, xticklabel1 + xticklabel2, fontsize=_fontsize )

        #--- x,y label; title
        if self.position is None:
            plt.xlabel("L", fontsize=_fontsize)
        plt.ylabel("E [eV]", fontsize=_fontsize, rotation=90)
        #plt.title(self.atom.Title, fontsize=_fontsize, y=1)

        #--- change x limit
        if len(xtick2) > 0:
            plt.xlim(-1, xtick2[-1]+2)

        #--- change y scale
        #if _forward is not None and _inverse is not None:
        ax = plt.gca()
        #def _forward(x):
        #    return x
        #def _inverse(x):
        #    return x
        #ax.set_yscale('function', functions=(_forward, _inverse))
        #ax.set_yscale("logit")

        #--- y ticklabel
        if self._scaleFunc is not None and self._scaleFunc_inv is not None:
            yticks = ax.get_yticks()
            yticks = [tick for tick in yticks if 1>=tick >=0]
            ax.set_yticks(yticks)
            ax.set_yticklabels( [ f"{self.atom.Level['erg'][:].max()/CST.eV2erg_*self._scaleFunc_inv(tick):.2f}" for tick in yticks  ] )
        #---------------------------------------------------------------------

        ylim = self.fig.gca().get_ylim()
        xlim = self.fig.gca().get_xlim()
        # ylim_max - ylim_min
        self.yr = ylim[1] - ylim[0]
        # xlim_max - xlim_min
        self.xr = xlim[1] - xlim[0]

        #---------------------------------------------------------------------
        # init line plot after self.pos_level has been defined
        #---------------------------------------------------------------------
        self.line_plot = self.init_line_plot( self.line_plot_defined )

    def save_fig(self, _filename, _dpi=120):
        r"""

        Parameters
        ----------

        _filename : str
            filename (including path) to save your figure

        _dpi : int
            "dot per inch" (resolution) to save figure, default : 120s
        """

        self.fig.savefig(_filename, dpi=_dpi)

    def show_fig(self):
        pass
        #plt.show()

    def _init_group_member2head(self):
        r""" """
        member_to_head = {}
        for _text_tuple, _val in self.group_dict.items():

            ## ignore group with only 1 level
            if len(_val["ctj_list"]) < 2:
                continue

            _head = _val["ctj_list"][0]
            for _ctj in _val["ctj_list"][1:]:
                member_to_head[_ctj] = _head
        return member_to_head


    def init_line_plot(self, line_plot_defined):
        r""" """

        line_plot = {}
        self.ctj_to_ctj_plot = {}

        for _ctj_ij in self.atom._ctj_table.Line + self.atom._ctj_table.Cont:

            ctj_i, ctj_j = _ctj_ij

            if ctj_i in self.member_to_head.keys():

                ## member -> member, ignore
                if ctj_j in self.member_to_head.keys():
                    continue
                ## member -> other  : head -> other
                else:
                    ctj_ij = ( self.member_to_head[ctj_i], ctj_j )

            ## other -> member : other -> head
            elif ctj_j in self.member_to_head.keys():

                ctj_ij = ( ctj_i, self.member_to_head[ctj_j] )

            ## other -> other
            else:
                ctj_ij = _ctj_ij

            ## a ctj map will be used in rate plot
            self.ctj_to_ctj_plot[ _ctj_ij ] = ctj_ij

            ## if defined in .Grotrian file, overwrite the parameters
            try:
                line_plot[ctj_ij] = line_plot_defined[ctj_ij]
            except KeyError:
                ## if not defined and already has parameter, ignore
                if ctj_ij not in line_plot.keys():
                    #r1, r2 = 0.5, 0.5
                    isGroud_lower = self.atom.Level["isGround"][ self.atom._ctj_table.Level.index(ctj_i) ]
                    isGroud_upper = self.atom.Level["isGround"][ self.atom._ctj_table.Level.index(ctj_j) ]
                    if isGroud_lower and isGroud_upper:
                        r1, r2 =  0.9, 0.9
                    else:
                        r1, r2 = init_vertical_r1r2_(self.pos_level[ctj_i]['xs'], self.pos_level[ctj_j]['xs'])
                    #print(ctj_ij, '  ', r1, r2)
                    line_plot[ctj_ij] = ["None", r1, r2]

        return line_plot

    def plot_transitions(self, _text_selection="wavelength", _hasText=True):
        r"""
        """

        _line_plot = self.line_plot


        for _ctj_ij in _line_plot.keys():
            _ctj1, _ctj2 = _ctj_ij
            wl, _r1, _r2 = _line_plot[_ctj_ij]


            if _text_selection == "all":
                if not _hasText:
                    wl = "None"
                self.connect_line(_cfj1=_ctj1, _cfj2=_ctj2, _r1=_r1, _r2=_r2, _c="black", _text=wl, _tsize=7, _r=0.4)

            elif _text_selection == "wavelength":
                if wl != "None":
                    if not _hasText:
                        wl = "None"
                    self.connect_line(_cfj1=_ctj1, _cfj2=_ctj2, _r1=_r1, _r2=_r2, _c="black", _text=wl, _tsize=7, _r=0.4)

            else:
                if wl in _text_selection:
                    if not _hasText:
                        wl = "None"
                    self.connect_line(_cfj1=_ctj1, _cfj2=_ctj2, _r1=_r1, _r2=_r2, _c="black", _text=wl, _tsize=7, _r=0.4)

        self.show_fig()

    def connect_line(self, _cfj1, _cfj2, _r1, _r2, _c, _text, _tsize=7, _r=0.4, _lcolor="black", _lwidth=0.7, _lstyle='-', _tcolor="black"):
        r"""

        Parameters
        ----------

        _cfj1 : tuple of str, (conf,term,J)
            to specify level 1

        _cfj2 : tuple of str, (conf,term,J)
            to specify level 2

        _r1 : float, in the range of [0 : 1]
            relative position along the line of level 1 starting from left point

        _r2 : float, in the range of [0 : 1]
            relative position along the line of level 2 starting from left point

        _c : str
            color of the line you're going to draw

        _text : str
            content of text

        _tsize : int
            fontsize of text, default : 7

        _r : float, in the range of [0 : 1]
            relative position along the line starting from left point, default : 0.4

        _lcolor : str
            color of the line, default : "black"

        _lwidth : float
            line width, default : 0.7

        _lstyle : str
            line style, default : '-'

        _tcolor : str
            color of the text, default : "black"
        """

        if _text == "None":

            _text = None

        pos_lvl1 = self.pos_level[_cfj1]
        pos_lvl2 = self.pos_level[_cfj2]

        # initialize left/right point of the line
        _lp = pos_lvl1["xs"][0] + (pos_lvl1["xs"][1] - pos_lvl1["xs"][0])*_r1 , pos_lvl1["ys"][0]
        _rp = pos_lvl2["xs"][0] + (pos_lvl2["xs"][1] - pos_lvl2["xs"][0])*_r2 , pos_lvl2["ys"][0]
        # if the initialization is not correct, swap them
        if _lp[0] > _rp[0] :
            _lp, _rp = _rp, _lp

        # compute the angle to rotate text
        _w, _h = self.fig.get_size_inches()
        _dx = (_rp[0]-_lp[0]) / self.xr * _w
        _dy = (_rp[1]-_lp[1]) / self.yr * _h
        _tangle = math.atan2( _dy, _dx )
        _tangle = math.degrees( _tangle )
        # plot line and text
        _ax = self.fig.gca()
        _line_obj, _text_obj = line_with_text_(_ax=_ax, _lp=_lp, _rp=_rp, _text=_text, _tsize=_tsize, _r=0.5,
                            _tangle=_tangle, _lcolor=_lcolor, _lwidth=_lwidth, _lstyle=_lstyle, _tcolor=_tcolor)



    def connect_arrow(self, _cfj_i, _cfj_j, _ri, _rj, _direction, _cmap, _norm, _v, _abserr, _asize=5, _lwidth=2):
        r"""
        """

        pos_lvl_i = self.pos_level[_cfj_i]
        pos_lvl_j = self.pos_level[_cfj_j]

        _pi = pos_lvl_i["xs"][0] + (pos_lvl_i["xs"][1] - pos_lvl_i["xs"][0])*_ri , pos_lvl_i["ys"][0]
        _pj = pos_lvl_j["xs"][0] + (pos_lvl_j["xs"][1] - pos_lvl_j["xs"][0])*_rj , pos_lvl_j["ys"][0]

        _ax = self.fig.gca()
        _annotation_obj = arrow_without_text_(_ax=_ax, _pi=_pi, _pj=_pj,
                            _direction=_direction, _cmap=_cmap, _norm=_norm,
                            _v=_v, _abserr=_abserr, _asize=_asize, _lwidth=_lwidth)

        return _annotation_obj

    def plot_transition_rate(self, _idxI, _idxJ, _rate, _direction, _cmap, _norm, _abserr=1E-5, _asize=5, _lwidth=2, _level_ctj_without_prefix=None):
        r"""
        """

        assert _idxI.shape[0] == _idxJ.shape[0] == _rate.shape[0]

        _line_plot = self.line_plot

        #---------------------------------------------------------------------
        # preprocess
        # (ctj_i,ctj_j) --> rate
        #---------------------------------------------------------------------
        ctj_rate_dict = {}
        for k in range(_idxI.shape[0]):
            _ctj_i = self.atom._ctj_table.Level[ _idxI[k] ]
            _ctj_j = self.atom._ctj_table.Level[ _idxJ[k] ]
            _ctj_ij = ( _ctj_i, _ctj_j )

            try:
                _ctj_ij_plot = self.ctj_to_ctj_plot[ _ctj_ij ]
            ## member/head -> member/head are not in self.ctj_to_ctj_plot
            except KeyError:
                continue

            try:
                ctj_rate_dict[_ctj_ij_plot] += _rate[k]
            except KeyError:
                ctj_rate_dict[_ctj_ij_plot]  = _rate[k]

        del _ctj_i, _ctj_j, _ctj_ij
        #---------------------------------------------------------------------
        # loop over ctj_rate_dict and plot arrow
        #---------------------------------------------------------------------
        #for key, value in _line_plot.items():
        #    print(key, '  ', value)

        for _ctj_ij_plot in ctj_rate_dict.keys():
            #print(_ctj_ij_plot)

            ## single level related
            if _level_ctj_without_prefix is not None:
                _level_ctj = ( self.prefix+_level_ctj_without_prefix[0],_level_ctj_without_prefix[1],_level_ctj_without_prefix[2] )
                if _level_ctj not in _ctj_ij_plot:
                    continue

            _, _ri, _rj = _line_plot[ _ctj_ij_plot ]
            _v       = ctj_rate_dict[ _ctj_ij_plot ]
            _ctj_i, _ctj_j = _ctj_ij_plot

            _annotation_obj = self.connect_arrow(_ctj_i, _ctj_j, _ri, _rj, _direction,
                                        _cmap, _norm, _v, _abserr=_abserr, _asize=_asize, _lwidth=_lwidth)

        # colorbar ax
        #_temp = [[1,1]]
        #_temp_ax = self.fig.add_axes([0.875, 0.2, 0.001, 0.001])
        #_img = _temp_ax.imshow(_temp, cmap=_cmap, norm=_norm)
        #self._temp_img = _img
        #_temp_ax.set_visible(False)
        #_cax = self.fig.add_axes([0.84, 0.15, 0.02, 0.7]) if _cax is None else _cax
        #self.fig.colorbar( _img, cax=_cax, orientation='vertical')

        self.show_fig()

    def add_colorbar(self, _cmap, _norm, _cax=None):

        # colorbar ax
        _temp = [[1,1]]
        _temp_ax = self.fig.add_axes([0.9, 9, 0.001, 0.001])
        _img = _temp_ax.imshow(_temp, cmap=_cmap, norm=_norm)
        _temp_ax.set_visible(False)
        if _cax is None:
            _cax = self.fig.add_axes([0.84, 0.15, 0.02, 0.7])
        self.fig.colorbar( _img, cax=_cax, orientation='vertical')




def init_vertical_r1r2_(xs_lower, xs_upper):

    #print(xs_lower, xs_upper)

    l_upper = xs_upper[1] - xs_upper[0]
    l_lower = xs_lower[1] - xs_lower[0]

    # no overlapping
    #         ----    or    ----
    #  ----                       ----
    #
    if xs_lower[1] <= xs_upper[0] or xs_lower[0] >= xs_upper[1]:
        return 0.5, 0.5

    # containing
    #     --------
    #       ----
    elif xs_lower[0] >= xs_upper[0] and xs_lower[1] <= xs_upper[1]:
        k = 0.3
        x = l_lower * k  + xs_lower[0]
        return k, (x-xs_upper[0])/l_upper
    # containing
    #       ----
    #      -------
    elif xs_lower[0] <= xs_upper[0] and xs_lower[1] >= xs_upper[1]:
        k = 0.7
        x = l_upper * k  + xs_upper[0]
        return (x-xs_lower[0])/l_lower,  k
    # partially overlap
    #       ------
    #     ----
    elif xs_lower[0] < xs_upper[0] < xs_lower[1]:
        x = 0.5* (xs_lower[1] + xs_upper[0])
        return (x-xs_lower[0])/l_lower, (x-xs_upper[0])/l_upper
    # partially overlap
    #       ------
    #           ----
    elif xs_lower[0] < xs_upper[1] < xs_lower[1]:
        x = 0.5 * (xs_upper[1] + xs_lower[0])
        return (x-xs_lower[0])/l_lower, (x-xs_upper[0])/l_upper


    else:
        print("here")
        return 0.5, 0.5
