
#-------------------------------------------------------------------------------
# function definition to "help" a python object
#-------------------------------------------------------------------------------
# VERSION
#
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

import numpy as _numpy
import enum as _enum

_LENGTH_BAR    = 90
_LENGTH_NAME   = 25
_LENGTH_TYPE   = 35
_LENGTH_VLS    = 15
_LENGTH_SPACE  = 2
_STR_PRINT_MAX_LENGTH = 10

_SEP_LINE = "-" * _LENGTH_BAR + '\n'


def format_type_tuple_list_(vtype_name : T_STR, value : T_UNION[T_LIST,T_TUPLE]) -> T_STR:

    if len(value) == 0:
        return f"{vtype_name} of empty"

    is_single_type = True

    type0 = type( value[0] )
    for val in value[1:]:
        if not isinstance(val, type0):
            is_single_type = False
            break
    
    if is_single_type:
        return f"{vtype_name} of {type0.__name__}"
    else:
        return f"{vtype_name} of Any"

def format_type_dict_(vtype_name : T_STR, value : T_DICT) -> T_STR:
    
    if len(value) == 0:
        return f"{vtype_name} of empty"

    is_single_type = True

    for i, ( key, val ) in value.items():
        if i == 0:
            kt0 = type( key )
            vt0 = type( val )
            continue
        
        if not isinstance( key, kt0 ):
            is_single_type = False
            break
        if not isinstance( val, kt0 ):
            is_single_type = False
            break

    if is_single_type:
        return f"{vtype_name}:{kt0}->{vt0}"
    else:
        return f"{vtype_name}:Any->Any"


def _prefix_( level : T_INT ) -> T_STR :

    return ' '*2*level + "|- "

def format_struct_array_(arr : T_ARRAY, level : T_INT) -> T_STR:

    prefix = _prefix_( level )
    dtype = arr.dtype
    s = ''
    for field_name in dtype.names:
        s += f"{prefix+' '+field_name:{_LENGTH_NAME}s}{' '*_LENGTH_SPACE}"
        s += f"{arr[field_name].dtype.name:{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}"
        s += f"s: {arr[field_name].shape}\n"
    
    return s


def _help_attribute_( obj : T_ANY, level : T_INT ):

    prefix = _prefix_( level )

    s = ''
    if level == 0:
        s += _SEP_LINE
        s += f"{'name':<{_LENGTH_NAME}s}{' '*_LENGTH_SPACE}{'type':<{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}{'value/len/shape':{_LENGTH_VLS}s}\n"
        s += _SEP_LINE
    s += f"{' '*(len(prefix)-3)}{ type(obj).__name__ }\n"
    print(s, end='')

    attrs = vars( obj )

    for name, value in attrs.items():

        vtype = type( value )

        # ignore double under variables
        if len(name) > 2 and name[:2] == '__':
            continue

        s = f"{prefix+name:{_LENGTH_NAME}s}{' '*_LENGTH_SPACE}"#"{format_type_(vtype.__name__, value):{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}"
        
        #if isinstance( type(value), _enum.EnumMeta ):
        if isinstance( value, _enum.IntEnum ):

            s += f"{vtype.__name__:{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}"
            s += f"v: {value.name}\n"

        elif isinstance( value, int ) or \
           isinstance( value, float ) or \
           isinstance( value, complex ) or \
           isinstance( value,  bool) :
            
            s += f"{vtype.__name__:{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}"
            s += f"v: {value}\n"
        
        elif isinstance( value, str ) :

            s += f"{vtype.__name__:{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}"
            if len(s) > _STR_PRINT_MAX_LENGTH:
                s += f"    v: {value[:10]+'...'}\n"
            else:
                s += f"    v: {value}\n"

        elif isinstance( value, tuple ) or isinstance( value, list ):

            s += f"{format_type_tuple_list_(vtype.__name__, value):{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}"
            s += f"l: {len(value)}\n"

        elif isinstance( value, dict ):

            s += f"{format_type_dict_(vtype.__name__, value):{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}"
            s += f"l: {len(value)}\n"

        elif isinstance( value, _numpy.ndarray ):

            # normal array
            if value.dtype.fields is None:
                s += f"{'ndarray':{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}"
                s += f"s: {value.shape}\n"
            # structured array
            else:                    
                s += f"{'struct array':{_LENGTH_TYPE}s}{' '*_LENGTH_SPACE}"
                s += f"s: {value.shape}\n"
                # print arributes
                s += format_struct_array_(value, level+1)

        # class/struct
        else:
            print(s)
            _help_attribute_( value , level+1 )
            continue


        print(s, end='')


def help_( obj : T_ANY ):

    _help_attribute_( obj, 0 )