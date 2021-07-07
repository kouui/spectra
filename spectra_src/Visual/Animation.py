
#-------------------------------------------------------------------------------
# definition of functions for animation
#-------------------------------------------------------------------------------
# VERSION
# 0.1.0
#    2021/07/04   u.k.   
# 0.1.0 
#    2021/05/18   u.k.   spectra-re
#-------------------------------------------------------------------------------

from ..ImportAll import *

import matplotlib.pyplot as plt
from matplotlib import animation
#from IPython.display import display, HTML



#def image_array_to_animation(image_array : T_ARRAY):
#    """return images sequence as an animation
#
#    Parameters
#    ----------
#    image_array : [type]
#        [description]
#
#    Returns
#    -------
#    [type]
#        [description]
#    """
#
#
#    ''' Display images sequence as an animation in jupyter notebook
#    
#    Args:
#        image_array(numpy.ndarray): image_array.shape equal to (num_images, height, width, num_channels)
#    '''
#
#
#    dpi = 72.0
#    xpixels, ypixels = image_array[0].shape[:2]
#    fig = plt.figure(figsize=(ypixels/dpi, xpixels/dpi), dpi=dpi)
#    im = plt.figimage(image_array[0])
#
#    def animate(i):
#        im.set_array(image_array[i])
#        return (im,)
#
#    anim = animation.FuncAnimation(fig, animate, frames=len(image_array), interval=33, repeat_delay=1, repeat=True)
#    #display(HTML(anim.to_html5_video()))