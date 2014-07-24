import numpy

# checks whether a string is a number
#@profile
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

def frexp10(x):
    exp = int(math.floor(math.log10(abs(x))))
    return x / 10**exp, exp
    
def trilinear(input_array,x,y,z):
    import numpy as np
    
    # check whether the three coordinate arrays are of the same shape:
    if x.shape != y.shape :
        print 'Error in trilinear'
        raise ValueError
    if y.shape != z.shape :
        print 'Error in trilinear'
        raise ValueError
    if z.shape != x.shape :
        print 'Error in trilinear'
        raise ValueError
    
    output = np.empty(x.shape)

    x0 = x.astype(np.integer)
    y0 = y.astype(np.integer)
    z0 = z.astype(np.integer)
    x1 = x0 + 1
    y1 = y0 + 1
    z1 = z0 + 1

    #Check if xyz1 is beyond array boundary:
    x1[np.where(x1==input_array.shape[0])] = x0.max()
    y1[np.where(y1==input_array.shape[1])] = y0.max()
    z1[np.where(z1==input_array.shape[2])] = z0.max()

    x = x - x0
    y = y - y0
    z = z - z0
    output = (input_array[x0,y0,z0]*(1-x)*(1-y)*(1-z) +
                 input_array[x1,y0,z0]*x*(1-y)*(1-z) +
                 input_array[x0,y1,z0]*(1-x)*y*(1-z) +
                 input_array[x0,y0,z1]*(1-x)*(1-y)*z +
                 input_array[x1,y0,z1]*x*(1-y)*z +
                 input_array[x0,y1,z1]*(1-x)*y*z +
                 input_array[x1,y1,z0]*x*y*(1-z) +
                 input_array[x1,y1,z1]*x*y*z)

    return output

def quadlinear(input_array,x,y,z,w):
    import numpy as np
    
    # check whether the three coordinate arrays are of the same shape:
    if x.shape != y.shape :
        print 'Error in quadlinear'
        raise ValueError
    if y.shape != z.shape :
        print 'Error in quadlinear'
        raise ValueError
    if z.shape != x.shape :
        print 'Error in quadlinear'
        raise ValueError
    if w.shape != x.shape :
        print 'Error in quadlinear'
        raise ValueError

    output = np.empty(x.shape)
    
    x0 = x.astype(np.integer)
    y0 = y.astype(np.integer)
    z0 = z.astype(np.integer)
    w0 = w.astype(np.integer)
    x1 = x0 + 1
    y1 = y0 + 1
    z1 = z0 + 1
    w1 = w0 + 1
    
    #Check if xyz1 is beyond array boundary:
    #x1[np.where(x1==input_array.shape[0])] = x0.max()
    #y1[np.where(y1==input_array.shape[1])] = y0.max()
    #z1[np.where(z1==input_array.shape[2])] = z0.max()
    #w1[np.where(w1==input_array.shape[3])] = w0.max()

    x = x - x0
    y = y - y0
    z = z - z0
    w = w - w0
    output = (input_array[x0,y0,z0,w0,:,:]*(1-x)*(1-y)*(1-z)*(1-w) +
              
              input_array[x1,y0,z0,w0,:,:]*x*(1-y)*(1-z)*(1-w) +
              input_array[x0,y1,z0,w0,:,:]*(1-x)*y*(1-z)*(1-w) +
              input_array[x0,y0,z1,w0,:,:]*(1-x)*(1-y)*z*(1-w) +
              input_array[x0,y0,z0,w1,:,:]*(1-x)*(1-y)*(1-z)*w +
              
              input_array[x0,y0,z1,w1,:,:]*(1-x)*(1-y)*z*w +
              input_array[x0,y1,z0,w1,:,:]*(1-x)*y*(1-z)*w +
              input_array[x0,y1,z1,w0,:,:]*(1-x)*y*z*(1-w) +
              input_array[x1,y0,z0,w1,:,:]*x*(1-y)*(1-z)*w +
              input_array[x1,y0,z1,w0,:,:]*x*(1-y)*z*(1-w) +
              input_array[x1,y1,z0,w0,:,:]*x*y*(1-z)*(1-w) +
              
              input_array[x0,y1,z1,w1,:,:]*(1-x)*y*z*w +
              input_array[x1,y0,z1,w1,:,:]*x*(1-y)*z*w +
              input_array[x1,y1,z0,w1,:,:]*x*y*(1-z)*w +
              input_array[x1,y1,z1,w0,:,:]*x*y*z*(1-w) +
              
              input_array[x1,y1,z1,w1,:,:]*x*y*z*w)
              
    return output


def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """
    import numpy as np
    
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out


def interp3(x, y, z, v, xi, yi, zi, **kwargs):
    """Sample a 3D array "v" with pixel corner locations at "x","y","z" at the
    points in "xi", "yi", "zi" using linear interpolation. Additional kwargs
    are passed on to ``scipy.ndimage.map_coordinates``."""
    import numpy as np
    from scipy.ndimage import map_coordinates
    
    if v.shape != (x.size,y.size,z.size):
        print 'size mismatch in interp3'
        raise ValueError
    
    def index_coords(corner_locs, interp_locs):
        index = np.arange(len(corner_locs))
        if np.all(np.diff(corner_locs) < 0):
            corner_locs, index = corner_locs[::-1], index[::-1]
        return np.interp(interp_locs, corner_locs, index)

    orig_shape = np.asarray(xi).shape
    xi, yi, zi = np.atleast_1d(xi, yi, zi)
    for arr in [xi, yi, zi]:
        arr.shape = -1

    output = np.empty(xi.shape, dtype=float)
    coords = [index_coords(*item) for item in zip([x, y, z], [xi, yi, zi])]

    map_coordinates(v, coords, order=3, output=output, **kwargs)

    return output.reshape(orig_shape)
    
# 5d and 6d interpolation
#@profile
def interp5(x, y, z, w, q, v, xi, yi, zi, wi, qi, **kwargs):
    """Sample a 5D array "v". Additional kwargs
    are passed on to ``scipy.ndimage.map_coordinates``."""
    import numpy as np
    from scipy.ndimage import map_coordinates
    
    if v.shape != (x.size,y.size,z.size, w.size, q.size):
        print 'size mismatch in interp3'
        raise ValueError
    
    def index_coords(corner_locs, interp_locs):
        index = np.arange(len(corner_locs))
        if np.all(np.diff(corner_locs) < 0):
            corner_locs, index = corner_locs[::-1], index[::-1]
        return np.interp(interp_locs, corner_locs, index)

    orig_shape = np.asarray(xi).shape
    xi, yi, zi, wi, qi = np.atleast_1d(xi, yi, zi, wi, qi)
    for arr in [xi, yi, zi, wi, qi]:
        arr.shape = -1

    output = np.empty(xi.shape, dtype=float)
    coords = [index_coords(*item) for item in zip([x, y, z, w, q], [xi, yi, zi, wi, qi])]

    map_coordinates(v, coords, order=3, output=output, **kwargs)

    return output.reshape(orig_shape)

#@profile
def interp6(x, y, z, w, q, qq, v, xi, yi, zi, wi, qi, qqi, **kwargs):
    """Sample a 5D array "v". Additional kwargs
    are passed on to ``scipy.ndimage.map_coordinates``."""
    import numpy as np
    from scipy.ndimage import map_coordinates
    
    if v.shape != (x.size,y.size,z.size, w.size, q.size, qq.size):
        print 'size mismatch in interp3'
        raise ValueError
    
    def index_coords(corner_locs, interp_locs):
        index = np.arange(len(corner_locs))
        if np.all(np.diff(corner_locs) < 0):
            corner_locs, index = corner_locs[::-1], index[::-1]
        return np.interp(interp_locs, corner_locs, index)

    orig_shape = np.asarray(xi).shape
    xi, yi, zi, wi, qi, qqi = np.atleast_1d(xi, yi, zi, wi, qi, qqi)
    for arr in [xi, yi, zi, wi, qi, qqi]:
        arr.shape = -1

    output = np.empty(xi.shape, dtype=float)
    coords = [index_coords(*item) for item in zip([x, y, z, w, q, qq], [xi, yi, zi, wi, qi, qqi])]

    map_coordinates(v, coords, order=3, output=output, **kwargs)

    return output.reshape(orig_shape)
