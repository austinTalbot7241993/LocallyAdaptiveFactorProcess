import ctypes
import os
import sys
import numpy as np

def _find_liblafp():
    """Locate liblafp shared library in build directory or package paths."""
    lib_names = ["liblafp.dylib", "liblafp.so", "lafp.dll"]
    
    # Search relative to workspace build directories
    possible_dirs = [
        os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "build")),
        os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "build_asan")),
        os.path.abspath(os.path.join(os.path.dirname(__file__))),
        "/usr/local/lib",
        "/opt/homebrew/lib",
    ]
    
    for d in possible_dirs:
        for name in lib_names:
            path = os.path.join(d, name)
            if os.path.exists(path):
                return path
                
    # Fallback to system loader
    for name in lib_names:
        try:
            return ctypes.CDLL(name)._name
        except OSError:
            pass
            
    raise FileNotFoundError("Could not locate liblafp shared library (liblafp.dylib/liblafp.so). Please build the project with CMake first.")

class LAFPBindings:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(LAFPBindings, cls).__new__(cls)
            cls._instance._init_ctypes()
        return cls._instance

    def _init_ctypes(self):
        lib_path = _find_liblafp()
        self.lib = ctypes.CDLL(lib_path)

        # Function signatures
        self.lib.NGPmcmc_New.restype = ctypes.c_void_p
        self.lib.NGPmcmc_New.argtypes = []

        self.lib.NGPmcmc_init.restype = ctypes.c_int
        self.lib.NGPmcmc_init.argtypes = [ctypes.c_void_p]

        self.lib.NGPmcmc_free.restype = ctypes.c_int
        self.lib.NGPmcmc_free.argtypes = [ctypes.c_void_p]

        self.lib.NGPmcmc_construct.restype = ctypes.c_int
        self.lib.NGPmcmc_construct.argtypes = [
            ctypes.c_void_p,       # self
            ctypes.c_void_p,       # y (gsl_matrix*)
            ctypes.c_void_p,       # tobs (gsl_vector*)
            ctypes.c_int,          # Niter
            ctypes.c_void_p,       # sigU (gsl_vector*)
            ctypes.c_void_p,       # sigA (gsl_vector*)
            ctypes.c_double,       # sigEps
            ctypes.c_double,       # sigMu
            ctypes.c_double,       # sigAlph
            ctypes.c_double,       # a
            ctypes.c_double,       # b
            ctypes.c_void_p,       # th (marray3d*)
            ctypes.c_void_p        # sig (gsl_matrix*)
        ]

        self.lib.NGPmcmc_operations.restype = ctypes.c_int
        self.lib.NGPmcmc_operations.argtypes = [ctypes.c_void_p]

        # GSL Helpers
        self.lib.gsl_matrix_alloc.restype = ctypes.c_void_p
        self.lib.gsl_matrix_alloc.argtypes = [ctypes.c_size_t, ctypes.c_size_t]

        self.lib.gsl_matrix_free.restype = None
        self.lib.gsl_matrix_free.argtypes = [ctypes.c_void_p]

        self.lib.gsl_matrix_set.restype = None
        self.lib.gsl_matrix_set.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_double]

        self.lib.gsl_matrix_get.restype = ctypes.c_double
        self.lib.gsl_matrix_get.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t]

        self.lib.gsl_vector_alloc.restype = ctypes.c_void_p
        self.lib.gsl_vector_alloc.argtypes = [ctypes.c_size_t]

        self.lib.gsl_vector_free.restype = None
        self.lib.gsl_vector_free.argtypes = [ctypes.c_void_p]

        self.lib.gsl_vector_set.restype = None
        self.lib.gsl_vector_set.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_double]

        self.lib.gsl_vector_get.restype = ctypes.c_double
        self.lib.gsl_vector_get.argtypes = [ctypes.c_void_p, ctypes.c_size_t]

        # marray3d Helpers
        self.lib.marray3d_alloc.restype = ctypes.c_void_p
        self.lib.marray3d_alloc.argtypes = [ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t]

        self.lib.marray3d_free.restype = None
        self.lib.marray3d_free.argtypes = [ctypes.c_void_p]

        self.lib.marray3d_get.restype = ctypes.c_double
        self.lib.marray3d_get.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t]

    def numpy_to_gsl_matrix(self, arr):
        """Convert a 2D numpy array to a GSL matrix pointer."""
        arr = np.asarray(arr, dtype=np.float64)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        rows, cols = arr.shape
        gsl_m = self.lib.gsl_matrix_alloc(rows, cols)
        for i in range(rows):
            for j in range(cols):
                self.lib.gsl_matrix_set(gsl_m, i, j, float(arr[i, j]))
        return gsl_m

    def numpy_to_gsl_vector(self, arr):
        """Convert a 1D numpy array to a GSL vector pointer."""
        arr = np.asarray(arr, dtype=np.float64).ravel()
        size = len(arr)
        gsl_v = self.lib.gsl_vector_alloc(size)
        for i in range(size):
            self.lib.gsl_vector_set(gsl_v, i, float(arr[i]))
        return gsl_v

    def gsl_matrix_to_numpy(self, gsl_m, rows, cols):
        """Convert a GSL matrix pointer to a 2D numpy array."""
        arr = np.zeros((rows, cols), dtype=np.float64)
        for i in range(rows):
            for j in range(cols):
                arr[i, j] = self.lib.gsl_matrix_get(gsl_m, i, j)
        return arr

    def marray3d_to_numpy(self, marray, d1, d2, d3):
        """Convert a 3D marray3d pointer to a 3D numpy array."""
        arr = np.zeros((d1, d2, d3), dtype=np.float64)
        for i in range(d1):
            for j in range(d2):
                for k in range(d3):
                    arr[i, j, k] = self.lib.marray3d_get(marray, i, j, k)
        return arr
