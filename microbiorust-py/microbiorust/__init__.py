# microbiorust/__init__.py
from .microbiorust import *

__all__ = list(globals().keys())

try:
    from . import microbiorust as _base
    
    #list of submodules importable
    submodules = ["gbk", "embl", "align", "seqmetrics", "blast"]

    for sub_name in submodules:
        #get the submodule from the base binary
        sub_module = getattr(_base, sub_name, None)
        
        if sub_module:
            #making the submodule available: microbiorust.gbk
            globals()[sub_name] = sub_module
            __all__.append(sub_name)
            
            #make functions available at top level: microbiorust.gbk_to_faa()
            for func_name in dir(sub_module):
                if not func_name.startswith('_'):
                    globals()[func_name] = getattr(sub_module, func_name)
                    if func_name not in __all__:
                        __all__.append(func_name)

except (ImportError, AttributeError):
    # This ensures the package still "exists" even if the binary 
    # hasn't been compiled yet (useful during certain build steps)
    pass
