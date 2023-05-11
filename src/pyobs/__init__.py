"""
   Python interface to several Level 2 observing systems, mostly from EOS.
This is a port to Python 3 from Python 2 sources. In doing the portion,

Packages that have been ignored
-------------------------------

   lidar_l2.py    (interface to Pete's older lidar sampling)
   oracles.py     (
   dragon.py
   fpl.py
   g5_icartt.py
   modis.py
   sknet.py

Packages that have been renamed
-------------------------------

   calipso.py       renamed  calipso_l1p5.py
   calipso_lev2.py  renamed  calipso_l2.py

"""

from .aeronet   import AERONET_L2
from .aura      import AURA_L2
from .aqs       import AQS
from .avhrr     import AVHRR_L2B
from .bits      import BITS

from .calipso_l1p5 import CALIPSO_L1p5
from .calipso_l2   import CALIPSO_L2

from .config    import Config
from .cpl       import CPL_L2
from .dial      import DIAL
from .gcs03     import GCS03
from .geo04     import GEO04_L2
from .goci      import GOCI
from .goesx     import GOESX_CLRSKY
from .hsrl      import HSRL
from .hsrl_pbl  import HSRL_PBL
from .icartt    import ICARTT
# from .igbp      import IGBP # requires f2py extension
from .improve   import IMPROVE
from .kde       import *
# from .mapss     import MAPSS # needs igbp

from .man       import MAN
from .mcd43     import McD43
from .mcd43gf   import MCD43GF # inconsistent class name
from .minx      import MINX
from .mxd03     import MxD03
from .mxd04     import MxD04_L2
from .mxd06     import MxD06_L2
from .mxd14     import MxD14_L2

from .naaps     import NAAPS
from .nc4ctl    import Dataset as ctlDataset   
from .npz       import NPZ

from .odsreader import ODSreader
from .omaeruv   import OMAERUV_L2
from .omso2     import OMSO2_L2
from .omno2     import OMNO2_L2
from .omps      import OMPS_L2
from .omps_ai   import OMPS_AI

from .sev03     import SEV03 # inconsistent class name; Level 2?
#from .sgp4      import getTrack, dayPeriod # requires f2py extension

from .toms      import TOMS_L2





