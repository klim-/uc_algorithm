# -*- coding: utf-8 -*-

# read system from file: -----------------------------------------------

import sympy as sp
diff_symbols = sp.Matrix([])


# Placeholder for original (implicit) equations 
F_eq_orig = None

# this might store some substitutions to replace abreviations
# which have been introduced to simplify the model
subs_tuples = [] 

sys_name = None
data = None




#from ex.franke_519 import * # **** einradfahrer (kompliziert!!)
#from ex.mechanik_generisch_4fg_2sg import *
#from ex.system_generisch_4zust_2sg import *
#from ex.martin_gegenbeispiel import *
#from ex.franke_513 import *  # einfaches Beispiel zum Testen


#from ex.fahrzeug_odeometrie_beide_raeder import *
#from ex.mechanik_generisch_4fg_2sg import *
#from ex.mechanik_2pendel import *
#from ex.mechanik_wagen_pendel import *
#from ex.mechanik_acrobot import *
#from ex.inertia_wheel_pendulum import *
#from ex.mechanik_inertia_wheel_symb import *
#from ex.mechanik_tora_symb import *
#from ex.mechanik_einfachpendel_np2_nq2 import *


#from ex.mechanik_satellit import *


#from ex.system_generisch_4zust_2sg_v2 import *

from ex.mechanik_ph_wagen_pendel import *

