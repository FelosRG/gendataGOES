import os
import sys
PATH_SCRIPT = os.path.realpath(__file__) 
DIR_SCRIPT  = "/".join(PATH_SCRIPT.split("/")[:-1])
DIR_REPO    = "/".join(DIR_SCRIPT.split("/")[:-1])

import datetime
from configparser import ConfigParser


# VARIABLES
# ---------

productos = {
    "1":"ABI-L1b-RadC"  , 
    "2":"ABI-L1b-RadC"  ,
    "3":"ABI-L1b-RadC"  ,
    "4":"ABI-L1b-RadC"  ,
    "5":"ABI-L1b-RadC"  ,
    "6":"ABI-L1b-RadC"  ,
    "7":"ABI-L1b-RadC"  ,
    "8":"ABI-L1b-RadC"  ,
    "9":"ABI-L1b-RadC"  ,
    "10":"ABI-L1b-RadC" ,
    "11":"ABI-L1b-RadC" ,
    "12":"ABI-L1b-RadC" ,
    "13":"ABI-L1b-RadC" ,
    "14":"ABI-L1b-RadC" ,
    "15":"ABI-L1b-RadC" ,
    "16":"ABI-L1b-RadC" ,
    "CTH":"ABI-L2-ACHAC",
    "CM" :"ABI-L2-ACMC" ,
    "CAPE":"ABI-L2-DSIC",
    "COD":"ABI-L2-CODC" ,
}

variables = {
    "1":"Rad",
    "2":"Rad",
    "3":"Rad",
    "4":"Rad",
    "5":"Rad",
    "6":"Rad",
    "7":"Rad",
    "8":"Rad",
    "9":"Rad",
    "10":"Rad",
    "11":"Rad",
    "12":"Rad",
    "13":"Rad",
    "14":"Rad",
    "15":"Rad",
    "16":"Rad",
    "CTH":"HT" ,
    "CM" :"BCM",
    "CAPE":"CAPE",
    "COD":"COD",
}

def string2bool(string):
    if string in ["True","true","Verdadero","verdadero"]:
        return True
    elif string in ["False","false","Falso","falso"]:
        return False
    else:
        raise ValueError("El string ingresado no pudo ser intepretado como bool.")

try:
    config = ConfigParser()
    config.read_file(open(f"{DIR_REPO}/gendata.config","r"))
except FileNotFoundError:
    print("Error: No se encontró archivo de configuración.")
    sys.exit()


# SECCIÓN "GENERAL"
# -----------------
NOMBRE_DATASET = config["General"]["nombre_dataset"]

# SECCIÓN "DESCARGA"
#--------------------
FORMATO_FECHA = "%Y/%m/%d"
str_fecha_inicio = config["Descarga"]["fecha_inicio"]
str_fecha_final  = config["Descarga"]["fecha_final"]
FECHA_INICIO = datetime.datetime.strptime(str_fecha_inicio, FORMATO_FECHA).date()
FECHA_FINAL  = datetime.datetime.strptime(str_fecha_final , FORMATO_FECHA).date()

FORMATO_HORA = "%H:%M"
str_hora_inicio = config["Descarga"]["hora_inicio_utc"]
str_hora_final  = config["Descarga"]["hora_final_utc"]
hora_inicio = datetime.datetime.strptime(str_hora_inicio, FORMATO_HORA).time()
hora_final  = datetime.datetime.strptime(str_hora_final , FORMATO_HORA).time()
HORA_INICIO_UTC, MIN_INICIO_UTC = hora_inicio.hour, hora_inicio.minute
HORA_FINAL_UTC , MIN_FINAL_UTC  = hora_final.hour , hora_final.minute
SALTO_ARCHIVOS = int(config["Descarga"]["saltear_archivos"])

DIAS_POR_TOMAR = int(config["Descarga"]["num_dias"])

# SECCIÓN "BANDAS"
#--------------------
BANDAS = config["Bandas"]["bandas"].split(",")
VENTANA = int(config["Bandas"]["ventana"])

# SECCIÓN "RESOLUCION"
#--------------------
REESCALADO_KM = float(config["Resolucion Km"]["reescalado"])

# SECCIÓN "Ubicaciones"
#---------------------
INF_LAT = float(config["Ubicaciones"]["limite_inferior_latitud"])
SUP_LAT = float(config["Ubicaciones"]["limite_superior_latitud"])
INF_LON = float(config["Ubicaciones"]["limite_inferior_longitud"])
SUP_LON = float(config["Ubicaciones"]["limite_superior_longitud"])

RESOLUCION_GRID = int(config["Ubicaciones"]["resolucion_grid"]) 

if __name__ == "__main__":
    print(DIR_SCRIPT)
    print(DIR_REPO)