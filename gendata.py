"""
Modulo con funciones útiles para generar un dataset
con los datos de GOES16 (CONUS)
"""
import os
import sys
PATH_SCRIPT = os.path.realpath(__file__) 
DIR_SCRIPT  = "/".join(PATH_SCRIPT.split("/")[:-1])
sys.path.append(f"{DIR_SCRIPT}/lib/")
print(f"{DIR_SCRIPT}/lib/")
import h5py
import datetime
import numpy as np
import pandas as pd
from pathlib import Path
import config
import GOES

DIR_DESCARGAS = f"{DIR_SCRIPT}/Descargas"
DIR_DATASETS  = f"{DIR_SCRIPT}/Datasets"

# Creamos directorio si no existen.
Path(DIR_DATASETS).mkdir(parents=True,exist_ok=True)


def generar_lista_dias(fecha_inicio,fecha_final,num_dias):
    """
    Genera una lista con el número de días especificado
    dado el intervalo de fechas proporcionado.

    Se escogen los días de forma que queden de forma equidistante.

    Input:
        * fecha_inicio: datetime.datetime
        * fecha_final : datetime.datetime
        * num_dias : int

    Return:
        * lista_dias (list)
    
        la lista_dias contendrá una lista de int, que marcarán
        los días transcurridos desde el primer día del año de la fecha
        de inicio. La intención es que se use este dato con los métodos
        de datetime.timedelta
    """
    lista_dias     = list(pd.date_range(fecha_inicio,fecha_final,num_dias))
    primer_día_año = datetime.datetime(fecha_inicio.year,1,1,0,0)

    lista_num_dias = []
    for dia in lista_dias:
        diferencia = (dia - primer_día_año).days
        lista_num_dias.append(diferencia)
    return lista_num_dias

def descargar_bandas(identificador,datetime_inicio,datetime_final,saltos=1):
    """
    Descarga los datos satélitales del  GOES16 en la región de
    CONUS. (Contiguos United States). La descarga se hará de los
    productos indicados con su identificador. (Ver readmi).

    Input:
        * identificador: list o str, ver readmi.
        * datetime_inicio: datetime.datetime
        * datetime_final : datetime.datetime
        * saltos : int (default = 1)

        Los saltos es usado cuando no quieres descargar la totalidad
        de los datos satélitales en el intervalo indicado. Por ejemplo,
        en general el satélite genera un nuevo dato cada 5 minutos, pero
        se quieren descargar los datos cada 15 minutos, para ello 
        pasamos como parámetro saltos=3, así solo se descargarán
        un archivo de cada 3 disponibles.

    Output:
        * None
    """

    # Si es un unico identificador lo pasamos a lista por consistencia.
    if type(identificador) == str:
        identificador = [identificador]

    for ident in identificador:
        # Obtenemos los detalles de los productos a descargar.
        producto = config.productos[ident]
        banda    = None
        if producto == "ABI-L1b-RadC": banda = int(ident)

        # Iniciamos descarga.
        GOES.descargaIntervaloGOES16(
            producto=producto,
            datetime_inicio=datetime_inicio,
            datetime_final=datetime_final,
            banda=banda,
            output_path=f"{DIR_DESCARGAS}/{ident}/",
            verbose=False,
            saltos=saltos 
        )

def generar_ubicaciones(lat_inf,lat_sup,lon_inf,lon_sup,res_grid,res_goes):
    """
    Dado los límites de coordenadas especificadas se genera una lista de
    ubicaciones en la imágen dentro de ese cuadrante.

    Input:
        * lat_inf: float, Coordenada de latitud inferior del cuadrante.
        * lat_sup: float, Coordenada de latitud superior del cuadrante.
        * lon_inf: float
        * lon_sup: float
        * res_grid: int, Resolución del grid con el que se dividirá la región.
        * res_goes: float, Resolución espacial de las imágenes del GOES (en km)

    Return:
        * coordenadas_px, np.array

        Las coordenadas_px consiste en las coordenadas (en pixeles) de las 
        localizaciones generadas dentro de cuadrante especificado.
        coodenadas_px.shape : (num_localizaciones,2) con axis 1 como (px_x,px_y)
        la ubicación del punto en pixeles. (int,int)
    """

    x = np.linspace(lat_inf,lat_sup,res_grid)
    y = np.linspace(lon_inf,lon_sup,res_grid)
    X , Y = np.meshgrid(x,y)
    coordenadas = np.stack([X.reshape(-1),Y.reshape(-1)],axis=1)

    # Transformamos la lista en coordendas en pixeles
    coordenadas_px = []
    for lat,lon in coordenadas:
        px_x,px_y = GOES.coordenadas2px(nc=res_goes,latitud=lat,longitud=lon)
        coordenadas_px.append([px_x,px_y])
    coordenadas_px = np.array(coordenadas_px)
    return coordenadas_px

def reescalar(producto,reescaldo_km):
    """
    Realiza un reescalado a la resolución espacial establecida.

    Input:
        * producto: GOES.Producto
        * reescalado_km: int o float
    
    Output:
        * producto
    """
    resolucion = producto.resolución
    escalado   = resolucion/reescaldo_km
    if escalado < 1:
        downscale = int(1/escalado)
        producto.downscale(downscale)
    elif escalado > 1:
        producto.zoom(int(escalado))
    return producto

def generar_dataset(identificador,lista_ubicaciones,escala_km,ventana,nombre_dataset="dataset"):
    """
    Pipeline para generar un dataset en formato .h5, ideal para el entrenamiento
    de un red neuronal.

    Input:
        * identificador: list o str
        * lista_ubicaciones: list, lista generada con la función generar_ubicaciones.
        * escala_km: float o int, se reescalarán los datos a esa escala. Debe de ser igual a res_goes de la lista_ubicaciones.
        * ventana: int, tamaño de las imágenes que comprenderá el dataset.
        * nombre_dataset: str, Nombre del archivo que contendrá el dataset. No debe de colocarse la extención.
    
    Output:
        * None
    """

    # Nombre de salida.
    output_name = f"{DIR_DATASETS}/{nombre_dataset}.h5"

    # Mantenemos consitencia.
    if type(identificador) == str:
        identificador = [identificador]

    with h5py.File(output_name,"w") as file:
        for ident in identificador:

            # Obtenemos la lista de archivos disponibles.
            path_archivos  = f"{DIR_DESCARGAS}/{ident}/"
            lista_archivos = os.listdir(path_archivos)
            lista_archivos.sort() # Importante para respetar el orden cronológico.
            num_archivos   = len(lista_archivos)

            # Lista para guardar todos los elementos de una variable.
            lista_datos_variable = []

            for index_archivo in range(num_archivos):
                path_archivo = path_archivos + lista_archivos[index_archivo]
                producto     = GOES.Producto(path_archivo)

                # Mandamos los fill_values a 0
                fill_value = producto.fill_value
                producto.array[producto.array == fill_value] = 0
                producto = reescalar(producto,escala_km)
                
                lista_datos_archivo = []
                for index_coord in range(len(lista_ubicaciones)):
                    px_x,px_y     = lista_ubicaciones[index_coord,0] , lista_ubicaciones[index_coord,1]
                    ventana_array = GOES.obtener_ventana(producto.array,px_x,px_y,ventana=ventana)
                    lista_datos_archivo.append(ventana_array)                    

                # Recolectamos los datos de cada archivo.
                lista_datos_variable.append(lista_datos_archivo)
        
            # Creamos el dataset en el archivo .h5
            lista_datos_variable = np.array(lista_datos_variable)
            file.create_dataset(data=lista_datos_variable,name=ident,dtype=np.float32,compression="gzip")

            print(f"Dataset creado de la variable {ident} con shape {lista_datos_variable.shape}.")


if __name__ == "__main__":
    print(f"""
--------------------------------
INICIANDO GENERACIÓN DEL DATASET
--------------------------------

PASO 1/2: DESCARGANDO LAS DATOS DE LAS BANDAS.

""")

    # Generamos la lista de días.
    lista_dias = generar_lista_dias(config.FECHA_INICIO,config.FECHA_FINAL,config.DIAS_POR_TOMAR)
    
    # Obtenemos la fecha de referencia.
    año_inicio   = config.FECHA_INICIO.year
    fecha_ref_inicio = datetime.datetime(año_inicio,1,1,config.HORA_INICIO_UTC,config.MIN_INICIO_UTC)
    fecha_ref_final  = datetime.datetime(año_inicio,1,1,config.HORA_FINAL_UTC,config.MIN_FINAL_UTC)

    dias_descargados = 0
    for num_dia in lista_dias:
        
        # Obtenemos el intervalo del dia a descargar.
        fecha_inicio = fecha_ref_inicio + datetime.timedelta(days=num_dia)
        fecha_final  = fecha_ref_final + datetime.timedelta(days=num_dia)

        print(f"Descargando dia {dias_descargados+1} de {config.DIAS_POR_TOMAR}...")
        descargar_bandas(
            identificador=config.BANDAS,
            datetime_inicio=fecha_inicio,
            datetime_final=fecha_final,
            saltos=config.SALTO_ARCHIVOS
        )
    
    print("Descarga finalizada!")

    print("""
PASO 2/2: GENERANDO DATASET.

""")

    # Generamos ubicaciones.
    coord_px = generar_ubicaciones(
        config.INF_LAT,
        config.SUP_LAT,
        config.INF_LON,
        config.SUP_LON,
        res_grid=config.RESOLUCION_GRID,
        res_goes=config.REESCALADO_KM,
    )

    generar_dataset(
        identificador=config.BANDAS,
        lista_ubicaciones=coord_px,
        escala_km=config.REESCALADO_KM,
        ventana=config.VENTANA,
        nombre_dataset=config.NOMBRE_DATASET
    )








    
