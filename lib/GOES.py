"""
Descripción:
Módulo con todo lo necesario para descargar
y trabajar con datos de los satélites GOES
párticularmente el satélite GOES-16.

Última modificación:
13 de Mayo del 2022

Autores/Fuentes:
Adrían Ramírez, Facultad de Ciencias, UNAM
felos@ciencias.unam.mx
FelosRG@github

Parte del código de este módulo fue tomado o 
modificado del proyecto GOES2GO:
https://github.com/blaylockbk/goes2go
"""

import os
_path_script = os.path.realpath(__file__) 
_path_script = "/".join(_path_script.split("/")[:-1])
_path_root = os.path.realpath(__file__)
_path_root = "/".join(_path_root.split("/")[:-2])

import os
import re
import s3fs
import math
import time
import netCDF4
import datetime

import numpy  as np
import pandas as pd
import scipy.ndimage
from matplotlib        import cm
from matplotlib.colors import ListedColormap

from pathlib import Path
import h5py


# Diccionario para dactor de reescalado.
# La resolución estandar son 2km. Por ejemplo, si tenemos banda de 0.5 km,
# hay que multipllicar la ventana por 4 para igualar la región de una banda
# de 2km.
reescalado = {
    2:1,
    1:2,
    0.5:4,
}

def degree2rad(degree):
    """
    Pasa de grados a radianes.
    """
    k   = math.pi / 180
    rad = degree*k
    return rad

def rad2degree(rad):
    """
    Pasa de radianes a grados.
    """
    k = 180 / math.pi
    degree = rad*k 
    return degree

def obtenerFecha(nc,return_datetime=False):
    """ 
    Devuelve un string con la fecha de la imágen.
    
    Parámetros:

    nc (netCDF4.Dataset ó int): Objeto de entrada o entero
        (int).

    return_datetime (bool): Si está en True, devuelve
        la fecha como objeto datetime. De lo contrario
        lo devuelve como string.
        ( Se deuelve en UTC, que es la zona horaria en 
        el que viene la información temporal de los 
        archivos nc )

        !! Falta hacer que devuelva la hora en formato de 24h

        !! Falta implementar que de la hora para otras zonas
        horarias.

    Es posible indicar en el parámetro nc el tiempo en unix que 
    será convertido en UTC.
    """
    if type(nc) == int: 
        t = nc
    else:
        # Obtiene los segundos desde 2000-01-01 12:00
        try:
            t = float(np.array(nc.variables["t"]))
        except:
            t = float(np.array(nc.variables["time"])[0])

    # Con ayuda de la librería "datetime" obtenemos la fecha actual.
    fecha_inicio = datetime.datetime(2000,1,1,12,0,0)
    time_delta   = datetime.timedelta(seconds=t)
    fecha        = fecha_inicio + time_delta
    if return_datetime:
        return fecha
    else:
        formato      = "%Y-%m-%d_%H-%M"
        return fecha.strftime(formato)

def obtenerBanda(nc):
    """ 
    Obtiene de que banda es el archivo nc correspondiente a el producto de
    "radiation".
    """
    id_banda = np.array(nc.variables["band_id"])
    id_banda = int(id_banda)
    return id_banda


def irradiancia2temperatura(array,nc):
    """
    Pasa de los valores de irradiancia a grados 
    centígrados para las bandas emisoras

    Las bandas emisoras son bandas como:
    6,7,8,9,10,11,12,13,14,15,16
    """
    fk1 = float(np.array(nc.variables["planck_fk1"]))
    fk2 = float(np.array(nc.variables["planck_fk2"]))
    bc1 = float(np.array(nc.variables["planck_bc1"]))
    bc2 = float(np.array(nc.variables["planck_bc2"]))
    a = np.log(1 + (fk1 / array))
    b = (fk2 / a - bc1)
    resultado  = b / bc2
    return resultado


def px2coordenadas(nc,px_x,px_y):
    """
    Pasa de pixeles en la imágen a coordenadas.
    """
    # Parámetros del satélite.
    # Fixed Grid scanning angles.
    X = nc.variables["x"]
    Y = nc.variables["y"]
    # Longitud of proyection of the origin
    lambda_o = nc.variables["goes_imager_projection"].longitude_of_projection_origin
    lambda_o = degree2rad(lambda_o)
    # Semi major axis value
    r_eq   = 6378137          # [m]
    # Semi minor axis value
    r_pool = 6356752.31414    # [m]
    # Satellite Hight from center of earth [m]
    H      = 42164160         # [m]
    
    # Cálculos previos.
    frac_r = r_eq / r_pool
    coef1  = frac_r**2
    x = X[px_x]
    y = Y[px_y]
    cosx = math.cos(x)
    cosy = math.cos(y)
    sinx = math.sin(x)
    siny = math.sin(y)
    a = sinx**2 + (cosx**2)*(cosy**2 + coef1*siny**2)
    b = -2*H*cosx*cosy
    c = H**2 - r_eq**2
    r_s = (-b-math.sqrt(b**2 - 4*a*c)) / 2*a
    s_x =  r_s*cosx*cosy
    s_y = -r_s*sinx
    s_z =  r_s*cosx*siny
    coef2 = s_z / math.sqrt((H-s_x)**2 + s_y**2)
    coef3 = s_y / (H-s_x)
    latitud  = math.atan(coef1*coef2)
    longitud = lambda_o  - math.atan(coef3)
    # Pasamos de rads a grados.
    latitud  = rad2degree(latitud )
    longitud = rad2degree(longitud)
    return latitud , longitud

def coordenadas2px(nc,latitud,longitud):
    """
    Pasa de coordenadas a localización en px.
    """
    # Parámetros del satélite.
    # Alternativa rápida a no tener que dar el nc.
    if nc==2:
        with h5py.File(_path_root + "/Recursos/CONUS/Georef_2km.h5") as dataset:
            X = dataset["x"][()]
            Y = dataset["y"][()]
            lambda_o = dataset["lambda_o"][()]
    else:
        try:
            X , Y , lambda_o = nc
        except TypeError:
            # Fixed Grid scanning angles.
            X = nc.variables["x"]
            Y = nc.variables["y"]
            # Longitud of proyection of the origin
            lambda_o = nc.variables["goes_imager_projection"].longitude_of_projection_origin
    lambda_o = degree2rad(lambda_o)
    # Semi major axis value
    r_eq   = 6378137          # [m]
    # Semi minor axis value
    r_pool = 6356752.31414    # [m]
    # Satellite Hight from center of earth [m]
    H      = 42164160         # [m]
    # exentricidad 
    e = 0.0818191910435
    # Pasamos de grados a radianes
    latitud  = degree2rad(latitud )
    longitud = degree2rad(longitud)
    # Cálculos intermedios
    coef1 = (r_pool / r_eq)**2
    phi_c = math.atan(coef1*math.tan(latitud))
    r_c   = r_pool / math.sqrt(1-(e*math.cos(phi_c))**2)
    s_x   = H - r_c*math.cos(phi_c)*math.cos(longitud-lambda_o)
    s_y   = -r_c*math.cos(phi_c)*math.sin(longitud -lambda_o)
    s_z   = r_c*math.sin(phi_c)
    # Revisamos la visiblidad desde el satélite.
    inequality1 = H*(H-s_x)
    inequality2 = s_y**2 + (s_z**2)*(r_eq/r_pool)**2
    message = f"Coordenada no visibles desde el satélite: {latitud},{longitud}"
    if inequality1 < inequality2:
        raise ValueError(message)
    
    # Obtenemos los ángulos delevación y escaneo N/S E/W.
    y = math.atan(s_z/s_x)
    x = math.asin(-s_y/math.sqrt(s_x**2 + s_y**2 + s_z**2))
    
    # De los ángulos de escaneo obtemos el pixel.
    
    # Si el array que contiene la variable X del .nc nos inidica que ángulo de escaneo le
    # ..  corresponde a cada pixel. ahora tenemos que encontrar "una situación inversa" , 
    # .. donde dado un ángulo de  escaneo en particular tenemos que encontrar su pixel. 
    # .. Esto no se puede hacer directo puesto que los ángulos de escaneo son números reales y la
    # .. posición de los pixeles se representa con enteros.
    # Para resolver este problema resto el ángulo de escaneo de nuestro interes con el array X, y
    # .. encuentro la posición o index del valor menor de esta diferencia.
    
    X_array = np.array(X)
    X_array = np.abs(X_array - x)
    px_x    = np.argmin(X_array)
    Y_array = np.array(Y)
    Y_array = np.abs(Y_array - y)
    px_y    = np.argmin(Y_array)
    return px_x , px_y


def cmap_banda13():
    """ 
    Obtiene un custom c_map adecuado para la banda 13, 
    (escala de temperaturas)
    """
    
    inicio = -110
    final  = 40
    dt     = final - inicio
    
    ini_gist_yarg = -110
    fin_gist_yarg = -78
    dy            = fin_gist_yarg - ini_gist_yarg
    
    ini_hsv = -78
    fin_hsv = -45
    dh = fin_hsv - ini_hsv
    
    ini_ocean = -45
    fin_ocean = -30
    do = fin_ocean - ini_ocean
    
    long_yarg = int(256*dy/dt)
    long_hsv  = int(256*dh/dt)
    long_do   = int(256*do/dt)
    long_db   = 256 - long_yarg - long_hsv - long_do
    
    gist_yarg = cm.get_cmap('gist_yarg', 256 )
    hsv       = cm.get_cmap('hsv'      , 256 )
    ocean     = cm.get_cmap("ocean"    , 256 )
    binary    = cm.get_cmap('binary'   , 256 )
    
    gist_yarg_parte  = gist_yarg(np.linspace(0,1,long_yarg  ))
    hsv_parte        =       hsv(np.linspace(0,0.29,long_hsv ))
    ocean_parte      =     ocean(np.linspace(0,1,long_do    ))
    binary_parte     =    binary(np.linspace(0.5,1,long_db    ))
    
    custom_cmap_array = np.concatenate([gist_yarg_parte,hsv_parte,ocean_parte,binary_parte])
    custom_cmap       = ListedColormap(custom_cmap_array)
    
    return custom_cmap
    
def obtener_ventana(topo,x,y,ventana=200):
    """
    Dado un par de pixeles (px_x , px_y) o coordenadas,
    obtiene un subarray cuadrado, de radio (ventana),
    a partir del array introducido (topo)
    """
    
    # Revisa que se respete los límites de la imágen.
    lim_izquierdo = max(x-ventana,0)
    lim_derecho   = min(x+ventana+1,topo.shape[1])
    lim_inferior  = max(y-ventana,0)
    lim_superior  = min(y+ventana+1,topo.shape[0])
    
    mensaje_aviso = "!! Aviso : Se ha alcanzado los límites de la imágen en el recorte, el resultado ya no será un array cuadrado."
    if lim_izquierdo == 0:
        lim_derecho = lim_izquierdo + ventana + 1
        print(mensaje_aviso)
        
    if lim_derecho == topo.shape[1]:
        lim_izquierdo = lim_derecho - ventana
        print(mensaje_aviso)
    if lim_inferior == 0:
        lim_superior = lim_inferior + ventana + 1
        print(mensaje_aviso)
    if lim_superior == topo.shape[0]:
        lim_inferior == lim_superior - ventana
        print(mensaje_aviso)
    if len(topo.shape) == 3:
        array = topo[ lim_inferior:lim_superior ,lim_izquierdo:lim_derecho,:]
    else:
        array = topo[ lim_inferior:lim_superior ,lim_izquierdo:lim_derecho]
    return array

def latlonArray(nc,data_path="",enviar_0_0 = True):
    """ 
        Cálcula las coordenadas de cada pixel de la variable principal
        del netCDF. El resultado son 2 arrays uno que contiene la latitud
        y otro que contiene la longitud, estos arrays estan en la forma
        de meshgrid.
        
        Función muy útil para la realización de proyecciones.
        
        Revisa si existen los arrays lat lon precalculados,
        de no existir los crea, en la carpeta especificada en 
        la variable como data_path.
        
        enviar_0_0 : Los pixeles fuera de la Tierra por default son
        ... marcados como np.inf , esto puede causar problemas. En particular
        ... cuando se usa para generar un pcolormesh con la libería Basemap de
        ... plt_toolkit.
        ... Esta opción si es marcada como True, mapea los np.inf a 0 y los envia a
        ... las coordenadas 0.0 , 0.0 (Fuera del continente americano)
        
        """
    
    # Obtenemos el tamaño del array.
    x_shape = np.array(nc.variables["x"]).shape[0]
    y_shape = np.array(nc.variables["y"]).shape[0]
    
    shape = (y_shape,x_shape)
    
    # Generamos nombre.
    if shape == (1500,2500):
        #nombre_h5 =  data_path + f"LAT_LON_G16_CONUS_{shape[0]}_{shape[1]}.h5"
        nombre_h5 = data_path +"CONUS/"+"Lat_Lon_CONUS_2km.h5"
    elif shape == (3000,5000):
        nombre_h5 = data_path +"CONUS/"+"Lat_Lon_CONUS_1km.h5"
    else:
        raise ValueError("Usar solo bandas de 1km y 2km de resolución, no la banda 2 con 0.5 km")
    
    # Creamos directorio si no existe.
    Path(data_path).mkdir(parents=True, exist_ok=True)
    
    # Buscamos en el path indicado.
    if os.path.exists(nombre_h5):
        with h5py.File(nombre_h5, "r") as file:
            lats = file["lats"][()]  # [()] <-- Para extraer los datos a np array.
            lons = file["lons"][()]
            file.close()
    # Si no existe,calculamos y guardamos.
    else:
        print("Aviso: No se encontró los arrays latlon, se calcularan unos nuevos.")
        # Obtenemos altura del satélite.
        sat_h = nc.variables['goes_imager_projection'].perspective_point_height
        # Obtenemos la longitud del satelite.
        sat_lon = nc.variables['goes_imager_projection'].longitude_of_projection_origin
        # Sweep del satélite ?? 
        sat_sweep = nc.variables['goes_imager_projection'].sweep_angle_axis
        
        # Obtenemos las coordenadas de la proyección geostacionaria.
        X = nc.variables['x']*sat_h
        Y = nc.variables['y']*sat_h
        # Calculamos los arrays.
        p = Proj(proj='geos', h=sat_h, lon_0=sat_lon, sweep=sat_sweep)
        XX, YY     = np.meshgrid(X, Y)
        lons, lats = p(XX, YY, inverse=True)

        # Mandamos a 0 los valores infinitos.
        # Dado que las coordendas 0.0,0.0 no estan en America, esta fuera de la vista,
        # ... y no afectara el cambio. Por otro lado si lo dejamos con valores infinitos
        # ... pcolormesh (para graficar los mapas) no funciona.
        if enviar_0_0:
            lons[lons == np.inf] = 0
            lats[lats == np.inf] = 0
        
        # Guardamos nuevos arrays antes de retornarlos.
        with h5py.File(nombre_h5, "w") as file:
            file.create_dataset("lats",data=lats,dtype=np.float32,compression="gzip")
            file.create_dataset("lons",data=lons,dtype=np.float32,compression="gzip")
        
        # Generamos los de méxico
        if shape == (1500,2500):
            lats_m = lats[625:,:1250]
            lons_m = lons[625:,:1250]
            nombre_h5 = data_path +"Mexico/"+"Lat_Lon_Mexico_2km.h5"
        elif shape == (3000,5000):
            lats_m = lats[1250:,:2500]
            lons_m = lons[1250:,:2500]
            nombre_h5 = data_path +"Mexico/"+"Lat_Lon_Mexico_1km.h5"
        with h5py.File(nombre_h5, "w") as file:
            file.create_dataset("lats",data=lats_m,dtype=np.float32,compression="gzip")
            file.create_dataset("lons",data=lons_m,dtype=np.float32,compression="gzip")

    return lats , lons


def _identificarBandas(df_files):
    """
    Le añade la información de a que banda pertence cada archivo, dado el nombre de un
    archivo netCDF diretamente descargado del los servidores usando regular expressions.
    """    
    Bandas = []
    for line in df_files["file"]:
        file_name = str(line)
        # Obtenemos los indices donde se encuentra la información de la banda. -M6C%%_
        # Nota, solo nos interesa las imágenes del Scan Mode 3 o 6, siendo el modo 6 "M6" el modelo por default del satélite.
        match  = re.search(r"-M6C\d\d_",file_name)
        if match == None:
            match = re.search(r"-M3C\d\d_",file_name)

        # debug
        if match is None:
            print("Ha salido None otravez...")
            print(file_name)
            print(df_files)
        
        span   = match.span()
        # Número de banda. (En string)
        banda  = file_name[span[1]-3:span[1]-1]
        Bandas.append(int(banda))
    df_files["Banda"] = Bandas
    return df_files


def _goes_file_df(satellite, product, start, end, refresh=True):
    """
    Get list of requested GOES files as pandas.DataFrame.
    Parameters
    ----------
    satellite : str
    product : str
    start : datetime
    end : datetime
    refresh : bool
        Refresh the s3fs.S3FileSystem object when files are listed.
        Default True will refresh and not use a cached list.
    """
    fs = s3fs.S3FileSystem(anon=True)
    
    DATES = pd.date_range(f"{start:%Y-%m-%d %H:00}", f"{end:%Y-%m-%d %H:00}", freq="1H")

    # List all files for each date
    # ----------------------------
    files = []
    for DATE in DATES:
        files += fs.ls(f"{satellite}/{product}/{DATE:%Y/%j/%H/}", refresh=refresh)

    # Build a table of the files
    # --------------------------
    df = pd.DataFrame(files, columns=["file"])
    df[["start", "end", "creation"]] = (
        df["file"].str.rsplit("_", expand=True, n=3).loc[:, 1:]
    )

    # Filter files by requested time range
    # ------------------------------------
    # Convert filename datetime string to datetime object
    df["start"] = pd.to_datetime(df.start, format="s%Y%j%H%M%S%f")
    df["end"] = pd.to_datetime(df.end, format="e%Y%j%H%M%S%f")
    df["creation"] = pd.to_datetime(df.creation, format="c%Y%j%H%M%S%f.nc")

    # Filter by files within the requested time range
    df = df.loc[df.start >= start].loc[df.end <= end].reset_index(drop=True)

    return df


def descargaIntervaloGOES16(producto,
                            datetime_inicio,
                            datetime_final,
                            banda=None,
                            output_path="NETCDF_DATA/",
                            verbose=False,
                            saltos=1):

    # Creamos el directorio si no existe.
    Path(output_path).mkdir(parents=True, exist_ok=True)
    
    # Nos conectamos a los servidores con credencial anónima. 
    fs = s3fs.S3FileSystem(anon=True)
    
    # Lista de productos
    lista_productos = fs.ls(f"noaa-goes16")

    # Asignamos fecha
    start = datetime_inicio
    end   = datetime_final

    # Obtenemos el dataframe con los elementos más recientes de cada banda.
    df = _goes_file_df(satellite="noaa-goes16",product=producto,start=start,end=end,refresh=True)

    
    # Identificamos cada archivo con la banda a la que corresponde.
    if banda != None:
        df = _identificarBandas(df) # Puse mas debug en la función.
        df = df[df["Banda"] == banda]
        
    # Aplicamos los saltos.
    df = df.iloc[::saltos]

    descargados = 0
    a_descargar = len(df)
    
    if a_descargar == 0:
        print("No se encontró ningun archivo por descargar.")
    else:
        # Descarga de los datos.
        for index in range(a_descargar):
            
            descarga_correcta = False

            file_name = df["file"].values[index]
            match  = re.search(r"OR_ABI.+",file_name)

            # debug
            if match is None:
                print("Ha salido None otravez...")
                print(producto)
                print(datetime_inicio)
                print(datetime_final)
                print(banda)
                print(file_name)
                print(df["file"])

            span   = match.span()
            output_name = file_name[span[0]:span[1]]

            # Si ya existe el archivo, continuamos.
            objeto_path = Path(output_path + output_name)
            if objeto_path.is_file():
                descargados += 1
                continue

            while descarga_correcta == False:
                try:
                    fs.get(file_name, output_path + output_name,)
                except KeyboardInterrupt:
                    raise
                except:
                    print("Error en la descarga, volviendo a intentar.")
                    time.sleep(5)
                else:
                    descarga_correcta = True
                    descargados += 1
            if verbose:
                print(f"Archivo descargado : \n{output_name}")
                print(f"Descargados {descargados} de {a_descargar}","\n")
        if verbose:
            print("Descargar completa.")



def datosActualesGOES16(producto,
                        banda=None,
                        output_name="GOES-descarga.nc"):
    """
    
    Descarga los datos más recientes de las categorias ingresadas, desde datos alojados en AWS.
    Los guarda en formato netCDF bajo el mismo nombre por los que se sobrescriben los datos.
    
    Cuando el producto es de clase ABI-L1b-RadC es necesario introducir la bada
    que se desea descargar.

    Basado en proyecto goes2go : https://github.com/blaylockbk/goes2go/
    
    LA FUNCIÓN SIGUE EN DESAROLLO, SOLO USAR CON PRODUCTOS EN EL DOMINIO DE CONUS.
    """
    
    # Nos conectamos a los servidores con credencial anónima. 
    fs = s3fs.S3FileSystem(anon=True)
    
    # Lista de productos
    lista_productos = fs.ls(f"noaa-goes16")
    
    # Revisamos si el producto solicitado está en la lista de productos disponibles.
    #assert(producto in lista_productos,f"El nombre del producto debe de ser {lista_productos}")
    
    # Obtenemos el intervalo de tiempo en el que buscaremos los archivos. (Hora UTC)
    start = datetime.datetime.utcnow() - datetime.timedelta(hours=1)
    end   = datetime.datetime.utcnow()
    

    # Obtenemos el dataframe con los elementos más recientes de cada banda.
    df = _goes_file_df(satellite="noaa-goes16",product=producto,start=start,end=end)
    df = df.loc[df.start == df.start.max()].reset_index(drop=True)

    # Identificamos cada archivo con la banda a la que corresponde.
    if banda != None:
        df = _identificarBandas(df)
        df = df[df["Banda"] == banda]
    
    #assert(len(df) > 0,"No se encontraron archivos.")
    if len(df) > 1 : 
        print("Aviso: Se encontró más de un archivo, solo se descargará el primero.")
        print(df,"\n")
        file_name = df["file"].values[0]
    else:
        # Obtenemos el nombre del archivo a descargar.
        file_name = df["file"].values[0]
    
    
    # Descarga de los datos.
    fs.get(file_name,output_name)
    print("Descargar completa.")


def estado_general(nc,variable):
    """
    Si el archivo netCDF tiene algún tipo de corrupción,
    se retorna False (archivo inválido), si no se detécta
    corrupción se marca como verdadero.

    Se revisan que los datos contengan al menos un valor que no sea
    fill_value.
    """
    fill_value = nc.variables[variable]._FillValue
    datos = np.array(nc.variables[variable])    
    if np.min(datos) == fill_value:
        archivo_valido = False
    else:
        archivo_valido = True
    return archivo_valido



def proyecciónMéxico():
    proyeccion = "+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.00335281068119356027 +lat_0=0.0 +lon_0=-89.5 +sweep=x +no_defs"
    p = pyproj.Proj(proyeccion)
    return p

def proyecciónMéxico_Coordenadas(lats,lons):
    p = proyecciónMéxico()
    lons , lats  = p(lons, lats)
    return lats,lons

# ------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
El código siguiente fue adaptado del proyecto Goes2Go.
"""

def normalize(value, lower_limit, upper_limit):
    norm = (value - lower_limit) / (upper_limit - lower_limit)
    norm = np.clip(norm, 0, 1)
    return norm

def gamma_correction(a, gamma):
    return np.power(a, 1 / gamma)

def breakpoint_stretch(C,breakpoint):
    """
    Contrast stretching by break point (number provided by Rick Kohrs)
    """
    lower = normalize(C, 0, 10)  # Low end
    upper = normalize(C, 10, 255)  # High end
    # Combine the two datasets
    # This works because if upper=1 and lower==.7, then
    # that means the upper value was out of range and the
    # value for the lower pass was used instead.
    combined = np.minimum(lower, upper)
    return combined

def GEOCOLOR(R,G,B,gamma=0.7):
    """
    Recibe como entrada los arrays de las bandas 1,2,3
    y genera una imágen en RGB.
    """

    # Reducimos la resolución de la banda roja.
    R = R[::2,::2]

    BASE = 400 # Base de la normalización
    
    R = normalize(R*1,0,BASE)
    B = normalize(B*0.7,0,BASE)
    G = normalize(G*1.2,0,BASE)
    
    print(np.max(R),np.max(G),np.max(G))

    R = np.clip(R, 0, 1)
    G = np.clip(G, 0, 1)
    B = np.clip(B, 0, 1)

    # Pseudo Green
    G = 0.45 * R + 0.1 * G + 0.45 * B
    G = np.clip(G, 0, 1)

    # Convert Albedo to Brightness, ranging from 0-255 K
    # (numbers based on email from Rick Kohrs)
    R = np.sqrt(R * 100) * 25.5
    G = np.sqrt(G * 100) * 25.5
    B = np.sqrt(B * 100) * 25.5

    R = breakpoint_stretch(R, 33)
    G = breakpoint_stretch(G, 40)
    B = breakpoint_stretch(B, 50) #50

    #B = B*0.6

    RGB = np.dstack([R, G, B])
    RGB = gamma_correction(RGB, gamma)
    
    return  RGB


class Producto:
    def __init__(self,path,QF_validos = None):
        # Abrimos el archivo nc y extraemos los datos
        self.path = path
        nc = netCDF4.Dataset(self.path)
        
        # Propiedades invariables
        self.producto = nc.title
        if self.producto == "ABI L2 Derived Stability Indices":
            self.variable = "CAPE"
        else:
            self.variable = list(nc.variables.keys())[0]
        self.banda    = self.obtener_banda(nc)
        self.t        = self._obtener_t(nc)
        self.datetime = self.obtener_datetime() # En UTC
        self.fill_value   = nc.variables[self.variable]._FillValue
        try:
            self.scale_factor = nc.variables[self.variable].scale_factor
        except AttributeError:
            self.scale_factor = None
        # Lista de QF_validos
        self.QF_validos = QF_validos 
        # Propiedades variables
        self.array = np.array(nc.variables[self.variable])
        try: 
            self.QF = np.array(nc.variables["DQF"])
        except KeyError:
            self.QF = np.array(nc.variables["DQF_Overall"])
        self.shape = self.array.shape
        self.resolución = self.obtener_resolución()

        # Primer chequeo general de los datos
        self.datos_validos = self.estado_general()

        nc.close

    def _obtener_t(self,nc):
        """
        Obtiene el tiempo en formato UNIX de los diferentes
        productos del GOES16.
        """
        try:
            t = float(np.array(nc.variables["t"]))
        except:
            t = float(np.array(nc.variables["time"])[0])
        return t
    
    def obtener_datetime(self):
        fecha_inicio = datetime.datetime(2000,1,1,12,0,0)
        time_delta   = datetime.timedelta(seconds=self.t)
        fecha        = fecha_inicio + time_delta
        return fecha

    def obtener_resolución(self):
        res_x_1km = 3000
        res_x = self.shape[0]
        return res_x_1km / res_x

    def zoom(self,factor):
        self.array = scipy.ndimage.zoom(self.array,factor,order=0)
        self.QF    = scipy.ndimage.zoom(self.QF,factor,order=0)
        self.shape = self.array.shape
        self.resolución = self.obtener_resolución()

    def downscale(self,factor):
        self.array = self.array[::factor,::factor]
        self.QF = self.QF[::factor,::factor]
        self.shape = self.shape,
        self.resolución = self.obtener_resolución()

    def obtener_banda(self,nc):
        """
        Obtiene la banda si está disponible. De lo contrario se
        returna None.
        """
        try:
            id_banda = np.array(nc.variables["band_id"])
            id_banda = int(id_banda)
        except KeyError:
            id_banda = None
        self.banda = id_banda
        return self.banda

    def estado_general(self):
        """
        Si el archivo netCDF tiene algún tipo de corrupción,
        se retorna False (archivo inválido), si no se detécta
        corrupción se marca como verdadero.
        Se revisan que los datos contengan al menos un valor que no sea
        fill_value.
        """

        test_fill_value = False
        if np.min(self.array) != self.fill_value:
            test_fill_value = True
        
        test_fecha = False
        if self.datetime.year > 2015:
            test_fecha = True

        resultado_test = test_fill_value and test_fecha
        return resultado_test

    def obtener_ventana(self,latitud,longitud,ventana):
        return Ventana(self,latitud,longitud,ventana)

class Ventana:
    """
    Toma como entrada un objeto de la clase Producto.
    Además de heredar la información básica de los datos
    que contiene le agrega la información de su ventana.
    """

    def __init__(self,Producto,latitud,longitud,ventana):

        self.producto = Producto.producto
        self.variable = Producto.variable
        self.banda    = Producto.banda
        self.datetime = Producto.datetime
        self.fill_value  = Producto.fill_value
        self.scale_value = Producto.scale_factor
        self.resolución  = Producto.resolución

        self.latitud  = latitud
        self.longitud = longitud
        self.ventana  = ventana

        # Extraemos la ventana
        px,py       = coordenadas2px(nc=self.resolución,latitud=latitud,longitud=longitud)
        ventana_px  = reescalado[self.resolución]*self.resolución
        array_datos = obtener_ventana(self.array,px,py,ventana_px)
        array_QF    = obtener_ventana(self.QF,px,py,ventana_px)

        self.array = array_datos
        self.QF   = array_QF

    def test_QF(self):
        """
        True si todos los pixeles del array pasan el test de QF
        False si no.
        """
        return None

