import streamlit as st
st.set_page_config(
   page_title="Prévision des phénomènes méteorologiques",
   page_icon=":wave:",
   layout="wide",
   initial_sidebar_state="expanded",
)
from typing import ClassVar
from PIL.Image import EXTENT
import numpy as np
from numpy.core.fromnumeric import clip
import pandas as pd
import glob
from pykrige.ok import OrdinaryKriging
from pykrige.kriging_tools import write_asc_grid
import pykrige.kriging_tools as kt
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Path, PathPatch
import geopandas as gpd
import datetime
import matplotlib as mpl
from shapely.geometry import Point
import matplotlib.cbook as cbook
import matplotlib.image as image

date = str(datetime.datetime.today()).split(' ', 1)[0].replace('-', '')

shapefile_path = r"C:\Users\hajar\Desktop\prediction\SIG complement\sous_bassin.shp"
gdf_shapefile = gpd.read_file(shapefile_path, epsg=4326)
datafile= r'C:\Users\hajar\Desktop\prediction\CSVTOTAL{}.csv'.format(date)
df=pd.read_csv(datafile,delimiter=',')
lons=np.array(df['Longitude'])
lats=np.array(df['Latitude'])
data=np.array(df['Total Precipitation'])
grid_space = 0.01
grid_lon = np.arange(np.amin(lons), np.amax(lons), grid_space) #grid_space is the desired delta/step of the output array
grid_lat = np.arange(np.amin(lats), np.amax(lats), grid_space)
OK = OrdinaryKriging(lons, lats, data, variogram_model='spherical', verbose=True, enable_plotting=False,nlags=20)
z1, ss1 = OK.execute('grid', grid_lon, grid_lat, backend="C", n_closest_points=12)
xintrp, yintrp = np.meshgrid(grid_lon, grid_lat)
fig, ax = plt.subplots(figsize=(10,10))
m = Basemap(llcrnrlon=lons.min()-0.1,llcrnrlat=lats.min()-0.1,urcrnrlon=lons.max()+0.1,urcrnrlat=lats.max()+0.1, projection='merc', resolution='h',area_thresh=1000.,ax=ax)
m.drawcoastlines() #draw coastlines on the map
x,y=m(xintrp, yintrp) # convert the coordinates into the map scales
ln,lt=m(lons,lats)
norm = mpl.colors.Normalize(vmin=0, vmax=50, clip=True)
cmap = LinearSegmentedColormap.from_list('my colormap',["#F8F9Fa","#96D2FA","#1E6EEB","#28BE28","#73F06E","#FFFAAA","#FFDE69","#FFA000","#E61900","#A50000"],N=100)




cs=m.pcolormesh(x, y, z1, cmap=cmap, norm=norm) #plot the data on the map.
cbar=m.colorbar(cs,location='right',pad="7%", label="Précipitations\nen mm") #plot the colorbar on the map
# draw parallels.
parallels = np.arange(21,40.0,2)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14, linewidth=0.0) #Draw the latitude labels on the map

# draw meridians
meridians = np.arange(-19,2,3)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14, linewidth=0.0)

##getting the limits of the map:
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
map_edges = np.array([[x0,y0],[x1,y0],[x1,y1],[x0,y1]])
##getting all polygons used to draw the coastlines of the map
polys = [p.boundary for p in m.landpolygons]

##combining with map edges
polys = [map_edges]+polys[:]
##creating a PathPatch
codes = [
[Path.MOVETO]+[Path.LINETO for p in p[1:]]
for p in polys
]

polys_lin = [v for p in polys for v in p]

codes_lin = [xx for cs in codes for xx in cs]

path = Path(polys_lin, codes_lin)
patch = PathPatch(path,facecolor='white', lw=0)
##masking the data outside of the limit
#ax.add_patch(patch)
#m.readshapefile('C:\Users\hajar\Desktop\prediction\SIG complement\Délimitation_administrative_des_ABH', 'Délimitation_administrative_des_ABH', linewidth= 1.5)
today = datetime.date.today()
tomorrow = today + datetime.timedelta(days = 1) 

min = 0.0
if z1.min() < 0:
    min == 0.0
else:
    min == z1.min()

title_text = f"Carte de vigilance de {tomorrow}  Etabli le {today} à 11h00"
m.drawmapscale(-5., 22, -2, 25.7, 1000, barstyle='fancy')

plt.title(title_text, size=15)
donnes=pd.read_csv(r'C:\Users\hajar\Desktop\prediction\CSVTOTAL{}.csv'.format(date))
geodonnees=gpd.GeoDataFrame(donnes,geometry=gpd.points_from_xy(donnes['Longitude'], donnes['Latitude']))
Projection = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'
geodonnees.to_file(filename = r'C:\Users\hajar\Desktop\prediction\CSVTOTAL.shp',driver='ESRI Shapefile',crs = Projection)
df_precip = gpd.read_file(r'C:\Users\hajar\Desktop\prediction\CSVTOTAL.shp')
df_precip["geometry"] = [Point(xy) for xy in zip(df_precip["Longitude"], df_precip["Latitude"])]
region_name = []
region_max =[]
region_mean =[]

for i, region in enumerate(gdf_shapefile["geometry"]):
 
    gdf_shapefile_region = gdf_shapefile[i:i+1]
 
    df_precip_region = df_precip.copy()
    region_poly = gdf_shapefile_region["geometry"].tolist()[0]
    df_precip_region.loc[:, "region_mask"] = df_precip_region["geometry"].apply(lambda x: x.within(region_poly))
    df_precip_region = df_precip_region[df_precip_region["region_mask"]]
 
    gdf_precip_region = gpd.GeoDataFrame(df_precip_region, crs=4326)
    
    region_name.insert(i,gdf_shapefile_region["NOM"].tolist()[0])
    region_max.insert(i, gdf_precip_region["Total Prec"].max())
    region_mean.insert(i, gdf_precip_region["Total Prec"].mean())
my_formatted_moy = [ '%.0f' % elem for elem in region_mean ]
my_formatted_max = [ '%.0f' % elem for elem in region_max ]
data = {'Nom du Bassin': region_name,
        'Maximum en mm': my_formatted_max,
        'Moyenne en mm': my_formatted_moy}
df = pd.DataFrame(data)

the_table = plt.table(cellText=df.values,
          colLabels=df.columns, colWidths=[0.24,0.15,0.15],
          cellLoc = 'left', rowLoc = 'center',
          transform=plt.gcf().transFigure, loc='bottom', 
          )
plt.subplots_adjust(bottom=0.05)



the_table.auto_set_font_size(False)
the_table.set_fontsize(12)
#logo = plt.imread('logo.png')

fig.subplots_adjust(0.05, 0.05, 0.97)
#ax.figure.figimage(logo, 20, 20, alpha=.9, zorder=10)
fig1 = plt.gcf()
 #to display the plot
plt.draw()
#plt.show()
#fig1.savefig('C:\Users\hajar\Desktop\prediction\{} Maroc.tif'.format(today),dpi=200,bbox_inches='tight') #to save the figure in tif format
plt.close('all')
original_title = '<p style="font-family:Arial; color:Gray; font-size: 30px;">Tableau de Bord</p>'
st.sidebar.markdown(original_title, unsafe_allow_html=True)
#st.sidebar.header("Tableau de Bord")
for i in range(5):
  st.sidebar.write('        \n         ')

def user_input():
    st.sidebar.write('        \n         ')
    max_depth=st.sidebar.slider('Choisir un nombre',1,3,1)
    st.sidebar.write('        \n         ')
    n=st.sidebar.selectbox('Number houres prediction',options=['24 heures','48 heures','72 heures'],index=0)
    return max_depth,n
st.sidebar.write('        \n         ')
increment = st.sidebar.button('Envoyer un e-mail au cas urgent')
S,hour=user_input()
header=st.container()
dataset=st.container()
model_training=st.container()
features=st.container()
st.title('Prévision des phénomènes méteorologiques')


    
with dataset:
    st.header('Hello')
    st.write(fig)

