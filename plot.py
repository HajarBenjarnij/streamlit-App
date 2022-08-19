from re import L
import streamlit as st
from PIL import Image 
st.set_page_config(
   page_title="Prévision des phénomènes méteorologiques",
   page_icon="🌦️",
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
import emoji
import smtplib
import ssl
from email.message import EmailMessage
import datetime

with st.spinner('Wait for it...'):
    def plot(data_file):
        date = str(datetime.datetime.today()).split(' ', 1)[0].replace('-', '')

        shapefile_path = r"C:\Users\hajar\Desktop\prediction\SIG complement\sous_bassin.shp"
        gdf_shapefile = gpd.read_file(shapefile_path, epsg=4326)
        datafile= data_file.format(date)
        df=pd.read_csv(datafile,delimiter=',')
        ordre=True
        for i in range(len(df)):
            if df['Total Precipitation'][i]>=40.0:
                ordre=False
                break
            else:
                continue
        
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
        #masking the data outside of the limit
        ax.add_patch(patch)
        m.readshapefile('C:\\Users\\hajar\\Desktop\\prediction\\SIG complement\\Délimitation_administrative_des_ABH', 'Délimitation_administrative_des_ABH', linewidth= 1.5)
        today = datetime.date.today()
        tomorrow = today + datetime.timedelta(days = 1) 

        min = 0.0
        if z1.min() < 0:
            min == 0.0
        else:
            min == z1.min()

        title_text = f"Carte de vigilance de {tomorrow}  Etabli le {today} à 11h00"
        m.drawmapscale(-5., 22, -2, 25.7, 1000, barstyle='fancy')

        #plt.title(title_text, size=15)
        donnes=pd.read_csv(data_file.format(date))
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

        #the_table = plt.table(cellText=df.values,
                #colLabels=df.columns, colWidths=[0.24,0.15,0.15],
                #cellLoc = 'left', rowLoc = 'center',
                #transform=plt.gcf().transFigure, loc='bottom', 
                #)
        #plt.subplots_adjust(bottom=0.05)



        #the_table.auto_set_font_size(False)
        #the_table.set_fontsize(12)
        #logo = plt.imread('logo.png')

        fig.subplots_adjust(0.05, 0.05, 0.97)
        #ax.figure.figimage(logo, 20, 20, alpha=.9, zorder=10)
        fig1 = plt.gcf()
        #to display the plot
        plt.draw()
        #plt.show()
        #fig1.savefig('C:\Users\hajar\Desktop\prediction\{} Maroc.tif'.format(today),dpi=200,bbox_inches='tight') #to save the figure in tif format
        
        #df["Moyenne en mm"] = df["Moyenne en mm"].astype(int)
        plt.close('all')
        return fig1,ordre,df

    original_title = '<p style="font-family:Arial; color:Gray; font-size: 30px;">Tableau de Bord</p>'
    st.sidebar.markdown(original_title, unsafe_allow_html=True)
    #st.sidebar.header("Tableau de Bord")
    for i in range(1):
        st.sidebar.write('                 ')
    def user_input():
        st.sidebar.write('        \n         ')
        st.sidebar.subheader('\nNombre des heurs de prédiction')
        danger=st.sidebar.slider('Selectionner',max_value=72,min_value=24,value=24,step=24)
        st.sidebar.subheader('\nEnvoyer Un mail au cas dangereus ⚠️!')
        n=st.sidebar.selectbox('Selectionner',options=['On','Off'],index=0)
        return danger,n


    image = Image.open(r'C:\\Users\\hajar\\Desktop\\Prediction\\Logo.png')

    st.sidebar.image(image,width=250)
    #st.sidebar.image(r'C:\\Users\\hajar\\Desktop\\Prediction\\logometo.png')
    #increment = st.sidebar.button('Envoyer un e-mail au cas urgent')
    danger,hour=user_input()
    header=st.container()
    dataset=st.container()
    model_training=st.container()
    features=st.container()
    def sendEmail():
        # Define email sender and receiver
        email_sender = 'hajarbenjarnij@gmail.com'
        email_password = 'usakqnhkibkdjnte'
        email_receiver = 'hajarbenjarnij@gmail.com'

        today = datetime.date.today()
        tomorrow = today + datetime.timedelta(days = 1) 
        # Set the subject and body of the email
        subject = emoji.emojize('Attention :zipper-mouth_face: !')
        body = "Bonjour\n\n Demain {} il y'a des points ou la précipitation dépasse 40 mm !\n\nService de Gestion des Phénomènes Extrêmes\n\nDivision de Gestion de l'Eau et des Phénomènes Extrêmes\n\nDirection de la Recherche et de la Planification de l'Eau".format(tomorrow)


        em = EmailMessage()
        em['From'] = email_sender
        em['To'] = email_receiver
        em['Subject'] = subject
        em.set_content(body)

        # Add SSL (layer of security)
        context = ssl.create_default_context()

        # Log in and send the email
        with smtplib.SMTP_SSL('smtp.gmail.com', 465, context=context) as smtp:
            smtp.login(email_sender, email_password)
            smtp.sendmail(email_sender, email_receiver, em.as_string())
    with dataset:
        tab1, tab2 = st.tabs(["🗺️ Carte de vigilance ", "🌧️ Des informations sur les Bassins"])
         
        date = str(datetime.datetime.today()).split(' ', 1)[0].replace('-', '')
        today = datetime.date.today()
        tomorrow = today + datetime.timedelta(days = 1) 
        figure,Ordre,df=plot(r'C:\\Users\\hajar\\Desktop\\prediction\\CSVTOTAL{}.csv'.format(date))
        
        #st.subheader("Carte de vigilance de demain établi le {} à 11h00\n\n".format(today))
        tab1.write(figure)
        st.text('\n\n')
        st.text('\n\n')
        if Ordre==False:
            sendEmail()
        with tab2.expander("Plus d'informations sur les Bassins"):
            #st.write('\n\n')
            col1,col2,col3,col4,col5=st.columns(5)
            col1.metric(label=df['Nom du Bassin'][0], value="{} mm".format(df['Moyenne en mm'][0]), delta="{} mm".format(df['Maximum en mm'][0]))
            col2.metric(label=df['Nom du Bassin'][1], value="{} mm".format(df['Moyenne en mm'][1]), delta="{} mm".format(df['Maximum en mm'][1]))
            col3.metric(label=df['Nom du Bassin'][2], value="{} mm".format(df['Moyenne en mm'][2]), delta="{} mm".format(df['Maximum en mm'][2]))
            col4.metric(label=df['Nom du Bassin'][3], value="{} mm".format(df['Moyenne en mm'][3]), delta="{} mm".format(df['Maximum en mm'][3]))
            with col5:
                st.metric(label=df['Nom du Bassin'][4], value="{} mm".format(df['Moyenne en mm'][4]), delta="{} mm".format(df['Maximum en mm'][4]))
            st.text('\t')
            col6,col7,col8,col9,col10=st.columns(5)
            with col6:
                st.metric(label=df['Nom du Bassin'][5], value="{} mm".format(df['Moyenne en mm'][5]), delta="{} mm".format(df['Maximum en mm'][5]))
            with col7:
                st.metric(label=df['Nom du Bassin'][6], value="{} mm".format(df['Moyenne en mm'][6]), delta="{} mm".format(df['Maximum en mm'][6]))
            with col8:
                st.metric(label=df['Nom du Bassin'][7], value="{} mm".format(df['Moyenne en mm'][7]), delta="{} mm".format(df['Maximum en mm'][7]))
            with col9:
                st.metric(label=df['Nom du Bassin'][8], value="{} mm".format(df['Moyenne en mm'][8]), delta="{} mm".format(df['Maximum en mm'][8]))
            with col10:
                st.metric(label=df['Nom du Bassin'][9], value="{} mm".format(df['Moyenne en mm'][9]), delta="{} mm".format(df['Maximum en mm'][9]))
            st.text('\t')
            col11,col12,col13,col14,col15=st.columns(5)
            with col11:
                st.metric(label=df['Nom du Bassin'][10], value="{} mm".format(df['Moyenne en mm'][10]), delta="{} mm".format(df['Maximum en mm'][10]))
            with col12:
                st.metric(label=df['Nom du Bassin'][11], value="{} mm".format(df['Moyenne en mm'][11]), delta="{} mm".format(df['Maximum en mm'][11]))
            with col13:
                st.metric(label=df['Nom du Bassin'][12], value="{} mm".format(df['Moyenne en mm'][12]), delta="{} mm".format(df['Maximum en mm'][12]))
            with col14:
                st.metric(label=df['Nom du Bassin'][13], value="{} mm".format(df['Moyenne en mm'][13]), delta="{} mm".format(df['Maximum en mm'][13]))
            with col15:
                st.metric(label=df['Nom du Bassin'][14], value="{} mm".format(df['Moyenne en mm'][14]), delta="{} mm".format(df['Maximum en mm'][14]))
            st.text('\t')
            col16,col17,col18,col19,col20=st.columns(5)
            with col16:
                st.metric(label=df['Nom du Bassin'][15], value="{} mm".format(df['Moyenne en mm'][15]), delta="{} mm".format(df['Maximum en mm'][15]))
            with col17:
                st.metric(label=df['Nom du Bassin'][16], value="{} mm".format(df['Moyenne en mm'][16]), delta="{} mm".format(df['Maximum en mm'][16]))
            with col18:
                st.metric(label=df['Nom du Bassin'][17], value="{} mm".format(df['Moyenne en mm'][17]), delta="{} mm".format(df['Maximum en mm'][17]))
            with col19:
                st.metric(label=df['Nom du Bassin'][18], value="{} mm".format(df['Moyenne en mm'][18]), delta="{} mm".format(df['Maximum en mm'][18]))
            with col20:
                st.metric(label=df['Nom du Bassin'][19], value="{} mm".format(df['Moyenne en mm'][19]), delta="{} mm".format(df['Maximum en mm'][19]))
            st.text('\t')
            col21,col22,s,r=st.columns(4)
            with col22:
                st.metric(label=df['Nom du Bassin'][20], value="{} mm".format(df['Moyenne en mm'][20]), delta="{} mm".format(df['Maximum en mm'][20]))
            with s:
                st.metric(label=df['Nom du Bassin'][21], value="{} mm".format(df['Moyenne en mm'][21]), delta="{} mm".format(df['Maximum en mm'][21]))
            
st.snow()




