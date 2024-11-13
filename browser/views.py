import os
import pandas as pd
import plotly.express as px
from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse, FileResponse


def search_view(request):
    return render(request, 'browser/search.html')


def results_view(request):
    fbpp_id = request.GET.get('q')
    if not fbpp_id:
        return HttpResponse("Please provide a FBpp ID.")

    path = f'jobs/{fbpp_id}/{fbpp_id}_normPred_evolCombi.txt'
    full_path = os.path.join(settings.BASE_DIR, 'browser', 'static', path)
    if not os.path.exists(full_path):
        return HttpResponse("File not found.")

    # Leer el archivo en un DataFrame
    df = pd.read_csv(full_path, skiprows=1, sep=' ', header=None, index_col=0)
    reversed_oranges = px.colors.sequential.Oranges[::-1]
    
    # Generar el heatmap usando Plotly
    fig = px.imshow(df, color_continuous_scale=reversed_oranges, title=f"GEMME mutational landscape for {fbpp_id}")
    fig.update_layout(title_x=0.5)
    fig.update_layout(
        autosize=False,
        width=1500,
        height=500,
    )
    fig.update_yaxes(tickfont=dict(size=14), dtick=1)  # Aumentar tamaño de los ticks en el eje Y
    heatmap_html = fig.to_html(full_html=False)

    # Añadir las URLs de las imágenes al contexto
    image_path_1 = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}/{fbpp_id}_GMM.jpg')
    image_path_2 = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}/GMM.png')

    image_url_1 = None
    image_url_2 = None
    if os.path.exists(image_path_1):
        image_url_1 = f'/static/jobs/{fbpp_id}/{fbpp_id}_GMM.jpg'
    if os.path.exists(image_path_2):
        image_url_2 = f'/static/jobs/{fbpp_id}/GMM.png'

    return render(request, 'browser/results.html', {
        'heatmap_html': heatmap_html,
        'query': fbpp_id,
        'file_path': full_path,
        'image_url_1': image_url_1,
        'image_url_2': image_url_2
    })



def download_file(request, fbpp_id):
    file_path = f"jobs/{fbpp_id}/{fbpp_id}_normPred_evolCombi.txt"
    if not os.path.exists(file_path):
        return HttpResponse("File not found.")
    
    response = FileResponse(open(file_path, 'rb'))
    response['Content-Disposition'] = f'attachment; filename="{fbpp_id}_normPred_evolCombi.txt"'
    return response
