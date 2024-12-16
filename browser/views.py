import os
import pandas as pd
import plotly.express as px
from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse, FileResponse
import numpy as np
import shutil
from io import BytesIO
import zipfile
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.core.files.storage import FileSystemStorage
from datetime import datetime

def contact_us(request):
    return render(request, 'browser/contact_us.html')

def documentation(request):
    return render(request, 'browser/documentation.html')

def citation(request):
    return render(request, 'browser/citation.html')

def search_view(request):
    return render(request, 'browser/search.html')

@csrf_exempt
def upload_file(request):
    if request.method == 'POST' and 'file' in request.FILES:
        uploaded_file = request.FILES['file']
        now = datetime.now()
        folder_name = 'jobs/'+now.strftime('%Y-%m-%d_%H-%M-%S')
        folder_path = os.path.join(settings.BASE_DIR, 'browser', 'static', folder_name)
        
        try:
            os.makedirs(folder_path, exist_ok=True)
            
            os.chmod(folder_path, 0o755)  
        except PermissionError:
            return JsonResponse({'error': 'Permission denied: unable to create directory.'}, status=500)
        except Exception as e:
            return JsonResponse({'error': f'Error creating directory: {str(e)}'}, status=500)
        
        file_storage = FileSystemStorage(location=folder_path)
        filename = file_storage.save(uploaded_file.name, uploaded_file)
        return JsonResponse({'message': 'File uploaded successfully', 'file_path': os.path.join(folder_name, filename)})
    
    return JsonResponse({'error': 'Invalid request'}, status=400)

def results_view(request):
    fbpp_id = request.GET.get('q')
    if not fbpp_id:
        return HttpResponse("Please provide a FBpp ID.")

    path = f'jobs/{fbpp_id}/Sensitivity_Confidence.csv'
    full_path = os.path.join(settings.BASE_DIR, 'browser', 'static', path)

    conf_df = pd.read_csv(full_path)
    confidence_values = conf_df['Confidence'].values.reshape(1, -1)

    seq = ''
    fasta = f'jobs/{fbpp_id}/{fbpp_id}.fasta'
    full_fasta = os.path.join(settings.BASE_DIR, 'browser', 'static', fasta)
    if os.path.exists(full_fasta):
        with open(full_fasta) as fst:
            for line in fst:
                line = line.rstrip('\n')
                if '>' not in line:
                    seq += line

    path = f'jobs/{fbpp_id}/{fbpp_id}_normPred_evolCombi.txt'
    full_path = os.path.join(settings.BASE_DIR, 'browser', 'static', path)
    if not os.path.exists(full_path):
        return HttpResponse("File not found.")

    df = pd.read_csv(full_path, skiprows=1, sep=' ', header=None, index_col=0)

    positions = list(range(1, len(seq) + 1))
    mutations = list(df.index)

    customdata = np.empty((df.shape[0], df.shape[1], 4), dtype=object) 
    for row_idx, mutation in enumerate(mutations):
        for col_idx in range(df.shape[1]):
            position = positions[col_idx]
            native_residue = seq[col_idx]
            gemme_score = df.iloc[row_idx, col_idx]
            customdata[row_idx, col_idx] = [position, native_residue, mutation.upper(), gemme_score]

    # Custom color scale for confidence
    confidence_colorscale = [
        [0, 'white'],  # Low confidence - white
        [1, 'darkblue']  # High confidence - dark blue
    ]

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        row_heights=[0.9, 0.1],
        vertical_spacing=0.02,
    )

    heatmap_main = go.Heatmap(
        z=df.values,  
        x=list(range(1, len(seq) + 1)),  
        y=list(df.index), 
        colorscale=px.colors.sequential.Oranges[::-1], 
        colorbar=dict(title="GEMME Score"), 
        customdata=customdata,  
        hovertemplate=(
            "Position: %{customdata[1]}%{x}%{customdata[2]}<br>"
            "Score: %{z:.2f}<extra></extra>"
        ),
    )
    fig.add_trace(heatmap_main, row=1, col=1)

    heatmap_confidence = go.Heatmap(
       z=confidence_values,
       colorscale=confidence_colorscale,
       showscale=False, 
       hoverinfo='skip',  
    )  
    fig.add_trace(heatmap_confidence, row=2, col=1)

    scatter_border = go.Scatter(
       x=[0, len(seq), len(seq), 0, 0],
       y=[-0.5, -0.5, 0.5, 0.5, -0.5],
       mode="lines",
       line=dict(color="black", width=2), 
       hoverinfo="skip",
       showlegend=False
    )
    fig.add_trace(scatter_border, row=2, col=1)

    fig.update_layout(
        title_x=1,
        autosize=False,
        width=1500,
        height=600,
        xaxis=dict(
            tickmode="array",
            tickvals=list(range(0, len(positions), 10)),  
            ticktext=[str(positions[i]) for i in range(0, len(positions), 10)],
        ),
        yaxis=dict(title="Mutations"),
        coloraxis_colorbar=dict(
            title="GEMME Score",
            tickvals=[-8, -4, 0], 
        ),
    )
    fig.update_yaxes(visible=False, row=2, col=1)
    heatmap_html = fig.to_html(full_html=False)

    image_path_1 = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}/{fbpp_id}_GMM.jpg')
    pdb_path_1 = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}/{fbpp_id}.pdb')
    fig_path_3  = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}/3.{fbpp_id}_msaRepresentation.jpg')
    fig_path_4  = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}/9.{fbpp_id}_SegProfile.png')
    image_url_1 = f'/static/jobs/{fbpp_id}/{fbpp_id}_GMM.jpg' if os.path.exists(image_path_1) else None
    pdb_url_1 = f'/static/jobs/{fbpp_id}/{fbpp_id}.pdb' if os.path.exists(pdb_path_1) else None
    fig_msarep = f'/static/jobs/{fbpp_id}/3.{fbpp_id}_msaRepresentation.jpg' if os.path.exists(fig_path_3) else None
    fig_segmentation = f'/static/jobs/{fbpp_id}/9.{fbpp_id}_SegProfile.png' if os.path.exists(fig_path_4) else None
    return render(request, 'browser/results.html', {
        'heatmap_html': heatmap_html,
        'query': fbpp_id,
        'file_path': full_path,
        'image_url_1': image_url_1,
        'pdb_url_1': pdb_url_1,
        'fig_msarep': fig_msarep,
        'fig_segmentation' : fig_segmentation,
    })

def download_folder(request, fbpp_id):
    folder_path = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}')
    if not os.path.exists(folder_path):
        return HttpResponse("Folder not found.")
    
    zip_path = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}.zip')
    shutil.make_archive(zip_path.replace('.zip', ''), 'zip', folder_path)
    
    response = FileResponse(open(zip_path, 'rb'))
    response['Content-Disposition'] = f'attachment; filename="{fbpp_id}.zip"' 

    # Delete the zip file after sending the response
    response['delete_zip'] = zip_path
    if os.path.exists(zip_path):
        os.remove(zip_path) 
    return response

