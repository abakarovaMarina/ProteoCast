from django.urls import reverse
from django.shortcuts import redirect
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
from django.views.decorators.csrf import csrf_exempt
from django.core.files.storage import FileSystemStorage
from datetime import datetime
import subprocess 
import uuid
from django.http import JsonResponse, HttpResponseRedirect

def contact_us(request):
    return render(request, 'browser/contact_us.html')

def documentation(request):
    return render(request, 'browser/documentation.html')

def citation(request):
    return render(request, 'browser/citation.html')

def search_view(request):
    return render(request, 'browser/search.html')

def drosophiladb(request):
    return render(request, 'browser/drosophiladb.html')

def job_running(request, job_id):
    return render(request, 'browser/job_running.html', {'job_id': job_id})


#DATA = '/Users/manchuta/Documents/GitHub/Droso_GEMMEwebsite/browser/static/jobs/Drosophila_ProteoCast/' #'/data/Drosophila_ProteoCast/'
DATA = '/data/Drosophila_ProteoCast/'

@csrf_exempt
def upload_file(request):
    if request.method == 'POST' and 'file' in request.FILES:
        uploaded_file = request.FILES['file']
        now = datetime.now()
        job_id = now.strftime('%Y-%m-%d_%H-%M-%S')
        folder_name = 'jobs/' + job_id
        folder_path = os.path.join(settings.BASE_DIR, 'browser', 'static', folder_name)

        os.makedirs(folder_path, mode=0o755, exist_ok=True)

        file_path = os.path.join(folder_path, uploaded_file.name)

        try:
            
            with open(file_path, 'wb+') as destination:
                for chunk in uploaded_file.chunks():
                    destination.write(chunk)

            with open(file_path, 'r') as file:
                first_line = file.readline().strip()

            new_folder_name = f"jobs/{first_line.lstrip('>')}"
            new_folder_path = os.path.join(settings.BASE_DIR, 'browser', 'static', new_folder_name)
            os.rename(folder_path, new_folder_path)

            docker_image = "elodielaine/gemme:gemme" 
            container_workdir = "/opt/job"

            subprocess.run(
                [
                    "docker", "run", "--rm",
                    "-v", f"{new_folder_path}:{container_workdir}", 
                    docker_image,
                    "bash", "-c", f"cd / && bash run.sh {first_line.lstrip('>')}"
                ],
                check=True  
            )

            return HttpResponseRedirect(f'/results/?q={uploaded_file.name[:-4]}')
        except subprocess.CalledProcessError as e:
            return JsonResponse({'error': f"Error running Docker: {e}"}, status=500)
        except Exception as e:
            return JsonResponse({'error': str(e)}, status=500)

    return JsonResponse({'error': 'Invalid request'}, status=400)

def results_view_forjobs(request):
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
    reversed_oranges = px.colors.sequential.Oranges[::-1]
    df.fillna(0.0, inplace=True)
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
        z=df.values[::-1],  # Reverse the order of rows
        x=list(range(1, df.shape[1] + 1)),
        y=alph,  # Reverse the order of y-axis labels
        colorscale=px.colors.sequential.Oranges[::-1],
        customdata=df_mut.values[::-1],  # Reverse the order of custom data
        hovertemplate=(
            "Position: %{customdata}<br>"
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
       x=[0, df.shape[1], df.shape[1], 0, 0],
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
            tickvals=list(range(0, df.shape[1]+1, 10)),
            ticktext=[str(i) for i in range(0,df.shape[1]+1, 10)],
        ),
        yaxis=dict(title="Mutations"),
        coloraxis_colorbar=dict(
            title="GEMME Score",
            tickvals=[-8, -4, 0],
        ),
    )

    fig.update_yaxes(visible=False, row=2, col=1)
    heatmap_html = fig.to_html(full_html=False)

    image_url_1 = f'{DATA}{id_folder}/6.{FBpp_id}_GMM.jpg'
    pdb_url_1 = f'{DATA}{id_folder}/AF-Q45VV3-F1-model_v4.pdb'
    fig_msarep = f'{DATA}{id_folder}/3.{FBpp_id}_msaRepresentation.jpg'
    fig_segmentation = f'{DATA}{id_folder}/9.{FBpp_id}_SegProfile.png'

    return render(request, 'browser/results.html', {
        'heatmap_html': heatmap_html,
        'query': FBpp_id,
        'image_url_1': image_url_1,
        'pdb_url_1': pdb_url_1,
        'fig_msarep': fig_msarep,
        'fig_segmentation': fig_segmentation,
    })

def results_view(request):
    prot_name = request.GET.get('q').lower()
    if not prot_name:
        return HttpResponse(f'Please provide a protein name.')

    alph = ["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"][::-1]
    alph = [i.upper() for i in alph]
    mapping_file_path = f'{DATA}mapping_database.csv'
    if not os.path.exists(mapping_file_path):
        return HttpResponse("Mapping file not found.")
    
    mapping_df = pd.read_csv(mapping_file_path, index_col=0)
    
    mapping_df['fbpp_low'] = mapping_df.index.str.lower()
    mapping_df['pr_sym'] = mapping_df['Protein_symbol'].str.lower()
    

    if prot_name[:4]=='fbpp':
        
        id_folder = mapping_df.loc[mapping_df['fbpp_low']==prot_name, 'id'].item()
        FBpp_id =mapping_df.loc[mapping_df['id']==id_folder].index[0]
    else:
        id_folder = mapping_df.loc[mapping_df['pr_sym']==prot_name, 'id']
        FBpp_id =mapping_df.loc[mapping_df['fbpp_low']==prot_name].index[0]
   
    
    ### Confidence values
    proteocast_path = f'{DATA}{id_folder}/4.{FBpp_id}_ProteoCast.csv'
    if not os.path.exists(proteocast_path):
        return HttpResponse("ProteoCast file not found.")
    
    
    df_proteocast = pd.read_csv(proteocast_path)
    df_proteocast['GEMME_LocalConfidence'] = df_proteocast['GEMME_LocalConfidence'].replace({True:1, False:0})
    confidence_values = np.array(df_proteocast.groupby('Residue')['GEMME_LocalConfidence'].apply(lambda x: x.iloc[0]).tolist()).reshape(1, -1)
    
    
    df = pd.DataFrame(np.array(df_proteocast['GEMME_score']).reshape(20, -1, order='F'))
    df_mut = pd.DataFrame(np.array(df_proteocast['Mutation']).reshape(20, -1, order='F'))

    

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
        z=df.values[::-1],  # Reverse the order of rows
        x=list(range(1, df.shape[1] + 1)),  
        y=alph,  # Reverse the order of y-axis labels
        colorscale=px.colors.sequential.Oranges[::-1], 
        customdata=df_mut.values[::-1],  # Reverse the order of custom data
        hovertemplate=(
            "Position: %{customdata}<br>"
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
       x=[0, df.shape[1], df.shape[1], 0, 0],
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
            tickvals=list(range(0, df.shape[1]+1, 10)),  
            ticktext=[str(i) for i in range(0,df.shape[1]+1, 10)],
        ),
        yaxis=dict(title="Mutations"),
        coloraxis_colorbar=dict(
            title="GEMME Score",
            tickvals=[-8, -4, 0], 
        ),
    )

    fig.update_yaxes(visible=False, row=2, col=1)
    heatmap_html = fig.to_html(full_html=False)

    image_url_1 = f'{DATA}{id_folder}/6.{FBpp_id}_GMM.jpg'
    pdb_url_1 = f'{DATA}{id_folder}/AF-Q45VV3-F1-model_v4.pdb'
    fig_msarep = f'{DATA}{id_folder}/3.{FBpp_id}_msaRepresentation.jpg'
    fig_segmentation = f'{DATA}{id_folder}/9.{FBpp_id}_SegProfile.png'

    if not os.path.exists(image_url_1):
        return HttpResponse(f'GMM file not found {image_url_1}.')
    if not os.path.exists(pdb_url_1):
        return HttpResponse("PDB file not found.")
    if not os.path.exists(fig_msarep):
        return HttpResponse("MSA file not found.")
    if not os.path.exists(fig_segmentation):
        return HttpResponse("Seg file  not found.")
    print(image_url_1)
    return render(request, 'browser/results.html', {
        'heatmap_html': heatmap_html,
        'query': FBpp_id,
        'image_url_1': image_url_1,
        'pdb_url_1': pdb_url_1,
        'fig_msarep': fig_msarep,
        'fig_segmentation': fig_segmentation,
    })

def download_folder(request, fbpp_id):
    folder_path = f'{DATA}{fbpp_id}'
    if not os.path.exists(folder_path):
        return HttpResponse("Folder not found.")
    
    zip_path = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}.zip')
    shutil.make_archive(zip_path.replace('.zip', ''), 'zip', folder_path)
    
    response = FileResponse(open(zip_path, 'rb'))
    #response['Content-Disposition'] = f'attachment; filename="{fbpp_id}.zip"' 

    # Delete the zip file after sending the response
    response['delete_zip'] = zip_path
    if os.path.exists(zip_path):
        os.remove(zip_path) 
    return response

