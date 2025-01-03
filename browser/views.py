from django.urls import reverse
# from django.shortcuts import redirect
import os
import pandas as pd
import plotly.express as px
from django.shortcuts import render
# from django.conf import settings
from django.http import HttpResponse, FileResponse
import numpy as np
import shutil
# from io import BytesIO
# import zipfile
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from django.views.decorators.csrf import csrf_exempt
# from django.core.files.storage import FileSystemStorage
from datetime import datetime
import subprocess
# from django.http import QueryDict
# import uuid
from django.http import JsonResponse


def contact_us(request):
    return render(request, 'browser/contact_us.html')

def documentation(request):
    return render(request, 'browser/documentation.html')

def citation(request):
    return render(request, 'browser/publications.html')

def search_view(request):
    return render(request, 'browser/search.html')

def drosophiladb(request):
    return render(request, 'browser/drosophiladb.html')
def download(request):
    return render(request, 'browser/download.html')

def job_running(request,job_id): 
    job_status_path = os.path.join('/data/jobs/', job_id, 'status.txt')

    if not os.path.exists(job_status_path):
        return render(request, 'browser/error.html', {'message': 'Status file not found for this job.'}, status=404)
    try:
        with open(job_status_path, 'r') as status_file:
            status = status_file.read().strip()
    except Exception as e:
        return render(request, 'browser/error.html', {'message': f'Error reading status file: {str(e)}'}, status=500)

    return render(request, 'browser/job_running.html', {'job_id': job_id, 'status': status})

def check_job_status(request):
    job_id = request.GET.get('job_id')
    if not job_id:
        return JsonResponse({'status': 'error', 'message': 'No job_id provided.'}, status=400)

    job_status_path = os.path.join('/data/jobs', job_id, 'status.txt')

    if not os.path.exists(job_status_path):
        return JsonResponse({'status': 'error', 'message': f'Status file not found at {job_status_path}'}, status=404)

    with open(job_status_path, 'r') as file:
        job_status = file.read().strip()

    if job_status == 'finished':
        # Construct the URL with the job ID as a query parameter
        job_id = 'job' + job_id
        results_url = f"{reverse('results')}?q={job_id}"
        return JsonResponse({'status': 'finished', 'redirect_url': results_url})
    else:
        return JsonResponse({'status': job_status, 'message': 'At the end of the process the results will be displayed.'})

DATA = '/data/Drosophila_ProteoCast/'

@csrf_exempt
def upload_file(request):
    if request.method == 'POST':
        if 'file' in request.FILES:
            uploaded_file = request.FILES['file']
            return handle_upload(request, uploaded_file)
    
    return JsonResponse({'error': 'Invalid request'}, status=400)

def handle_upload(request, uploaded_file):
    now = datetime.now()
    global job_id
    job_id = now.strftime('%Y%m%d%H%M%S')
    folder_path = os.path.join('/data/jobs/', job_id)
    os.makedirs(folder_path, mode=0o755, exist_ok=True)

    file_path = os.path.join(folder_path, uploaded_file.name)

    try:
        with open(file_path, 'wb+') as destination:
            for chunk in uploaded_file.chunks():
                destination.write(chunk)
        with open(file_path, 'r') as file:
            first_line = file.readline().strip()
        prot_name = first_line.lstrip('>')
        #new_folder_path = '/data/jobs/' + prot_name
        #job_id = prot_name
        #if not os.path.exists(new_folder_path):
        #    os.rename(folder_path, new_folder_path)
        #    os.chdir(new_folder_path)
        run_docker_script = os.path.join(folder_path, 'run_docker.sh')
        with open(run_docker_script, 'w') as script:
            script.write(f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --job-name=test
#SBATCH --mail-type=END
#SBATCH --mail-user=abakamarina@gmail.com
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
                         
docker run --rm -v "/data/jobs/{job_id}:/opt/job" elodielaine/gemme:gemme /bin/bash -c "cd / && bash run.sh {uploaded_file.name}"
""")
        os.chmod(run_docker_script, 0o755)
        os.chdir(folder_path)
        subprocess.run(['sbatch', run_docker_script], check=True)

        
        job_status_path = os.path.join(folder_path, 'status.txt')
        with open(job_status_path, 'w') as status_file:
            status_file.write('in_progress')
        return JsonResponse({
                'status': 'in_progress',
                'job_id': job_id,
                'redirect_url': f'/results/job_running/{job_id}'})
    except Exception as e:
        return JsonResponse({'status': 'error', 'message': str(e)}, status=500)


def serve_file(request, folder, filename):
    file_path = os.path.join(DATA, folder, filename)
    if os.path.exists(file_path):
        response = FileResponse(open(file_path, 'rb'))
        return response
    else:
        return HttpResponse(f"File not found: {filename}", status=404)


def results_view(request):
    prot_name = request.GET.get('q').lower()

    if not prot_name:
        return HttpResponse(f'Please provide a protein name.')
    
    if prot_name[:3] == 'job':
        data_path = '/data/jobs/'
        alias_dir = 'job'
        prot_name = prot_name[3:]
        id_folder = prot_name
        path = os.path.join(data_path, id_folder)
        proteocast_path = f'{data_path}{id_folder}/4.FBpp0070001_ProteoCast.csv'
        if not os.path.exists(proteocast_path):
            return HttpResponse("ProteoCast file not found.")
        df = pd.read_csv(proteocast_path)
        path = os.listdir('data/jobs/')
        if prot_name:
           return HttpResponse(f'Read succ.{path}') 
        '''# Loop through filenames to find the first one with 'FBpp'
        prot_id = None
        for file_name in files:
            if "ProteoCast" in file_name:
                prot_id = file_name.split('.')[1].split('_')[0]  # Extract protein ID before the first dot
                break'''
    else:
        # Only for the fly
        data_path = '/data/Drosophila_ProteoCast/'
        alias_dir = 'data'
        mapping_file_path = f'{data_path}mapping_database.csv'
        if not os.path.exists(mapping_file_path):
            return HttpResponse("Mapping file not found.")

        mapping_df = pd.read_csv(mapping_file_path, index_col=0)
        
        mapping_df['fbpp_low'] = mapping_df.index.str.lower()
        mapping_df['pr_sym'] = mapping_df['Protein_symbol'].str.lower()
        if prot_name[:4] == 'fbpp':
            id_folder = mapping_df.loc[mapping_df['fbpp_low'] == prot_name, 'id'].item()
            prot_id = mapping_df.loc[mapping_df['id'] == id_folder].index[0]
        else:
            id_folder = mapping_df.loc[mapping_df['pr_sym'] == prot_name, 'id'].item()
            prot_id = mapping_df.loc[mapping_df['id'] == id_folder].index[0]
    
    alph = ["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"][::-1]
    alph = [i.upper() for i in alph]
    
    # Confidence values
    proteocast_path = f'{data_path}{id_folder}/4.{prot_id}_ProteoCast.csv'
    if not os.path.exists(proteocast_path):
        return HttpResponse("ProteoCast file not found.")
    
    df_proteocast = pd.read_csv(proteocast_path)
    df_proteocast['GEMME_LocalConfidence'] = df_proteocast['GEMME_LocalConfidence'].replace({True: 1, False: 0})
    confidence_values = np.array(df_proteocast.groupby('Residue')['GEMME_LocalConfidence'].apply(lambda x: x.iloc[0]).tolist()).reshape(1, -1)
    
    df = pd.DataFrame(np.array(df_proteocast['GEMME_score']).reshape(20, -1, order='F'))
    df_classes = pd.DataFrame(np.array(df_proteocast['Variant_class'].replace({'neutral': 1, 'uncertain': 2, 'impactful': 3})).reshape(20, -1, order='F'))
    df_classesStr = pd.DataFrame(np.array(df_proteocast['Variant_class']).reshape(20, -1, order='F'))
    df_mut = pd.DataFrame(np.array(df_proteocast['Mutation']).reshape(20, -1, order='F'))

    # Custom color scale for confidence
    variantClasses_colorscale = [
        [0, '#3688ED'],  # neutral- white
        [0.5, '#E097CE'],  # uncertain - purple
        [1, '#F25064']      # impactful - red
    ]
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

    fig_VariantClasses = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        row_heights=[0.9, 0.1],
        vertical_spacing=0.02,
    )
    heatmap_main = go.Heatmap(
        z=df.values[::-1],  # Reverse the order of rows
        x=list(range(1, df.shape[1])),  
        y=alph,  # Reverse the order of y-axis labels
        colorscale=px.colors.sequential.Oranges[::-1],
        showscale=False, 
        customdata=df_mut.values[::-1],  # Reverse the order of custom data
        hovertemplate=(
            "Position: %{customdata}<br>"
            "Score: %{z:.2f}<extra></extra>"
        ),
    )
    heatmap_classes = go.Heatmap(
        z=df_classes.values[::-1],
        x=list(range(1, df_classes.shape[1])),
        y=alph,
        customdata=np.dstack([df_mut.values[::-1], df_classesStr.values[::-1]]),
        colorscale=variantClasses_colorscale,
        showscale=False,
        hovertemplate=(
            "Position: %{customdata[0]}<br>"
            "Class: %{customdata[1]}<extra></extra>"
        ),
        # These just control the gap size; their color is determined by the figure background.
        xgap=0.3,
        ygap=0.3,
    )

    fig.add_trace(heatmap_main, row=1, col=1)
    fig_VariantClasses.add_trace(heatmap_classes, row=1, col=1)
    
    heatmap_confidence = go.Heatmap(
       z=confidence_values,
       x=list(range(1, df_classes.shape[1])),  
       colorscale=confidence_colorscale,
       showscale=False, 
       hovertemplate=(
            "%{z}<extra></extra>"
        ),
    )  
    
    fig.add_trace(heatmap_confidence, row=2, col=1)
    fig_VariantClasses.add_trace(heatmap_confidence, row=2, col=1)
    
    scatter_border = go.Scatter(
       x=[0, df.shape[1]+1, df.shape[1]+1, 0, 0],
       y=[-0.5, -0.5, 0.5, 0.5, -0.5],
       mode="lines",
       line=dict(color="darkblue", width=2), 
       hoverinfo="skip",
       showlegend=False,
    )
    fig.add_trace(scatter_border, row=2, col=1)
    fig_VariantClasses.add_trace(scatter_border, row=2, col=1)
    fig.update_layout(title_x=1, autosize=False, width=1500, height=600, xaxis=dict(tickmode="array", tickvals=list(range(1, df.shape[1], 10)), ticktext=[str(i) for i in range(1, df.shape[1], 10)]), yaxis=dict(title="Substituting amino acid"), xaxis2=dict(title="Residue"))
    fig_VariantClasses.update_layout(title_x=1, autosize=False, width=1500, height=600, xaxis=dict(tickmode="array", tickvals=list(range(1, df.shape[1]+1, 10)), ticktext=[str(i) for i in range(1, df.shape[1], 10)]), yaxis=dict(title="Substituting amino acid"),  xaxis2=dict(title="Residue"))

    fig.update_yaxes(visible=False, row=2, col=1)
    fig_VariantClasses.update_yaxes(visible=False, row=2, col=1)
    heatmap_html = fig.to_html(full_html=False)
    heatmapClasses_html = fig_VariantClasses.to_html(full_html=False)

    image_url_1 = f'/{alias_dir}/{id_folder}/6.{prot_id}_GMM.jpg'
    fig_msarep = f'/{alias_dir}/{id_folder}/3.{prot_id}_msaRepresentation.jpg'
    fig_segmentation = f'/{alias_dir}/{id_folder}/9.{prot_id}_SegProfile.png'

    # Validate existence of files
    for file_path in [image_url_1, fig_msarep, fig_segmentation]:
        if not os.path.exists(file_path.replace('/data/', data_path)):
            return HttpResponse(f"File not found: {file_path}")
        
    pdb_url_1 = f'/{alias_dir}/{id_folder}/AF-Q45VV3-F1-model_v4.pdb'
    if not os.path.exists(pdb_url_1.replace('/data/', data_path)):
        pdb_url_1 = None
    
    return render(request, 'browser/results.html', {
        'heatmap_html': heatmap_html,
        'heatmapClasses_html': heatmapClasses_html,
        'query': id_folder,
        'prot_name': prot_id,
        'image_url_1': image_url_1,
        'pdb_url_1': pdb_url_1,
        'fig_msarep': fig_msarep,
        'fig_segmentation': fig_segmentation,
    })

def download_folder(request, fbpp_id):
    # Path to the folder to be downloaded
    folder_path = os.path.join('/data/Drosophila_ProteoCast', fbpp_id)
    
    if not os.path.exists(folder_path):
        return HttpResponse("Folder not found.", status=404)
    
    # Create a temporary ZIP file in the system's temporary directory
    temp_zip_path = os.path.join('/tmp', f'{fbpp_id}.zip')
    
    try:
        # Create a ZIP archive of the folder
        shutil.make_archive(temp_zip_path.replace('.zip', ''), 'zip', folder_path)
        
        # Serve the ZIP file
        response = FileResponse(open(temp_zip_path, 'rb'), as_attachment=True)
        response['Content-Disposition'] = f'attachment; filename="{fbpp_id}.zip"'
        
        return response
    except Exception as e:
        return HttpResponse(f"An error occurred: {e}", status=500)
    finally:
        # Ensure the temporary ZIP file is deleted after serving
        if os.path.exists(temp_zip_path):
            os.remove(temp_zip_path)
