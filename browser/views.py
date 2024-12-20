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

def job_running(request):
    return render(request, 'browser/job_running.html')


#DATA = 'browser/static/jobs/Drosophila_ProteoCast/' #'/data/Drosophila_ProteoCast/'
DATA = '/data/Drosophila_ProteoCast/'

@csrf_exempt
def upload_file(request):
    if request.method == 'POST' and 'file' in request.FILES:

        uploaded_file = request.FILES['file']
        now = datetime.now()
        job_id = now.strftime('%Y-%m-%d_%H-%M-%S')
        #folder_name = 'jobs/' + job_id
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
            new_folder_path = '/data/jobs/' + prot_name
            os.rename(folder_path, new_folder_path)
            os.chdir(new_folder_path)
            
            run_docker_script = os.path.join(new_folder_path, 'run_docker.sh')
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


docker run --rm -v "/data/jobs/{prot_name}:/opt/job" elodielaine/gemme:gemme /bin/bash -c "cd / && bash run.sh {uploaded_file.name}"
                            """)
            os.chmod(run_docker_script, 0o755)
            subprocess.run(['sbatch', run_docker_script], check=True)
            job_status_path = os.path.join(new_folder_path, 'status.txt')
            with open(job_status_path, 'w') as status_file:
                status_file.write('in_progress') 
            return redirect('job_running')

        except subprocess.CalledProcessError as e:
            return JsonResponse({'error': f"Error running Docker: {e}"}, status=500)
        except Exception as e:
            return JsonResponse({'error': str(e)}, status=500)

    return JsonResponse({'error': 'oups Invalid request'}, status=400)


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
        id_folder = mapping_df.loc[mapping_df['pr_sym']==prot_name, 'id'].item()
        FBpp_id =mapping_df.loc[mapping_df['id']==id_folder].index[0]
   
    
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

    image_url_1 = f'/data/{id_folder}/6.{FBpp_id}_GMM.jpg'
    pdb_url_1 = f'/data/{id_folder}/AF-Q45VV3-F1-model_v4.pdb'
    fig_msarep = f'/data/{id_folder}/3.{FBpp_id}_msaRepresentation.jpg'
    fig_segmentation = f'/data/{id_folder}/9.{FBpp_id}_SegProfile.png'

    # Validate existence of files
    for file_path in [image_url_1, pdb_url_1, fig_msarep, fig_segmentation]:
        if not os.path.exists(file_path.replace('/data/', DATA)):
            return HttpResponse(f"File not found: {file_path}")
    
    return render(request, 'browser/results.html', {
        'heatmap_html': heatmap_html,
        'query': id_folder,
        'prot_name': FBpp_id,
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

