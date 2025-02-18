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
from django.core.mail import send_mail
from django.shortcuts import render

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

def past_jobs(request):
    return render(request, 'browser/past_jobs.html')

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
        job_status = 'computing conservation levels...'
        for i in os.listdir(f'/data/jobs/{job_id}/'):
            if '.out' in i:
               if os.path.exists('/data/jobs/'+job_id+'/'+str(i)):
                  checks = open('/data/jobs/'+job_id+'/'+str(i))
                  lines = checks.readlines()[-1]
                  checks.close()
                  if 'calculating' in lines or 'Running' in lines or 'Preparing' in lines or 'Generating' in lines:
                    job_status = lines
               break
        return JsonResponse({'status': job_status, 'message': 'At the end of the process the results will be displayed.'})

DATA = '/data/Drosophila_ProteoCast/'

@csrf_exempt
def upload_file(request):
#    if request.method == 'POST':
#        if 'file' in request.FILES:
#            uploaded_file = request.FILES['file']
#            return handle_upload(request, uploaded_file)
    if request.method == 'POST':
        files = request.FILES
        if 'file' in files and 'pdbFile' in files:
            uploaded_file = files['file']
            pdb_file = files['pdbFile']
            return handle_upload(request, uploaded_file, pdb_file)
        elif 'file' in files:
            main_file = files['file']
            return handle_upload(request, main_file, None)
    return JsonResponse({'error': 'Invalid request'}, status=400)

def handle_upload(request, uploaded_file, pdb_file):
    uniprot_id = request.POST.get('uniprotId')
    #if not uniprot_id:
    #    return JsonResponse({'status': 'error', 'message': 'No UniProt ID provided.'}, status=400)
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
        if pdb_file:
            pdb_file_path = os.path.join(folder_path, pdb_file.name)
            with open(pdb_file_path, 'wb+') as destination:
               for chunk in pdb_file.chunks():
                   destination.write(chunk)
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
#SBATCH --mail-user=marina.abakarova@sorbonne-universite.fr
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
                         
docker run --rm -v "/data/jobs/{job_id}:/opt/job" marinaabakarova/proteocast /bin/bash -c "cd / && bash run.sh {uploaded_file.name} {uniprot_id}"
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
    
def segmentation_dico(path_segF,path_bfactors):

    if not os.path.exists(path_segF):
        return None #HttpResponse("Segmentation file does not exist", status=404)
    if not os.path.exists(path_bfactors):
        return None #HttpResponse("B-factors file does not exist", status=404)
    
    dico_colors = {2:{'r': 182, 'g': 132, 'b': 187 },1:{'r': 243, 'g': 119, 'b': 140 }}
    df_segmentation = pd.read_csv(path_segF)
    df_bfactors = pd.read_csv(path_bfactors)
    if 'pLDDT' not in df_bfactors.columns:
        return None

    resi_70 = np.array(df_bfactors[df_bfactors['pLDDT'] <= 0.7].index)+1
    seg_dico=[]
    for _, row in df_segmentation.loc[df_segmentation['type']=='GEMME'].iterrows():

        start = int(row['start'])
        end = int(row['end'])
        state = row['state']
        
        if state == 1 or state == 2:
            segment = list(range(start, end + 1))
            sub_segments = []
            sub_segment = []
            
            for res in segment:
                if res in resi_70:
                    sub_segment.append(res)
                else:
                    if sub_segment:
                        sub_segments.append(sub_segment)
                        sub_segment = []
            
            if sub_segment:
                sub_segments.append(sub_segment)
            
            for sub_segment in sub_segments:
                dico = {
                    'start_residue_number': sub_segment[0],
                    'end_residue_number': sub_segment[-1],
                    'color': dico_colors[state],
                    'representation': 'putty',
                    'representationColor': dico_colors[state]
                }

                seg_dico.append(dico)
    return seg_dico


def results_view(request):
    prot_name = request.GET.get('q').lower()
 
    if not prot_name:
        return HttpResponse(f'Please provide a protein name.')
    
    pdb_id = ''
    msa_file_job = ''
    ## job
    if (prot_name[:3] == 'job'):
        data_path = '/data/jobs/'
        alias_dir = 'jobs'
        prot_name = prot_name[3:]
        id_folder = prot_name
        files = os.listdir(f'/data/jobs/{id_folder}')
        # Loop through filenames to find the first one with 'FBpp'
        prot_id = None
        for file_name in files:
            if "ProteoCast" in file_name:
                prot_id = '_'.join(file_name.split('.')[1].split('_')[:-1])  # Extract protein ID 
            if ('10.' in file_name): #and ('ResClass' not in file_name or 'Sensitivity' not in file_name): #10. in filname
                pdb_id = file_name.split('.')[1]
            if ('a3m' in file_name):
                msa_file_job = file_name
            elif ('2.ali' in file_name):
                msa_file_job = file_name[5:] 
    ## drosophila db
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
        ## getting pdb_id for the fly
        try:
            pdb_id = mapping_df.loc[prot_id, 'Structure_3D_file'].split('.')[0]
        except:
            pdb_id = ''

    ## reading SNPs file
    try:
        snps_file = f'{data_path}{id_folder}/7.{prot_id}_SNPs.csv'
        df_snps = pd.read_csv(snps_file)
        if df_snps.shape[0] == 0:
            df_snps = None
    except:
        df_snps = None

    # Basic checks
    if not data_path or not id_folder or not prot_id:
        message = 'ProteoCast file not found: Protein ID.'
        return render(request, 'browser/error.html', {'prot_id': prot_id,'message': message}, status=500)
        #return HttpResponse("Missing required path or protein ID.", status=500)

    alph = ["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"][::-1]
    alph = [i.upper() for i in alph]

    proteocast_path = f'{data_path}{id_folder}/4.{prot_id}_ProteoCast.csv'
    if not os.path.exists(proteocast_path):
        message = 'ProteoCast file not found'
        return render(request, 'browser/error.html', {'prot_id': prot_id,'message': message}, status=500)
#        return HttpResponse(f"ProteoCast file not found: {proteocast_path}", status=404)

    try:
        df_proteocast = pd.read_csv(proteocast_path)
    except Exception as e:
        message = 'ProteoCast could not read ProteoCast.csv file'
        return render(request, 'browser/error.html', {'prot_id': prot_id,'message': message}, status=500)
#        return HttpResponse(f"Error reading ProteoCast CSV: {e}", status=500)

    try:
        df_proteocast['LocalConfidence'] = df_proteocast['LocalConfidence'].replace({True: 1, False: 0})
        confidence_values = np.array(df_proteocast.groupby('Residue')['LocalConfidence']
                                     .apply(lambda x: x.iloc[0]).tolist()).reshape(1, -1)
    except Exception as e:
        confidence_values = None
        message = 'ProteoCast encountered a problem processing the confidence metric'
        return render(request, 'browser/error.html', {'prot_id': prot_id,'message': message}, status=500)
        
    try:
        df = pd.DataFrame(np.array(df_proteocast['Variant_score']).reshape(20, -1, order='F'))
        df_mut = pd.DataFrame(np.array(df_proteocast['Mutation']).reshape(20, -1, order='F'))
        
    except Exception as e:
        df = None
        message = 'ProteoCast has not found Variant_score or Mutation values in the ProteoCast.csv file'
        return render(request, 'browser/error.html', {'prot_id': prot_id,'message': message}, status=500)
    
    try:
        df_classes = pd.DataFrame(np.array(df_proteocast['Variant_class'].replace({'neutral': 1, 'uncertain': 2, 'impactful': 3})).reshape(20, -1, order='F'))
        df_classesStr = pd.DataFrame(np.array(df_proteocast['Variant_class']).reshape(20, -1, order='F'))
    except Exception as e:
        df_classes = None
        message = 'Variant classes are not found in the ProteoCast.csv file'
        return render(request, 'browser/error.html', {'prot_id': prot_id,'message': message}, status=500)

    variantClasses_colorscale = [

        [0, '#3688ED'],
        [0.5, '#E097CE'],
        [1, '#F25064']
    ]

    confidence_colorscale = [
        [0, 'white'],
        [1, 'darkblue']
    ]
    
    heatmap_html = ""
    heatmapClasses_html = ""
    heatmapSNPs_html = ""
    image_url_1 = ""
    fig_msarep = ""
    fig_segmentation = ""
    pdb_url_1 = ""
    pdb_url_2 = ""
    pdb_url_3 = ""

    ## GENERATING HEATMAPS
        #--- GEMME heatmap
    if df is not None:
        fig = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            row_heights=[0.9, 0.1],
            vertical_spacing=0.02,
        )
        heatmap_main = go.Heatmap(
            z=df.values[::-1],
            x=list(range(1, df.shape[1])),
            y=alph,
            colorscale=px.colors.sequential.Oranges[::-1],
            showscale=False,
            customdata=df_mut.values[::-1],
            hovertemplate=("Mutation: %{customdata}<br>"
                           "Score: %{z:.2f}<extra></extra>")
        )
        fig.add_trace(heatmap_main, row=1, col=1)
        
        #--- Variant classes heatmap
    if df_classes is not None:
        fig_VariantClasses = make_subplots(
                rows=2, cols=1,
                shared_xaxes=True,
                row_heights=[0.9, 0.1],
                vertical_spacing=0.02,
            )
        heatmap_classes = go.Heatmap(
            z=df_classes.values[::-1],
            x=list(range(1, df_classes.shape[1])),
            y=alph,
            customdata=np.dstack([df_mut.values[::-1], df_classesStr.values[::-1]]),
            colorscale=variantClasses_colorscale,
            showscale=False,
            hovertemplate=("Mutation: %{customdata[0]}<br>"
                           "Class: %{customdata[1]}<extra></extra>"),
            xgap=0.3,
            ygap=0.3,
        )
        fig_VariantClasses.add_trace(heatmap_classes, row=1, col=1)

        # --- SNPs heatmap
    if df_snps is not None:
        df_snps = df_snps.loc[df_snps['Set_name']!='Hypomorphic'].copy()
        fig_SNPs = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            row_heights=[0.9, 0.1],
            vertical_spacing=0.02,
        )

        # Prepare the SNPs data structure
        df_snps_STR = pd.DataFrame(columns=df_classesStr.columns, index=df_classesStr.index)
        for snp in df_snps['Mutation'].unique():
            ind_mut = alph.index(snp[-1])
            position = int(snp[1:-1])
            df_snps_STR.loc[ind_mut, position - 1] = '/'.join(df_snps.loc[df_snps['Mutation'] == snp, 'Set_name'].tolist())

        df_snps_STR = df_snps_STR.fillna('-')

        # Create a numerical mask for highlights
        highlight_mask = np.zeros(df_snps_STR.shape)  # Default is no highlight (0)
        highlight_mask[df_snps_STR.isin(['Lethal'])] = 1  # Red for Lethal
        highlight_mask[df_snps_STR.isin(['DEST2', 'DGRP', 'DEST2/DGRP', 'DGRP/DEST2'])] = 2  # Blue for DEST or DGRP

        #highlighted_positions = (highlight_mask > 0)  # Cells with highlights (red or blue)
        #df_modified = df.values[::-1].copy()  # Create a copy to avoid modifying the original dataframe
        #df_modified[highlighted_positions] = 0
        
        # Main heatmap (greyscale background)
        heatmap_snps = go.Heatmap(
            z=df.values[::-1],
            x=list(range(1, df.shape[1] + 1)),
            y=alph,
            customdata=np.dstack([df_mut.values[::-1], df_classesStr.values[::-1], df_snps_STR.values]),
            colorscale=px.colors.sequential.Greys[::-1][4:],
            showscale=False,
            hovertemplate=("Mutation: %{customdata[0]}<br>"
                        "Score: %{z:.2f}<br>"
                        "Class: %{customdata[1]}<br>"
                        "SNPs: %{customdata[2]}<extra></extra>")
        )

        # Highlight heatmap (overlay with colors for specific SNPs)
        highlight_layer = go.Heatmap(
            z=highlight_mask,  # Use the mask to determine colors
            x=list(range(1, df.shape[1] + 1)),
            y=alph,
            colorscale=[
                [0, "rgba(0,0,0,0)"],  # Transparent for no highlight
                [ 1 / 2, "rgba(255,50,50,1)"],  # Red for 'Lethal'[2 / 3, "rgba(0,0,255,0.8)"],  # Blue for 'DEST2' or 'DGRP'
                [1, "rgba(0,0,255,0.7)"],  # Blue continued
            ],
            showscale=False,
            hoverinfo="skip",  # Skip hover info for the highlight layer
            zmin=0, zmax=2
        )

        # Add traces to the figure
        fig_SNPs.add_trace(heatmap_snps, row=1, col=1)
        fig_SNPs.add_trace(highlight_layer, row=1, col=1)


    if confidence_values is not None:
        heatmap_confidence = go.Heatmap(
            z=confidence_values,
            x=list(range(1, df.shape[1])),
            colorscale=confidence_colorscale,
            showscale=False,
            hovertemplate="%{z}<extra></extra>",
            xgap=0.2, 
            zmin=0, zmax=1
        )
        scatter_border = go.Scatter(
            x=[0.5, df.shape[1]+0.5, df.shape[1]+0.5, 0.5, 0.5],
            y=[-0.5, -0.5, 0.5, 0.5, -0.5],
            mode="lines",
            line=dict(color="darkblue", width=2),
            hoverinfo="skip",
            showlegend=False,
        )
        if df is not None:
            fig.add_trace(heatmap_confidence, row=2, col=1)
            fig.add_trace(scatter_border, row=2, col=1)
            fig.update_layout(
                title_x=1,
                autosize=False,
                width=1500,
                height=600,
                xaxis=dict(
                    tickmode="array",
                    tickvals=list(range(1, df.shape[1], 10)),
                    ticktext=[str(i) for i in range(1, df.shape[1], 10)]
                ),
                yaxis=dict(title="Substituting amino acid"),
                xaxis2=dict(title="Residue")
            )
            fig.update_yaxes(visible=False, row=2, col=1)
            heatmap_html = fig.to_html(full_html=False)


        if df_classes is not None:
            fig_VariantClasses.add_trace(heatmap_confidence, row=2, col=1)
            fig_VariantClasses.add_trace(scatter_border, row=2, col=1)
            fig_VariantClasses.update_layout(
                title_x=1,
                autosize=False,
                width=1500,
                height=600,
                xaxis=dict(
                    tickmode="array",
                    tickvals=list(range(1, df.shape[1]+1, 10)),
                    ticktext=[str(i) for i in range(1, df.shape[1], 10)]
                ),
                yaxis=dict(title="Substituting amino acid"),
                xaxis2=dict(title="Residue")
            )
            fig_VariantClasses.update_yaxes(visible=False, row=2, col=1)
            heatmapClasses_html = fig_VariantClasses.to_html(full_html=False)

        if df_snps is not None:
            fig_SNPs.add_trace(heatmap_confidence, row=2, col=1)
            fig_SNPs.add_trace(scatter_border, row=2, col=1)
            fig_SNPs.update_layout(
                title_x=1,
                autosize=False,
                width=1500,
                height=600,
                xaxis=dict(
                    tickmode="array",
                    tickvals=list(range(1, df.shape[1]+1, 10)),
                    ticktext=[str(i) for i in range(1, df.shape[1], 10)]
                ),
                yaxis=dict(title="Substituting amino acid"),
                xaxis2=dict(title="Residue")
            )
            fig_SNPs.update_yaxes(visible=False, row=2, col=1)
            heatmapSNPs_html = fig_SNPs.to_html(full_html=False)


    image_url_1 = f'/{alias_dir}/{id_folder}/6.{prot_id}_GMM.jpg'
    fig_msarep = f'/{alias_dir}/{id_folder}/3.{prot_id}_msaRepresentation.jpg'
    #fig_segmentation = f'/{alias_dir}/{id_folder}/9.{prot_id}_SegProfile.png'

    for file_path in [image_url_1, fig_msarep]:#, fig_segmentation]:
        check_path = file_path.replace(alias_dir, data_path)
        if not os.path.exists(check_path):
            return HttpResponse(f"File not found: {check_path}", status=404)

    pdb_url_1 = f'/{alias_dir}/{id_folder}/10.{pdb_id}.pdb'
    pdb_url_2 = f'/{alias_dir}/{id_folder}/12.{pdb_id}_Sensitivity.pdb'
    pdb_url_3 = f'/{alias_dir}/{id_folder}/11.{pdb_id}_ResClass.pdb'
    pdb_check = None
    if pdb_url_1:
        pdb_check = pdb_url_1.replace(alias_dir, data_path)
        if not os.path.exists(pdb_check):
            pdb_url_1 = None
    if pdb_url_2:
        pdb_check = pdb_url_2.replace(alias_dir, data_path)
        if not os.path.exists(pdb_check):
            pdb_url_2 = None

    if pdb_url_3:
        pdb_check = pdb_url_3.replace(alias_dir, data_path)
        if not os.path.exists(pdb_check):
            pdb_url_3 = None




    ###### SEGMENTATION ######

    if os.path.exists(f'{data_path}{id_folder}/13.{prot_id}_GEMME_pLDDT.csv'):
        df_segmentation = pd.read_csv(f'{data_path}{id_folder}/8.{prot_id}_Segmentation.csv')
        df_sefPrep = pd.read_csv(f'{data_path}{id_folder}/13.{prot_id}_GEMME_pLDDT.csv')
        
        
        # Initialize plot
        fig_Seg = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            row_heights=[0.1, 0.9],
            vertical_spacing=0.05,
        )
        
        if 'pLDDT' in df_sefPrep.columns:
            
            # Prepare the one-row heatmap data
            bfactors = df_sefPrep['pLDDT'].tolist()
            # Define palette and bin ranges
            bins = [0, 0.5, 0.7, 0.9, 1]
            df_plddt = pd.DataFrame(bfactors, columns=['value'])

            # Bin the values into groups
            df_plddt['group'] = pd.cut(df_plddt['value'], bins=bins, labels=[0, 1, 2, 3], right=False)
            grouped_values_numeric = np.array(df_plddt['group'])

            # Reshape to a one-row matrix
            data_matrix = grouped_values_numeric.reshape(1, -1)

            # Map colors
            heatmap_colorscale = [
                [0, '#FF7D45'],  # Bin 0
                [0.33, '#FFDB13'],  # Bin 1
                [0.66, '#65CBF3'],  # Bin 2
                [1, '#0053D6']  # Bin 3
            ]

            # Add the one-row heatmap
            heatmapSeg = go.Heatmap(
                z=data_matrix,
                colorscale=heatmap_colorscale,
                showscale=False,
                customdata=(df_sefPrep['pLDDT'].values * 100).astype(int).reshape(1, -1),
                x=list(range(1, len(bfactors) + 1)),
                y=[''],
                hovertemplate="Residue %{x}<br>pLDDT: %{customdata}<extra></extra>",
                xgap=0.15, zmin=0, zmax=3
            )
            fig_Seg.add_trace(heatmapSeg, row=1, col=1)

        
        
        # Prepare the trace plot data
        GEMME_mean = df_sefPrep['GEMME_mean'].tolist()
        n_res = len(GEMME_mean)
        signal = pd.DataFrame({
            'x': list(range(1, n_res + 1)),
            'y': list(GEMME_mean),
            'type': ['GEMME'] * n_res
        })
        # Segmentation data
        gemme = df_segmentation[df_segmentation['type'] == 'GEMME']
        cp = pd.DataFrame({
            'x': list(gemme['end'][:-1]),
            'type': ['GEMME'] * (len(gemme) - 1)
        })
        
        mean = pd.DataFrame({
            'x': list(gemme['start'] - 1),
            'xend': list(gemme['end']),
            'y': list(gemme['mean']),
            'type': ['GEMME'] * len(gemme),
            'state': list(gemme['state'])
        })

        # Add the GEMME trace plot
        trace = go.Scatter(
            x=signal['x'],
            y=signal['y'],
            mode='lines',
            line=dict(color='#2777B4'),  # Changed color to green in hexadecimal
            showlegend=False,
            hovertemplate="Residue %{x}<br>GEMME Score: %{y:.2f}<extra></extra>"
        )
        fig_Seg.add_trace(trace, row=2, col=1)

        # Add segmentation regions (axvspan equivalent)
        for _, row in mean.iterrows():
            fig_Seg.add_shape(
                type="rect",
                x0=row['x']+0.5, x1=row['xend']+0.5,
                y0=0, y1=max(GEMME_mean),
                line=dict(width=0),
                fillcolor='white' if row['state'] == 0 else 'red' if row['state'] == 1 else 'purple',
                opacity=0.1,
                row=2, col=1
            )
            
            fig_Seg.add_shape(
                type="line",
                x0=row['x']+0.5, x1=row['xend']+0.5,
                y0=row['y'], y1=row['y'],
                line=dict(color="darkblue", width=2),
                row=2, col=1
            )

        # Add segmentation change points (axvline equivalent)
        for _, row in cp.iterrows():
            fig_Seg.add_shape(
                type="line",
                x0=row['x']+0.5, x1=row['x']+0.5,
                y0=0, y1=max(GEMME_mean),
                line=dict(color="red", width=2),
                row=2, col=1
            )
            # Add horizontal lines for each segment
        
        # Define the frame for the second plot
        scatter_frame = go.Scatter(
            x=[0.5, len(GEMME_mean) + 0.5, len(GEMME_mean) + 0.5, 0.5, 0.5],
            y=[0, 0, 1, 1, 0],  # The y-range matches the axis limits (0 to 1)
            mode="lines",
            line=dict(color="black", width=2),  # Adjust the color and width of the frame
            hoverinfo="skip",  # No hover effect for the frame
            showlegend=False  # No legend for the frame
        )

        # Add the frame to the second subplot
        fig_Seg.add_trace(scatter_frame, row=2, col=1)
        # Update layout
        fig_Seg.update_layout(
            width=1500,
            height=600,
            plot_bgcolor="white",  # Set the background color to white
            paper_bgcolor="white",
            xaxis2=dict(title="Residue Index"),
            yaxis2=dict(title="Average GEMME Score", range=[0, 1]),
            hovermode="x unified"
        )
        # Generate HTML for Django
        fig_segmentation = fig_Seg.to_html(full_html=False)


    
    warning_message = ''
    if set(confidence_values[0]) == {0}:
        warning_message = 'Unreliable Predictions Warning: Evolutionary information in the provided MSA is scarse (<200 sequences) and insufficient to provide reliable predictions.'

    ## segmentation data for 3D
    seg_dico = segmentation_dico(f'{data_path}{id_folder}/8.{prot_id}_Segmentation.csv', f'{data_path}{id_folder}/13.{prot_id}_GEMME_pLDDT.csv') 


    return render(request, 'browser/results.html', {
        'heatmap_html': heatmap_html,
        'heatmapClasses_html': heatmapClasses_html,
        'heatmapSNPs_html': heatmapSNPs_html,
        'query': id_folder,
        'prot_name': prot_id,
        'image_url_1': image_url_1,
        'pdb_url_1': pdb_url_1,
        'pdb_url_2': pdb_url_2,
        'pdb_url_3': pdb_url_3,
        'fig_msarep': fig_msarep,
        'fig_segmentation': fig_segmentation,
        'select_segments': seg_dico,
        'warning_message':warning_message,
        'pdb_file': pdb_id+'.pdb',
        'msa_file_job':msa_file_job
    })

def download_folder(request, fbpp_id):
    # Path to the folder to be downloaded
    if fbpp_id[:4] =='2025':
        folder_path = os.path.join('/data/jobs', fbpp_id)
    else:
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
