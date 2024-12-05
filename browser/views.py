import os
import pandas as pd
import plotly.express as px
from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse, FileResponse
import numpy as np



def contact_us(request):
    return render(request, 'browser/contact_us.html')

def documentation(request):
    return render(request, 'browser/documentation.html')

def citation(request):
    return render(request, 'browser/citation.html')

def search_view(request):
    return render(request, 'browser/search.html')

def results_view(request):
    fbpp_id = request.GET.get('q')
    if not fbpp_id:
        return HttpResponse("Please provide a FBpp ID.")

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

    fig = px.imshow(
        df,
        color_continuous_scale=reversed_oranges,
        title=f"GEMME mutational landscape for {fbpp_id}",
        labels={'x': 'Position', 'y': 'Mutations', 'color': 'GEMME Score'},
    )

    fig.update_traces(
        customdata=customdata,
        hovertemplate=(
            "%{customdata[1]}%{customdata[0]}%{customdata[2]}<br>"
            "%{customdata[3]:.2f}<extra></extra>"
        )
    )

    fig.update_layout(
        title_x=0.5,
        autosize=False,
        width=1500,
        height=500,
        xaxis=dict(
            tickmode='array',
            tickvals=list(range(0, len(positions), 10)),
            ticktext=[str(positions[i]) for i in range(0, len(positions), 10)],
        ),
    )

    heatmap_html = fig.to_html(full_html=False)

    image_path_1 = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}/{fbpp_id}_GMM.jpg')
    pdb_path_1 = os.path.join(settings.BASE_DIR, 'browser', 'static', f'jobs/{fbpp_id}/{fbpp_id}.pdb')
    image_url_1 = f'/static/jobs/{fbpp_id}/{fbpp_id}_GMM.jpg' if os.path.exists(image_path_1) else None
    pdb_url_1 = f'/static/jobs/{fbpp_id}/{fbpp_id}.pdb' if os.path.exists(pdb_path_1) else None
    print(pdb_url_1)
    return render(request, 'browser/results.html', {
        'heatmap_html': heatmap_html,
        'query': fbpp_id,
        'file_path': full_path,
        'image_url_1': image_url_1,
        'pdb_url_1': pdb_url_1,
    })


def download_file(request, fbpp_id):
    file_path = f"jobs/{fbpp_id}/{fbpp_id}_normPred_evolCombi.txt"
    if not os.path.exists(file_path):
        return HttpResponse("File not found.")
    
    response = FileResponse(open(file_path, 'rb'))
    response['Content-Disposition'] = f'attachment; filename="{fbpp_id}_normPred_evolCombi.txt"'
    return response
