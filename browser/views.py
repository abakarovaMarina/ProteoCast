import os
import pandas as pd
import plotly.express as px
from django.shortcuts import render
from django.http import HttpResponse, FileResponse


dir ="FBpp_example/"  

def search_view(request):
    return render(request, 'browser/search.html')


def results(request):
    fbpp_id = request.GET.get('q')
    if not fbpp_id:
        return HttpResponse("Please provide a FBpp ID.")

    file_path = f"{dir}{fbpp_id}/{fbpp_id}_normPred_evolCombi.txt"
    if not os.path.exists(file_path):
        return HttpResponse("File not found.")

    # Read the file into a DataFrame
    df = pd.read_csv(file_path, skiprows=1, sep=' ', header=None, index_col=0)
    reversed_oranges = px.colors.sequential.Oranges[::-1]
    # Generate the heatmap using Plotly
    fig = px.imshow(df, color_continuous_scale=reversed_oranges, title=f"GEMME mutational landscape for {fbpp_id}")
    fig.update_layout(title_x=0.5)
    fig.update_layout(
        autosize=False,
        width=1500,  # Adjust the width as needed
        height=500,  # Adjust the height as needed
    )
    # Convert the plotly figure to HTML
    fig.update_yaxes(tickfont=dict(size=14), dtick=1)  # Increase the yticks size and ensure all ticks are visible
    heatmap_html = fig.to_html(full_html=False)

    return render(request, 'browser/results.html', {'heatmap_html': heatmap_html, 'query': fbpp_id, 'file_path': file_path})

def download_file(request, fbpp_id):
    file_path = f"{dir}{fbpp_id}/{fbpp_id}_normPred_evolCombi.txt"
    if not os.path.exists(file_path):
        return HttpResponse("File not found.")
    
    response = FileResponse(open(file_path, 'rb'))
    response['Content-Disposition'] = f'attachment; filename="{fbpp_id}_normPred_evolCombi.txt"'
    return response