import os
import pandas as pd
import plotly.express as px
from django.shortcuts import render
from django.http import HttpResponse


def search_view(request):
    return render(request, 'browser/search.html')


def results(request):
    fbpp_id = request.GET.get('q')
    if not fbpp_id:
        return HttpResponse("Please provide a FBpp ID.")

    file_path = f"/Users/manchuta/Documents/GitHub/Droso_GEMMEwebsite/FBpp_example/{fbpp_id}/{fbpp_id}_normPred_evolCombi.txt"
    if not os.path.exists(file_path):
        return HttpResponse("File not found.")

    # Read the file into a DataFrame
    df = pd.read_csv(file_path, skiprows=1, sep=' ', header=None, index_col=0)

    # Generate the heatmap using Plotly
    fig = px.imshow(df, color_continuous_scale='Oranges', title=f"Heatmap for {fbpp_id}")
    fig.update_layout(
        autosize=False,
        width=800,  # Adjust the width as needed
        height=400,  # Adjust the height as needed
        margin=dict(l=50, r=50, b=100, t=100, pad=4)
    )
    # Convert the plotly figure to HTML
    heatmap_html = fig.to_html(full_html=False)

    return render(request, 'browser/results.html', {'heatmap_html': heatmap_html, 'query': fbpp_id})
