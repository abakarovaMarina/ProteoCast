from django.shortcuts import render
from django.http import HttpResponse

def search_view(request):
    return render(request, 'browser/search.html')

def results_view(request):
    query = request.GET.get('q') 
    
    results = [
        {"title": "Resultado 1", "description": "Description 1."},
        {"title": "Resultado 2", "description": "Description 2."},
        {"title": "Resultado 3", "description": "Description 3."},
    ]

    if query:
        filtered_results = [result for result in results if query.lower() in result["title"].lower()]
    else:
        filtered_results = results 

    context = {
        "query": query,
        "results": filtered_results,
    }
    
    return render(request, 'browser/results.html', context)
