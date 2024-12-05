# browser/views.py
from django.shortcuts import render

def search_view(request):
    return render(request, 'browser/search.html')

def contact_us(request):
    return render(request, 'browser/contact_us.html')

def documentation(request):
    return render(request, 'browser/documentation.html')

def citation(request):
    return render(request, 'browser/citation.html')
