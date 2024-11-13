# browser/views.py
from django.shortcuts import render

def search_view(request):
    return render(request, 'browser/search.html')
