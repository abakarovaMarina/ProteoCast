from django.shortcuts import render
# nombre_de_la_aplicacion/views.py
from django.http import HttpResponse

def search_view(request):
    return render(request, 'browser/search.html')