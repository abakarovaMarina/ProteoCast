from django.urls import path
from .views import search_view, results

urlpatterns = [
    path('', search_view, name='search'),
    path('search/', search_view, name='search'),
    path('results/', results, name='results'),
]