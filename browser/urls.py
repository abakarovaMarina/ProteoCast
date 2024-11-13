from django.urls import path
from .views import search_view, results_view, download_file

urlpatterns = [
    path('results/download/<str:fbpp_id>/', download_file, name='download_file'),
    path('search/', search_view, name='search'),
    path('results/', results_view, name='results'),
]