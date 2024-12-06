from django.urls import path
from .views import search_view, results_view, download_folder

urlpatterns = [
    path('results/download/<str:fbpp_id>/', download_folder, name='download_file'),
    path('search/', search_view, name='search'),
    path('results/', results_view, name='results'),
]